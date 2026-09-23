/* Coordinated, world-scaled PBR channels and shared landscape/surface state. */
import { OM } from './om.js';
import L from '../world/landscape.js';
OM.loadSurfaceAssets = async function (T) {
  const spec = OM.surfaceSpec,
    loader = new T.TextureLoader();
  const load = (uri) =>
    new Promise((resolve, reject) => loader.load(uri, resolve, undefined, reject));
  const [albedo, packed] = await Promise.all([load(spec.albedo), load(spec.packed)]);
  // Each channel is uploaded as an independently mipmapped, repeating layer.
  // Array layers prevent distant atlas mipmaps from blending adjacent materials.
  const makeArray = (image, colourSpace) => {
    const { tile_pixels: size, padding_pixels: pad, layout, layers } = spec.manifest;
    const canvas = document.createElement('canvas');
    canvas.width = image.width;
    canvas.height = image.height;
    const ctx = canvas.getContext('2d', { willReadFrequently: true });
    ctx.drawImage(image, 0, 0);
    const bytes = new Uint8Array(size * size * 4 * layers.length);
    for (const layer of layers) {
      const x = (layer.layer % layout[0]) * (size + 2 * pad) + pad,
        y = Math.floor(layer.layer / layout[0]) * (size + 2 * pad) + pad;
      bytes.set(ctx.getImageData(x, y, size, size).data, layer.layer * size * size * 4);
    }
    const t = new T.DataArrayTexture(bytes, size, size, layers.length);
    t.colorSpace = colourSpace;
    t.wrapS = t.wrapT = T.RepeatWrapping;
    t.generateMipmaps = true;
    t.minFilter = T.LinearMipmapLinearFilter;
    t.magFilter = T.LinearFilter;
    t.anisotropy = 8;
    t.unpackAlignment = 1;
    t.needsUpdate = true;
    return t;
  };
  OM.surfaceAssets = {
    albedo: makeArray(albedo.image, T.SRGBColorSpace),
    packed: makeArray(packed.image, T.NoColorSpace),
    manifest: spec.manifest,
  };
  albedo.dispose();
  packed.dispose();
};
OM.createFieldUniforms = function (T, water) {
  const grid = L.GRID;
  const texture = (bytes, n) => {
    const t = new T.DataTexture(bytes, n, n, T.RGBAFormat, T.FloatType);
    t.minFilter = t.magFilter = T.LinearFilter;
    t.generateMipmaps = false;
    t.needsUpdate = true;
    return t;
  };
  const bathymetry = new Float32Array(grid.n * grid.n * 4);
  for (let i = 0; i < grid.n * grid.n; i++) bathymetry[i * 4] = L.data.elevation[i];
  const uniforms = {
    uBathymetry: { value: texture(bathymetry, grid.n) },
    uMaterialFields: { value: texture(L.materialFieldRGBA(), grid.n) },
    uEnvironmentFields: { value: texture(L.environmentFieldRGBA(), grid.n) },
    uSurfaceState: { value: texture(water.textureData(), water.n) },
    uFieldBounds: { value: new T.Vector4(grid.minX, grid.minZ, grid.cell, grid.n) },
    uStateBounds: { value: new T.Vector4(water.minX, water.minZ, water.cell, water.n) },
    uSurfaceAlbedo: { value: OM.surfaceAssets.albedo },
    uSurfacePacked: { value: OM.surfaceAssets.packed },
    uPondLevels: { value: texture(new Float32Array(4), 1) },
    uPondBounds: { value: new T.Vector4(0, 0, 1, 1) },
    uLandscapeDebug: { value: 0 },
  };
  return uniforms;
};
OM.FIELD_GLSL = `
uniform sampler2D uMaterialFields,uEnvironmentFields,uSurfaceState,uPondLevels;
uniform highp sampler2DArray uSurfaceAlbedo,uSurfacePacked;
uniform vec4 uFieldBounds,uStateBounds,uPondBounds;uniform float uLandscapeDebug;
vec2 omGridUV(vec2 p,vec4 bounds){return ((p-bounds.xy)/bounds.z+.5)/bounds.w;}
float omGridInside(vec2 p,vec4 b){vec2 q=(p-b.xy)/b.z;vec2 a=smoothstep(vec2(0),vec2(6),q),c=1.-smoothstep(vec2(b.w-7.),vec2(b.w-1.),q);return a.x*a.y*c.x*c.y;}
vec4 omEnvironment(vec2 p){return texture2D(uEnvironmentFields,omGridUV(p,uFieldBounds));}
vec4 omSurfaceStateAt(vec2 p){return texture2D(uSurfaceState,omGridUV(p,uStateBounds))*omGridInside(p,uStateBounds);}
float omNoise2(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.-2.*f);vec4 v=fract(sin(vec4(dot(i,vec2(127.1,311.7)),dot(i+vec2(1,0),vec2(127.1,311.7)),dot(i+vec2(0,1),vec2(127.1,311.7)),dot(i+1.,vec2(127.1,311.7))))*43758.5453);return mix(mix(v.x,v.y,f.x),mix(v.z,v.w,f.x),f.y);}
// Unwrapped UV derivatives and independently repeating array layers.
vec4 omAtlas(highp sampler2DArray tex,vec2 uv,float layer){return textureGrad(tex,vec3(uv,layer),dFdx(uv),dFdy(uv));}
struct OMSurface{vec3 albedo;vec2 slope;float rough;float height;float ao;};
OMSurface omLayer(vec2 uv,float layer){vec4 a=omAtlas(uSurfaceAlbedo,uv,layer),p=omAtlas(uSurfacePacked,uv,layer);OMSurface m;m.albedo=a.rgb;m.slope=(p.rg*2.-1.)/max(.20,sqrt(max(.02,1.-dot(p.rg*2.-1.,p.rg*2.-1.))));m.rough=p.b;m.height=p.a-.5;m.ao=a.a;return m;}
OMSurface omBlend(OMSurface a,OMSurface b,float t){OMSurface m;m.albedo=mix(a.albedo,b.albedo,t);m.slope=mix(a.slope,b.slope,t);m.rough=mix(a.rough,b.rough,t);m.height=mix(a.height,b.height,t);m.ao=mix(a.ao,b.ao,t);return m;}
vec3 omMappedNormal(vec3 normal,vec3 gradient){vec3 g=mat3(viewMatrix)*gradient;return normalize(normal+g-normal*dot(normal,g));}
`;
OM.material = function (T, atm, state, kind, options = {}) {
  const leaf = kind === 'leaf' || kind === 'grass',
    mat = new T.MeshStandardMaterial({ color: 0xffffff, roughness: 0.9, ...options });
  OM.prepareMaterial(T, mat, atm, state, {
    ground: kind === 'terrain',
    sway: kind === 'leaf' ? 0.015 : kind === 'grass' ? 0.025 : 0,
  });
  const previous = mat.onBeforeCompile;
  mat.onBeforeCompile = (shader) => {
    previous(shader);
    shader.fragmentShader = shader.fragmentShader.replace(
      'mix(1.,.58,vOMCanopy)',
      'mix(1.,.68,habitat.g)',
    );
    Object.assign(shader.uniforms, state.fieldUniforms);
    shader.vertexShader =
      'varying vec2 vOMTex;varying vec2 vOMSize;' +
      (kind === 'terrain' ? 'attribute vec4 omWeights;varying vec4 vOMWeights;' : '') +
      '\n' +
      shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace(
      '#include <begin_vertex>',
      `#include <begin_vertex>\nvOMTex=uv;vOMSize=vec2(1);\n#ifdef USE_INSTANCING\nvOMSize=vec2(length(instanceMatrix[0].xyz),length(instanceMatrix[1].xyz));\n#endif\n${kind === 'terrain' ? 'vOMWeights=omWeights;' : ''}`,
    );
    shader.fragmentShader =
      'varying vec2 vOMTex;varying vec2 vOMSize;' +
      (kind === 'terrain' ? 'varying vec4 vOMWeights;' : '') +
      '\n' +
      OM.FIELD_GLSL +
      shader.fragmentShader;
    let detail = '';
    if (kind === 'terrain')
      detail = `
   vec2 p=vOMWorld.xz;float elevation=vOMWorld.y+dot(p,p)/(2.*uR);
   vec4 habitat=omEnvironment(p),surfaceState=omSurfaceStateAt(p);
   vec4 w=mix(vOMWeights,texture2D(uMaterialFields,omGridUV(p,uFieldBounds)),omGridInside(p,uFieldBounds));w=max(w,vec4(0));w/=max(.001,dot(w,vec4(1)));
   OMSurface sand=omLayer(p/.75,0.),gravel=omLayer(p/1.5,1.),rock=omLayer(p/2.8,2.),soil=omLayer(p/1.1,3.);
   float litter=habitat.g*omGridInside(p,uFieldBounds)*smoothstep(.08,.42,w.a);
   soil=omBlend(soil,omLayer(p/1.1,4.),litter);
   // Mean unresolved herb cover blends in beyond the explicit foreground tufts.
   float coverDistance=smoothstep(28.,85.,length(p-cameraPosition.xz));
   float herbCover=coverDistance*(1.-litter)*smoothstep(1.3,4.,elevation)*(.45+.38*omNoise2(p*.05));
   soil.albedo=mix(soil.albedo,vec3(.102,.153,.035)*(.78+.38*omNoise2(p*.13)),herbCover);
   vec3 albedo=sand.albedo*w.r+gravel.albedo*w.g+rock.albedo*w.b+soil.albedo*w.a;
   vec2 microSlope=sand.slope*w.r+gravel.slope*w.g+rock.slope*w.b+soil.slope*w.a;
   float microRough=sand.rough*w.r+gravel.rough*w.g+rock.rough*w.b+soil.rough*w.a;
   float microAO=sand.ao*w.r+gravel.ao*w.g+rock.ao*w.b+soil.ao*w.a;
   float broad=omNoise2(p*.035);albedo*=.90+.20*broad;
   float dampContact=1.-smoothstep(-.03,.24,elevation);
   vec4 pondHead=texture2D(uPondLevels,omGridUV(p,uPondBounds));
   float basinWet=pondHead.a*omGridInside(p,uPondBounds)*(1.-smoothstep(-.001,.025,elevation-pondHead.r));
   float spatialWet=max(max(surfaceState.r,basinWet),dampContact),soilWet=surfaceState.g;
   albedo*=mix(1.,.66,soilWet*(w.a+w.r*.65));
   diffuseColor.rgb*=albedo;vec3 microGradient=vec3(microSlope.x,0,microSlope.y)*.72*(1.-smoothstep(.035,.25,max(length(dFdx(p)),length(dFdy(p)))));
  `;
    else if (kind === 'rock')
      detail = `
   vec3 p=vOMWorld;vec3 wn=abs(normalize(cross(dFdx(p),dFdy(p))));wn=pow(wn,vec3(4));wn/=max(.001,wn.x+wn.y+wn.z);
   OMSurface top=omLayer(p.xz/2.8,2.),side=omLayer(p.zy/2.8,2.),front=omLayer(p.xy/2.8,2.);
   vec3 albedo=top.albedo*wn.y+side.albedo*wn.x+front.albedo*wn.z;
   float lichen=smoothstep(.64,.83,omNoise2(p.xz*.6))*smoothstep(.35,.9,wn.y)*.16;albedo=mix(albedo,vec3(.19,.21,.09),lichen);
   diffuseColor.rgb*=albedo;float microRough=top.rough*wn.y+side.rough*wn.x+front.rough*wn.z,microAO=top.ao*wn.y+side.ao*wn.x+front.ao*wn.z;
   vec3 microGradient=(vec3(top.slope.x,0,top.slope.y)*wn.y+vec3(0,side.slope.y,side.slope.x)*wn.x+vec3(front.slope.x,front.slope.y,0)*wn.z)*.7;
   vec4 surfaceState=omSurfaceStateAt(p.xz);float elevation=p.y+dot(p.xz,p.xz)/(2.*uR),spatialWet=max(surfaceState.r,1.-smoothstep(-.03,.16,elevation));
  `;
    else if (kind === 'wood')
      detail = `
   vec2 barkUV=vec2(vOMTex.x*6.283185*vOMSize.x,vOMTex.y*vOMSize.y)/1.1;
   OMSurface bark=omLayer(barkUV,5.);diffuseColor.rgb*=bark.albedo;
   float microRough=bark.rough,microAO=bark.ao;vec3 microGradient=vec3(bark.slope.x*.25,0,bark.slope.y*.10);
   vec4 surfaceState=omSurfaceStateAt(vOMWorld.xz);float spatialWet=surfaceState.r*.7;
  `;
    else if (kind === 'timber')
      detail = `
   float stripe=omNoise2(vec2(vOMWorld.x*12.+vOMWorld.z*7.,vOMWorld.y*.65));diffuseColor.rgb*=mix(vec3(.085,.055,.029),vec3(.20,.135,.071),stripe);
   float microRough=.85,microAO=1.;vec3 microGradient=vec3(0);vec4 surfaceState=omSurfaceStateAt(vOMWorld.xz);float spatialWet=surfaceState.r;
  `;
    else
      detail = `
   float midrib=exp(-abs(vOMTex.x-.5)*85.);float vein=pow(.5+.5*cos((vOMTex.y*13.+abs(vOMTex.x-.5)*8.)*6.283185),18.)*smoothstep(.015,.18,abs(vOMTex.x-.5));
   diffuseColor.rgb*=mix(.83,1.05,vOMTex.y)*(1.-.11*midrib-.08*vein);
   float microRough=.76,microAO=1.;vec3 microGradient=vec3(0);vec4 surfaceState=omSurfaceStateAt(vOMWorld.xz);float spatialWet=surfaceState.a;
  `;
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <color_fragment>',
      '#include <color_fragment>\n' + detail,
    );
    // Replace the legacy representative-category wetness at its existing hook.
    const oldWet = /float omWet=mix\([\s\S]*?diffuseColor\.rgb\*=mix\(1\.,\.70,omWet\);/;
    shader.fragmentShader = shader.fragmentShader.replace(
      oldWet,
      `float omWet=spatialWet*(1.-max(omRoof,uProtected));roughnessFactor=mix(microRough,max(.055,microRough*.20),omWet);diffuseColor.rgb*=mix(1.,.72,omWet);`,
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <normal_fragment_maps>',
      '#include <normal_fragment_maps>\nnormal=omMappedNormal(normal,microGradient);',
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <lights_physical_fragment>',
      '#include <lights_physical_fragment>\nmaterial.specularColor=mix(material.specularColor,vec3(.02037),omWet);',
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <aomap_fragment>',
      '#include <aomap_fragment>\nreflectedLight.indirectDiffuse*=microAO;reflectedLight.indirectSpecular*=mix(1.,microAO,.4);',
    );
    if (leaf) {
      // Shadowed direct radiance is already in directLight.color at RE_Direct.
      const needle =
        'RE_Direct( directLight, geometryPosition, geometryNormal, geometryViewDir, geometryClearcoatNormal, material, reflectedLight );';
      shader.fragmentShader = shader.fragmentShader.replaceAll(
        needle,
        needle +
          `\nreflectedLight.directDiffuse+=directLight.color*diffuseColor.rgb*(.22/3.14159265)*max(0.,dot(-geometryNormal,directLight.direction));`,
      );
    }
    if (kind === 'terrain') {
      shader.fragmentShader = shader.fragmentShader.replace(
        '#include <opaque_fragment>',
        `#include <opaque_fragment>
    if(uLandscapeDebug>.5){vec3 debug=vec3(0);if(uLandscapeDebug<1.5)debug=mix(vec3(.06,.20,.36),vec3(.8,.76,.5),smoothstep(-2.,12.,elevation));else if(uLandscapeDebug<2.5)debug=w.r*vec3(.83,.7,.4)+w.g*vec3(.46,.34,.24)+w.b*vec3(.36,.43,.53)+w.a*vec3(.18,.4,.17);else if(uLandscapeDebug<3.5)debug=mix(vec3(.15),vec3(.1,.65,1.),habitat.b);else if(uLandscapeDebug<4.5)debug=vec3(surfaceState.g,surfaceState.r,smoothstep(0.,4.,surfaceState.b));else debug=vec3(habitat.g,habitat.r,1.-habitat.g);gl_FragColor.rgb=debug;}
   `,
      );
    }
  };
  mat.customProgramCacheKey = () => `open-moon-landscape-2-${kind}`;
  return mat;
};
