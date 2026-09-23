/* Bright stars in the lunar sky, from the illumination domain's catalogue (BSC5).
 *
 * Each star is a point of light at its J2000 position, turned into the scene by the
 * site's celestial rotation (engine/sky-bodies.js). Brightness follows V magnitude
 * (V = 0 gives 2.54e-6 lux above the atmosphere), colour follows B-V, and light is
 * dimmed by the atmosphere's RGB transmission along the line of sight, by cloud cover
 * and by fog. Stars share the sky's display scale, so exposure decides which show.
 */
import O from './cloud-optics.js';

const LUX_AT_V0 = 2.54e-6;
const SPRITE_PIXELS = 2;
// Stars are drawn at the eye's resolution (about one arcminute), not the screen's:
// a screen pixel spans several arcminutes and would dilute a point of light far below
// what the eye sees against the same sky. A perceptual convention, stated here.
const EYE_SOLID_ANGLE = (Math.PI / 180 / 60) ** 2;
const MU_MIN = -0.05;

/** Approximate linear-sRGB colour of a star with colour index B-V, unit luminance. */
export function starTint(bv) {
  const t = 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62));
  const k = t / 100;
  const r = k <= 66 ? 255 : 329.698727446 * Math.pow(k - 60, -0.1332047592);
  const g =
    k <= 66
      ? 99.4708025861 * Math.log(k) - 161.1195681661
      : 288.1221695283 * Math.pow(k - 60, -0.0755148492);
  const b = k >= 66 ? 255 : k <= 19 ? 0 : 138.5177312231 * Math.log(k - 10) - 305.0447927307;
  const lin = [r, g, b].map((v) => {
    const c = Math.min(255, Math.max(0, v)) / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  });
  const y = lin[0] * 0.2126 + lin[1] * 0.7152 + lin[2] * 0.0722;
  return lin.map((v) => v / Math.max(y, 1e-6));
}

export class Stars {
  constructor(T, catalogue, atmosphere) {
    if (catalogue?.schema !== 'terluna.illumination.bright-stars/1')
      throw new Error('Unexpected star catalogue schema: ' + catalogue?.schema);
    this.T = T;
    this.atmosphere = atmosphere;
    const n = catalogue.stars.length,
      direction = new Float32Array(n * 3),
      light = new Float32Array(n * 3);
    catalogue.stars.forEach(([ra, dec, v, bv], i) => {
      direction.set(
        [Math.cos(dec) * Math.cos(ra), Math.cos(dec) * Math.sin(ra), Math.sin(dec)],
        i * 3,
      );
      const flux = (LUX_AT_V0 * Math.pow(10, -0.4 * v)) / 8500; // the sky's display scale
      light.set(
        starTint(bv).map((c) => c * flux),
        i * 3,
      );
    });
    const geometry = new T.BufferGeometry();
    geometry.setAttribute('position', new T.BufferAttribute(direction, 3));
    geometry.setAttribute('light', new T.BufferAttribute(light, 3));
    this.extinction = new T.DataTexture(
      new Float32Array(128 * 4),
      128,
      1,
      T.RGBAFormat,
      T.FloatType,
    );
    this.extinction.minFilter = this.extinction.magFilter = T.LinearFilter;
    const au = atmosphere.uniforms;
    this.uniforms = {
      uCelestial: { value: new T.Matrix3() },
      uExtinction: { value: this.extinction },
      uShow: { value: 1 },
      uWhiteBalance: au.uWhiteBalance,
      uCover: au.uCover,
      uTau: au.uTau,
      uVisibility: au.uVisibility,
      uEarth: au.uEarth,
      uEarthCos: au.uEarthCos,
      uEarthShow: au.uEarthShow,
    };
    this.material = new T.ShaderMaterial({
      uniforms: this.uniforms,
      vertexShader: `
        attribute vec3 light;
        uniform mat3 uCelestial;
        uniform float uShow, uCover, uTau, uVisibility, uEarthCos, uEarthShow;
        uniform vec3 uEarth;
        uniform sampler2D uExtinction;
        varying vec3 vRadiance;
        void main() {
          vec3 d = uCelestial * position;
          float behindEarth = uEarthShow > .5 && dot(d, uEarth) > uEarthCos ? 0. : 1.;
          float mu = clamp((d.y - ${MU_MIN.toFixed(2)}) / (1. - ${MU_MIN.toFixed(2)}), 0., 1.);
          vec3 through = texture(uExtinction, vec2(mu, .5)).rgb;
          float clouds = 1. - uCover * (1. - exp(-uTau));
          float fog = uVisibility < 35000. ? exp(-3.912 / max(80., uVisibility) * (uVisibility < 500. ? 55. : 350.) / max(.03, d.y)) : 1.;
          float up = d.y > ${MU_MIN.toFixed(2)} ? 1. : 0.;
          vRadiance = light * through * clouds * fog * behindEarth * up * uShow / ${EYE_SOLID_ANGLE.toExponential(6)};
          vec4 p = projectionMatrix * mat4(mat3(viewMatrix)) * vec4(d, 1.);
          gl_Position = p.xyww;
          gl_PointSize = ${SPRITE_PIXELS.toFixed(1)};
        }`,
      fragmentShader: `
        uniform vec3 uWhiteBalance;
        varying vec3 vRadiance;
        void main() {
          gl_FragColor = vec4(vRadiance, 1.);
          #ifdef TONE_MAPPING
          gl_FragColor.rgb *= uWhiteBalance;
          #endif
          #include <tonemapping_fragment>
          #include <colorspace_fragment>
        }`,
      blending: T.AdditiveBlending,
      depthWrite: false,
      depthTest: true,
      transparent: true,
      toneMapped: true,
    });
    this.points = new T.Points(geometry, this.material);
    this.points.frustumCulled = false;
    this.points.renderOrder = -99;
    this.world = null;
    this.eyeHeight = null;
  }
  /** Line-of-sight RGB transmission from the eye to space, by elevation. */
  updateExtinction(world, eyeHeight) {
    const data = this.extinction.image.data;
    for (let i = 0; i < 128; i++) {
      const mu = MU_MIN + ((1 - MU_MIN) * i) / 127;
      const t = O.transmission(world, Math.max(0, eyeHeight), mu);
      data.set([t[0], t[1], t[2], 1], i * 4);
    }
    this.extinction.needsUpdate = true;
    this.world = world;
    this.eyeHeight = eyeHeight;
  }
  update(camera, renderer, world) {
    const earth = this.atmosphere.earth;
    this.points.visible = world !== 'earth' && !!earth;
    if (!this.points.visible) return;
    const eye = Math.round(this.atmosphere.uniforms.uEyeHeight.value);
    if (world !== this.world || eye !== this.eyeHeight) this.updateExtinction(world, eye);
    const m = earth.celestialToScene;
    this.uniforms.uCelestial.value.set(...m[0], ...m[1], ...m[2]);
    const height = renderer.getDrawingBufferSize(new this.T.Vector2()).y,
      pixel = (2 * Math.tan((camera.fov * Math.PI) / 360)) / Math.max(1, height);
    const diameter = (2 * Math.asin(this.atmosphere.uniforms.uEarthSin.value)) / pixel;
    this.atmosphere.setEarthPixels(diameter);
  }
}
