/* Renderer lab entry point: the shoreline scene through Three.js WebGPU/TSL.
 *
 * An engine experiment, kept optional until it shows value over the WebGL path
 * (see docs/renderer-roadmap.md). It reuses the shoreline engine modules and
 * swaps material factories through the shared OM registry.
 */
import * as webgl from 'three';
import * as webgpu from 'three/webgpu';
import * as tsl from 'three/tsl';
import { OM } from '../../engine/om.js';
import '../../engine/atmosphere.js';
import '../../engine/materials.js';
import '../../engine/terrain.js';
import '../../engine/surface-water.js';
import '../../engine/ponds.js';
import '../../engine/ecology.js';
import '../../engine/scene.js';
import '../../engine/water.js';
import './node-materials.js';
import './renderer-lab.js';
import surfaceManifest from '../../assets/surfaces/manifest.json';

const SKY_ATLAS = new URL('../../assets/sky/atmosphere.json', import.meta.url);
const ALBEDO = new URL('../../assets/surfaces/albedo-ao.png', import.meta.url);
const PACKED = new URL('../../assets/surfaces/normal-roughness-height.png', import.meta.url);

OM.downloadBlob = function (blob, name) {
  const url = URL.createObjectURL(blob),
    a = document.createElement('a');
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 30000);
};

async function start() {
  const status = document.getElementById('boot-status');
  try {
    if (status) status.textContent = 'Loading the sky atlas…';
    const response = await fetch(SKY_ATLAS);
    if (!response.ok) throw Error(`Sky atlas unavailable (${response.status}).`);
    OM.skyData = await response.json();
    OM.surfaceSpec = { manifest: surfaceManifest, albedo: ALBEDO.href, packed: PACKED.href };
    const manifest = { name: 'three', version: '0.186.0', revision: webgl.REVISION };
    await OM.bootLab({ webgl, webgpu, tsl, manifest });
  } catch (e) {
    console.error(e);
    globalThis.openMoonShorelineError = e.message;
    if (status) status.textContent = e.message;
    document.getElementById('boot-actions')?.removeAttribute('hidden');
  }
}

start();
