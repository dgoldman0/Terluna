/* Shoreline experience entry point.
 *
 * Loads the baked assets this experience reads (the clear-sky atlas packed from
 * illumination/sky, and the procedural surface textures from bake/surfaces.py),
 * then boots the app. Engine modules register their systems on the shared OM
 * registry when imported; the order below matches the original script order.
 */
import * as THREE from 'three';
import { OM } from '../../engine/om.js';
import '../../engine/atmosphere.js';
import '../../engine/materials.js';
import '../../engine/terrain.js';
import '../../engine/surface-water.js';
import '../../engine/ponds.js';
import '../../engine/ecology.js';
import '../../engine/scene.js';
import '../../engine/water.js';
import '../../engine/audio.js';
import './app.js';
import CloudRenderer from '../../engine/cloud-renderer.js';
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
    status.textContent = 'Loading the sky atlas…';
    const response = await fetch(SKY_ATLAS);
    if (!response.ok) throw Error(`Sky atlas unavailable (${response.status}).`);
    OM.skyData = await response.json();
    OM.surfaceSpec = { manifest: surfaceManifest, albedo: ALBEDO.href, packed: PACKED.href };
    // Measurement probes (illumination/references) read these from the page, as
    // they did when every module was a global script.
    Object.assign(globalThis, { THREE, OM, OpenMoonCloudRenderer: CloudRenderer });
    await OM.boot(THREE);
  } catch (e) {
    console.error(e);
    globalThis.openMoonShorelineError = e.message;
    status.textContent = e.message;
    document.getElementById('boot-actions')?.removeAttribute('hidden');
  }
}

document.getElementById('retry')?.addEventListener('click', () => location.reload());
document
  .getElementById('boot-about')
  ?.addEventListener('click', () => document.getElementById('help').showModal());
start();
