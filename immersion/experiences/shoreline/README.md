# Shoreline

The walkable development scene: an authored cove with clipmap terrain out to about
16 km, 173 procedurally built trees, rocks and ground cover placed by habitat
rules, a sea with lunar-gravity waves, rain and surface water, the clear-sky
atlas, and clouds drawn from atmospheric column soundings. It exists to develop
and test the engine. Its landscape and ecology are placeholders, not proposals
for an Open Moon shore.

Run it from `immersion/` with `npm run dev` and open `/experiences/shoreline/`.
Drag to look, WASD to walk, Space to hop (lunar gravity), H to hide the controls.
**Show controls** opens viewpoints, the Sun slider, weather, and the column
studies.

`main.js` loads the baked sky atlas and surface textures and boots `app.js`,
which holds the render loop, controls and panels. For browser harnesses the page
exposes `window.openMoonShoreline` (state, `setPhase`, `setLocation`,
`setWeather`, `setColumn`, `renderOnce`, `snapshotState`, `benchmark`). The
measurement probes in [illumination/references](../../../illumination/references/)
also read `THREE`, `OM` and `OpenMoonCloudRenderer` from the page.
