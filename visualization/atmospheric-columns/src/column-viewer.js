/* The atmospheric column viewer: for each preset sounding, the Open Moon column
 * beside the Earth column, with temperature, humidity and condensate profiles, a
 * readout at any altitude and the model's diagnostics. Profiles are drawn from every
 * row of the product; the readout samples rows by the product's own interpolation
 * rule (sampleAt, inlined before this script by the build). */
const DATA = JSON.parse(document.getElementById('column-data').textContent);
const NAMES = { moon: 'Open Moon', moon_no_ozone: 'Open Moon, zero ozone', earth: 'Earth' };
const WORLDS = Object.keys(DATA.worlds).sort((a, b) => (a === 'earth') - (b === 'earth'));
const state = { preset: Object.keys(DATA.presets)[0], range: 'cloud', hover: {} };
const $ = (id) => document.getElementById(id);
const km = (m) => (m === null || m === undefined ? '—' : (m / 1000).toFixed(m < 10000 ? 2 : 1));
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

const CHARTS = [
  {
    id: 'temperature',
    title: 'Temperature',
    unit: '°C',
    lines: [
      { field: 'T', label: 'Air', colour: '--air', map: (v) => v - 273.15 },
      { field: 'parcelT', label: 'Lifted parcel', colour: '--parcel', map: (v) => v - 273.15 },
    ],
  },
  {
    id: 'humidity',
    title: 'Relative humidity',
    unit: '%',
    fixed: [0, 105],
    lines: [{ field: 'rh', label: 'Relative humidity', colour: '--air', map: (v) => v * 100 }],
  },
  {
    id: 'condensate',
    title: 'Retained condensate',
    unit: 'g/m³',
    floor: 0,
    lines: [
      { field: 'liquid', label: 'Liquid', colour: '--liquid', map: (v) => v * 1000 },
      { field: 'ice', label: 'Ice', colour: '--ice', map: (v) => v * 1000 },
    ],
  },
];

function column(world) {
  return DATA.worlds[world][state.preset];
}

function altitudeTop(c) {
  const top = c.series.z.at(-1);
  if (state.range === 'column') return top;
  if (state.range === 'low') return Math.min(top, 2000);
  return Math.min(top, Math.max(1000, (c.summary.cloudTop_m ?? 2000) * 1.15));
}

function draw(canvas, c, spec, hoverZ) {
  const dpr = window.devicePixelRatio || 1,
    w = canvas.clientWidth,
    h = canvas.clientHeight;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  const g = canvas.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  const L = 46,
    R = 12,
    T = 26,
    B = 34,
    zTop = altitudeTop(c),
    zs = c.series.z;
  let lo = Infinity,
    hi = -Infinity;
  for (const line of spec.lines)
    c.series[line.field].forEach((v, i) => {
      if (zs[i] > zTop) return;
      const x = line.map(v);
      lo = Math.min(lo, x);
      hi = Math.max(hi, x);
    });
  if (spec.fixed) [lo, hi] = spec.fixed;
  if (!(hi > lo)) hi = lo + 1;
  const pad = spec.fixed ? 0 : (hi - lo) * 0.06;
  lo = spec.floor ?? lo - pad;
  hi += pad;
  const X = (v) => L + ((v - lo) / (hi - lo)) * (w - L - R),
    Y = (z) => h - B - (z / zTop) * (h - T - B);
  g.font = '11px system-ui, sans-serif';
  const { cloudBase_m: base, cloudTop_m: top } = c.summary;
  if (base !== null && top !== null && base < zTop) {
    g.fillStyle = css('--band');
    g.fillRect(L, Y(Math.min(top, zTop)), w - L - R, Y(base) - Y(Math.min(top, zTop)));
  }
  g.strokeStyle = css('--grid');
  g.fillStyle = css('--muted');
  g.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const z = (zTop * i) / 4,
      y = Math.round(Y(z)) + 0.5;
    g.beginPath();
    g.moveTo(L, y);
    g.lineTo(w - R, y);
    g.stroke();
    g.fillText(km(z), 4, y + 4);
  }
  for (let i = 0; i <= 3; i++) {
    const v = lo + ((hi - lo) * i) / 3,
      text = Math.abs(hi - lo) < 3 ? v.toFixed(2) : v.toFixed(0);
    g.fillText(text, X(v) - g.measureText(text).width / 2, h - B + 15);
  }
  g.fillText('km', 4, T - 8);
  g.fillText(`${spec.title} (${spec.unit})`, L, h - 6);
  let legend = L;
  for (const line of spec.lines) {
    g.strokeStyle = css(line.colour);
    g.lineWidth = 2;
    g.beginPath();
    let started = false;
    c.series[line.field].forEach((v, i) => {
      if (zs[i] > zTop && started) return;
      const x = X(line.map(v)),
        y = Y(Math.min(zs[i], zTop));
      if (!started) g.moveTo(x, y);
      else g.lineTo(x, y);
      started = true;
    });
    g.stroke();
    if (spec.lines.length > 1) {
      g.fillStyle = css(line.colour);
      g.fillText(line.label, legend, T - 8);
      legend += g.measureText(line.label).width + 14;
    }
  }
  if (hoverZ !== undefined && hoverZ <= zTop) {
    g.strokeStyle = css('--ink');
    g.setLineDash([3, 3]);
    g.beginPath();
    g.moveTo(L, Y(hoverZ) + 0.5);
    g.lineTo(w - R, Y(hoverZ) + 0.5);
    g.stroke();
    g.setLineDash([]);
  }
  canvas.onpointermove = (e) => {
    const r = canvas.getBoundingClientRect(),
      z = ((h - B - (e.clientY - r.top)) / (h - T - B)) * zTop;
    state.hover[c.world] = Math.max(0, Math.min(zTop, z));
    renderWorld(c.world);
  };
}

const READOUT = [
  ['Altitude', (s) => km(s.z) + ' km'],
  ['Pressure', (s) => (s.p / 1000).toFixed(2) + ' kPa'],
  ['Air temperature', (s) => (s.T - 273.15).toFixed(2) + ' °C'],
  ['Parcel temperature', (s) => (s.parcelT - 273.15).toFixed(2) + ' °C'],
  ['Parcel buoyancy', (s) => s.buoyancy.toExponential(2) + ' m/s²'],
  ['Relative humidity', (s) => (s.rh * 100).toFixed(1) + ' %'],
  ['Water-vapour mixing ratio', (s) => (s.r * 1000).toFixed(3) + ' g/kg'],
  ['Liquid / ice', (s) => `${(s.liquid * 1000).toFixed(4)} / ${(s.ice * 1000).toFixed(4)} g/m³`],
  ['Cloud extinction', (s) => (s.extinction * 1000).toFixed(3) + ' /km'],
  ['Wind (x, z)', (s) => `${s.windX.toFixed(2)}, ${s.windZ.toFixed(2)} m/s`],
  ['Air density', (s) => s.rho.toFixed(4) + ' kg/m³'],
  ['Gravity', (s) => s.g.toFixed(4) + ' m/s²'],
];

function renderWorld(world) {
  const c = column(world),
    section = document.querySelector(`[data-world="${world}"]`),
    hoverZ = state.hover[world];
  for (const spec of CHARTS) draw(section.querySelector(`canvas.${spec.id}`), c, spec, hoverZ);
  const s = sampleAt(c, hoverZ ?? 0);
  section.querySelector('.readout').innerHTML =
    `<caption>At ${hoverZ === undefined ? 'the surface (point at a chart to move)' : km(s.z) + ' km'}</caption>` +
    READOUT.map(([k, f]) => `<tr><th>${k}</th><td>${f(s)}</td></tr>`).join('');
}

const DIAGNOSTICS = [
  ['Surface pressure', (q) => `${(q.surfacePressure_Pa / 1000).toFixed(1)} kPa (${(q.surfacePressure_Pa / 101325).toFixed(2)} atm)`],
  ['Surface temperature / humidity', (q) => `${(q.surfaceT_K - 273.15).toFixed(1)} °C / ${(q.surfaceRH * 100).toFixed(0)} %`],
  ['Surface scale height', (q) => km(q.surfaceScaleHeight_m) + ' km'],
  ['Dry adiabatic lapse rate', (q) => q.dryLapse_K_per_km.toFixed(2) + ' K/km'],
  ['Lifting condensation level', (q) => km(q.lcl_m) + ' km'],
  ['Cloud base / top', (q) => `${km(q.cloudBase_m)} / ${km(q.cloudTop_m)} km`],
  ['Level of free convection', (q) => km(q.lfc_m) + ' km'],
  ['Equilibrium level', (q) => km(q.equilibriumLevel_m) + ' km'],
  ['Parcel reaches cloud', (q) => (q.parcelReachesCloud ? 'yes' : 'no')],
  ['Positive buoyancy integral', (q) => q.positiveBuoyancyIntegral_J_per_kg.toFixed(1) + ' J/kg'],
  ['Negative work below the LFC', (q) => q.preLFCNegativeWork_J_per_kg.toFixed(1) + ' J/kg'],
  ['Liquid / ice water path', (q) => `${(q.liquidWaterPath_kg_m2 * 1000).toFixed(1)} / ${(q.iceWaterPath_kg_m2 * 1000).toFixed(1)} g/m²`],
  ['In-cloud optical depth', (q) => q.inCloudOpticalDepth.toFixed(2)],
  ['Hydrostatic weight error', (q) => q.hydrostaticWeightRelativeError.toExponential(1)],
  ['Pressure at column top', (q) => (q.pressureFloor_Pa / 1000).toFixed(2) + ' kPa'],
];

function renderDiagnostics() {
  const cols = WORLDS.map(column);
  $('diagnostics').innerHTML =
    `<thead><tr><th></th>${WORLDS.map((w) => `<th>${NAMES[w]}</th>`).join('')}</tr></thead><tbody>` +
    [
      ['Planet radius / surface gravity', (c) => `${(c.planet.radius / 1000).toFixed(1)} km / ${c.planet.g0} m/s²`],
      ...DIAGNOSTICS.map(([k, f]) => [k, (c) => f(c.summary)]),
    ]
      .map(([k, f]) => `<tr><th>${k}</th>${cols.map((c) => `<td>${f(c)}</td>`).join('')}</tr>`)
      .join('') +
    '</tbody>';
  const i = cols[0].inputs,
    m = i.morphology;
  $('inputs').textContent =
    `Prescribed for both worlds: surface ${(i.surfaceT - 273.15).toFixed(1)} °C at ${(i.rh * 100).toFixed(0)} % humidity, ` +
    `parcel heating ${i.parcelHeating} K, launch speed ${i.launchSpeed} m/s, condensate retention ${i.retention}, ` +
    `cloud cover ${m.coverage}, droplet / crystal radius ${m.liquidRadius_um} / ${m.iceRadius_um} µm. ` +
    'Temperature, humidity and wind are specified against ln(surface pressure / pressure).';
  $('assumptions').innerHTML = (cols[0].summary.assumptions || []).map((a) => `<li>${a}</li>`).join('');
}

function render() {
  document
    .querySelectorAll('[data-preset]')
    .forEach((b) => b.classList.toggle('selected', b.dataset.preset === state.preset));
  document
    .querySelectorAll('[data-range]')
    .forEach((b) => b.classList.toggle('selected', b.dataset.range === state.range));
  WORLDS.forEach(renderWorld);
  renderDiagnostics();
}

function download(world) {
  const c = column(world);
  const text = JSON.stringify(
    {
      schema: c.schema,
      source: DATA.source,
      world: c.world,
      preset: c.key,
      planet: c.planet,
      inputs: c.inputs,
      summary: c.summary,
      units: DATA.units,
      interpolation: DATA.interpolation,
      rows: rowsOf(c),
    },
    null,
    1,
  );
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type: 'application/json' }));
  a.download = `terluna-column-${c.world}-${c.key}.json`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function build() {
  $('presets').innerHTML = Object.keys(DATA.presets)
    .map((k) => `<button data-preset="${k}">${DATA.worlds[WORLDS[0]][k].summary.name}</button>`)
    .join('');
  $('worlds').innerHTML = WORLDS.map(
    (w) => `<section data-world="${w}">
      <h2>${NAMES[w]}</h2>
      ${Object.entries(DATA.aliases)
        .filter(([, to]) => to === w)
        .map(([from]) => `<p class="note">${NAMES[from] || from}: identical columns.</p>`)
        .join('')}
      ${CHARTS.map((c) => `<canvas class="${c.id}" aria-label="${c.title} against altitude, ${NAMES[w]}"></canvas>`).join('')}
      <table class="readout"></table>
      <button data-download="${w}">Download this column (JSON)</button>
    </section>`,
  ).join('');
  document.querySelectorAll('[data-preset]').forEach(
    (b) =>
      (b.onclick = () => {
        state.preset = b.dataset.preset;
        state.hover = {};
        render();
      }),
  );
  document.querySelectorAll('[data-range]').forEach(
    (b) =>
      (b.onclick = () => {
        state.range = b.dataset.range;
        state.hover = {};
        render();
      }),
  );
  document.querySelectorAll('[data-download]').forEach((b) => (b.onclick = () => download(b.dataset.download)));
  const p = DATA.source.producer;
  $('provenance').innerHTML =
    `<p>${DATA.evidence}</p><p>Units: ${DATA.units}. Interpolation: ${DATA.interpolation}</p>` +
    `<p>Product <code>${DATA.source.schema}</code>, sha256 <code>${DATA.source.sha256}</code>; ` +
    `model <code>${p.model}</code> (sha256 <code>${p.model_sha256}</code>), exported by <code>${p.exporter}</code>.</p>`;
  render();
  window.addEventListener('resize', () => WORLDS.forEach(renderWorld));
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener?.('change', render);
  globalThis.columnViewer = { data: DATA, state, render, sampleAt, rowsOf };
}

build();
