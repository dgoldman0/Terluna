'use strict';
/* Export the column-set product: every preset sounding for each world, with its
 * diagnosed rows, summary and the renderer texture sampling.
 *
 * Usage: node atmosphere/column/export_columns.cjs [output.json]
 * Consumers read rows by linear interpolation in z (sampleRows); they do not run
 * the column model. Values are JSON doubles, so the product round-trips exactly.
 */
const fs = require('node:fs'),
  path = require('node:path'),
  crypto = require('node:crypto');
const W = require('./weather-column.js');

const WORLDS = ['moon', 'moon_no_ozone', 'earth'];
const sha256 = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');

function finite(value, where) {
  if (typeof value === 'number' && !Number.isFinite(value))
    throw new Error(`Non-finite value at ${where}; JSON would not round-trip it.`);
  if (value && typeof value === 'object')
    for (const [k, v] of Object.entries(value)) finite(v, `${where}.${k}`);
}

const columns = [];
for (const world of WORLDS)
  for (const key of Object.keys(W.PRESETS)) {
    const model = W.create(key, world),
      texture = W.textureData(model);
    const column = {
      schema: model.schema,
      world,
      key: model.key,
      planet: model.planet,
      inputs: model.inputs,
      summary: model.summary,
      rows: model.rows,
      texture: { base: texture.base, top: texture.top, n: texture.n, data: Array.from(texture.data) },
    };
    finite(column, `${world}.${key}`);
    columns.push(column);
  }

const product = {
  schema: 'terluna.atmosphere.column-set/1',
  producer: {
    domain: 'atmosphere',
    model: 'atmosphere/column/weather-column.js',
    model_sha256: sha256(path.join(__dirname, 'weather-column.js')),
    exporter: 'atmosphere/column/export_columns.cjs',
  },
  evidence:
    'Selected column experiments from prescribed soundings; not forecasts or climatology. ' +
    'Retained condensate is an optical closure, not a precipitation budget.',
  units: 'SI: z in m, pressure in Pa, temperature in K, extinction in 1/m, wind in m/s',
  interpolation: 'Linear in z between adjacent rows (sampleRows in weather-column.js).',
  presets: W.PRESETS,
  columns,
};

const out = path.resolve(process.argv[2] || path.join(__dirname, 'products', 'columns.json'));
fs.mkdirSync(path.dirname(out), { recursive: true });
fs.writeFileSync(out, JSON.stringify(product));
console.log(`${columns.length} columns -> ${out} (${fs.statSync(out).size.toLocaleString()} bytes)`);
