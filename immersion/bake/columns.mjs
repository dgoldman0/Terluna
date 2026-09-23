/* Bake the atmosphere domain's column-set product into the experience's asset.
 *
 * Input: atmosphere/column/products/columns.json (node atmosphere/column/export_columns.cjs).
 * Output: assets/columns/columns.json, smaller but exact where the experience samples it:
 *   - texture samples (the cloud field) are kept whole;
 *   - fog rows at or below NEAR_M keep full resolution, because the fog regime samples
 *     extinction at eye height; all other rows are thinned for the sounding plot only;
 *   - a world whose columns match another's exactly (the no-ozone Moon) is aliased.
 * The script fails if thinning changes any fog-extinction sample.
 */
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { sampleRows } from '../engine/columns.js';

const here = path.dirname(new URL(import.meta.url).pathname);
const source = path.resolve(
  process.argv[2] ?? path.join(here, '../../atmosphere/column/products/columns.json'),
);
const target = path.resolve(process.argv[3] ?? path.join(here, '../assets/columns/columns.json'));
const NEAR_M = 500,
  STRIDE = 10;

const bytes = fs.readFileSync(source);
const product = JSON.parse(bytes);
if (product.schema !== 'terluna.atmosphere.column-set/1')
  throw new Error('Unexpected column product schema: ' + product.schema);

const thin = (rows, near) =>
  rows.filter((r, i) => (near && r.z <= NEAR_M) || i % STRIDE === 0 || i === rows.length - 1);
const worlds = {},
  aliases = {};
for (const column of product.columns) {
  const rows = thin(column.rows, column.key === 'fog');
  if (column.key === 'fog') {
    const top = column.rows.at(-1).z;
    for (let z = 0; z <= top; z += top / 20000) {
      if (sampleRows(rows, z).extinction !== sampleRows(column.rows, z).extinction)
        throw new Error(`Thinning changed fog extinction at z=${z} (${column.world})`);
    }
  }
  (worlds[column.world] ??= {})[column.key] = { ...column, rows };
}
const same = (a, b) =>
  JSON.stringify(Object.values(a).map(({ world, ...c }) => c)) ===
  JSON.stringify(Object.values(b).map(({ world, ...c }) => c));
for (const world of Object.keys(worlds))
  for (const other of Object.keys(worlds))
    if (other < world && !aliases[other] && same(worlds[world], worlds[other])) {
      aliases[world] = other;
      delete worlds[world];
      break;
    }

const asset = {
  schema: 'terluna.immersion.columns/1',
  source: {
    schema: product.schema,
    producer: product.producer,
    sha256: crypto.createHash('sha256').update(bytes).digest('hex'),
  },
  evidence: product.evidence,
  interpolation: product.interpolation,
  thinning: `Fog rows up to ${NEAR_M} m keep full resolution; other rows keep every ${STRIDE}th (sounding plot only). Fog sampling is exact.`,
  presets: product.presets,
  aliases,
  worlds,
};
fs.mkdirSync(path.dirname(target), { recursive: true });
fs.writeFileSync(target, JSON.stringify(asset));
console.log(
  `Baked ${Object.values(worlds).reduce((n, w) => n + Object.keys(w).length, 0)} columns ` +
    `(aliases ${JSON.stringify(aliases)}) -> ${path.relative(process.cwd(), target)} ` +
    `(${fs.statSync(target).size.toLocaleString()} bytes)`,
);
