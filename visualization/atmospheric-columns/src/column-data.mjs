/* Display preparation for the atmosphere domain's column-set product
 * (atmosphere/column/export_columns.cjs). Nothing here changes a value: the
 * renderer's texture samples are dropped, each column's rows become one array per
 * field, and a world whose columns equal another's is listed once, as an alias.
 * The build inlines this file into the page, so it has no imports.
 */
export const PRODUCT_SCHEMA = 'terluna.atmosphere.column-set/1';
export const DISPLAY_SCHEMA = 'terluna.visualization.atmospheric-columns/1';

export function prepare(product, sourceSha256) {
  if (product.schema !== PRODUCT_SCHEMA)
    throw new Error(`Unexpected column product schema: ${product.schema}`);
  const worlds = {};
  for (const { texture, rows, ...column } of product.columns) {
    const fields = Object.keys(rows[0]);
    const series = Object.fromEntries(fields.map((f) => [f, rows.map((r) => r[f])]));
    (worlds[column.world] ??= {})[column.key] = { ...column, fields, series };
  }
  // The same rule as the immersion bake: a later world equal to an earlier one is an alias.
  const same = (a, b) =>
    JSON.stringify(Object.values(a).map(({ world, ...c }) => c)) ===
    JSON.stringify(Object.values(b).map(({ world, ...c }) => c));
  const aliases = {};
  for (const world of Object.keys(worlds))
    for (const other of Object.keys(worlds))
      if (other < world && !aliases[other] && worlds[other] && same(worlds[world], worlds[other])) {
        aliases[world] = other;
        delete worlds[world];
        break;
      }
  return {
    schema: DISPLAY_SCHEMA,
    source: { schema: product.schema, sha256: sourceSha256, producer: product.producer },
    evidence: product.evidence,
    units: product.units,
    interpolation: product.interpolation,
    presets: product.presets,
    worlds,
    aliases,
  };
}

/** The column's rows, rebuilt from its per-field arrays. */
export function rowsOf(column) {
  const n = column.series.z.length;
  return Array.from({ length: n }, (_, i) =>
    Object.fromEntries(column.fields.map((f) => [f, column.series[f][i]])),
  );
}

/** Every field at altitude z by the product's rule: linear in z between adjacent rows,
 * held at the end rows outside the column (sampleRows in weather-column.js). */
export function sampleAt(column, z) {
  const zs = column.series.z,
    last = zs.length - 1,
    at = (i) => Object.fromEntries(column.fields.map((f) => [f, column.series[f][i]]));
  if (!Number.isFinite(z)) throw new TypeError('Altitude must be finite');
  if (z <= zs[0]) return at(0);
  if (z >= zs[last]) return at(last);
  let lo = 0,
    hi = last;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (zs[m] <= z) lo = m;
    else hi = m;
  }
  const t = (z - zs[lo]) / (zs[hi] - zs[lo]),
    out = {};
  for (const f of column.fields) {
    const a = column.series[f][lo];
    out[f] = a + (column.series[f][hi] - a) * t;
  }
  out.z = z;
  return out;
}
