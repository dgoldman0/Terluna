/* Atmospheric columns for the experience, read from the baked column set.
 *
 * The atmosphere domain computes the columns (atmosphere/column); bake/columns.mjs
 * packs its product into assets/columns/columns.json. This module only reads that
 * data: rows are interpolated linearly in z, the product's stated rule.
 */
export const columns = { presets: null, worlds: null, aliases: {}, source: null };

export function loadColumns(asset) {
  if (asset?.schema !== 'terluna.immersion.columns/1')
    throw new Error('Unexpected column asset schema: ' + asset?.schema);
  Object.assign(columns, {
    presets: asset.presets,
    worlds: asset.worlds,
    aliases: asset.aliases ?? {},
    source: asset.source,
  });
}

/** Linear interpolation of every row field at height z (metres). */
export function sampleRows(rows, z) {
  let lo = 0,
    hi = rows.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (rows[m].z <= z) lo = m;
    else hi = m;
  }
  const a = rows[lo],
    b = rows[hi],
    t = (z - a.z) / (b.z - a.z),
    out = {};
  for (const k of Object.keys(a)) out[k] = a[k] + (b[k] - a[k]) * t;
  out.z = z;
  return out;
}

const views = new Map();
/** One column by world ('moon', 'moon_no_ozone', 'earth') and preset key. */
export function getColumn(world, key) {
  if (!columns.worlds) throw new Error('Column asset not loaded.');
  const id = world + ':' + key;
  if (!views.has(id)) {
    const data = columns.worlds[columns.aliases[world] ?? world]?.[key];
    if (!data) throw new RangeError(`No baked column for ${world}/${key}`);
    const texture = { ...data.texture, data: Float32Array.from(data.texture.data) };
    views.set(id, { ...data, world, texture, sample: (z) => sampleRows(data.rows, z) });
  }
  return views.get(id);
}
