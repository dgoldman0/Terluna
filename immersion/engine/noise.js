/* Seeded value noise shared by terrain authoring and vegetation geometry. */
const mix = (a, b, t) => a + (b - a) * t;
export function hash(x, z, s = 0) {
  let n = Math.imul(x | 0, 374761393) + Math.imul(z | 0, 668265263) + Math.imul(s | 0, 1442695041);
  n = Math.imul(n ^ (n >>> 13), 1274126177);
  return ((n ^ (n >>> 16)) >>> 0) / 4294967295;
}
export function noise(x, z, s = 0) {
  const i = Math.floor(x),
    j = Math.floor(z);
  let u = x - i,
    v = z - j;
  u = u * u * (3 - 2 * u);
  v = v * v * (3 - 2 * v);
  return mix(
    mix(hash(i, j, s), hash(i + 1, j, s), u),
    mix(hash(i, j + 1, s), hash(i + 1, j + 1, s), u),
    v,
  );
}
export function fbm(x, z, s = 0) {
  return (
    0.58 * noise(x, z, s) +
    0.28 * noise(x * 2.03 + 17.4, z * 2.03 + 8.1, s) +
    0.14 * noise(x * 4.13 + 31.7, z * 4.13 + 19.6, s)
  );
}
