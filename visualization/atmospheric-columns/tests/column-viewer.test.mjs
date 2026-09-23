/* The viewer shows the column-set product without changing it. */
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { createRequire } from 'node:module';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { prepare, rowsOf, sampleAt } from '../src/column-data.mjs';
import { build } from '../build.mjs';

const here = path.dirname(fileURLToPath(import.meta.url)),
  repo = path.resolve(here, '../../..'),
  model = createRequire(import.meta.url)(path.join(repo, 'atmosphere/column/weather-column.js'));
const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'column-viewer-')),
  productPath = path.join(dir, 'columns.json');
execFileSync(process.execPath, [path.join(repo, 'atmosphere/column/export_columns.cjs'), productPath]);
const product = JSON.parse(fs.readFileSync(productPath)),
  data = prepare(product, 'test');
const shown = (c) => data.worlds[data.aliases[c.world] ?? c.world][c.key];

test('every row, input and diagnostic is shown unchanged', () => {
  for (const c of product.columns) {
    const s = shown(c);
    assert.deepEqual(rowsOf(s), c.rows);
    assert.deepEqual(s.summary, c.summary);
    assert.deepEqual(s.inputs, c.inputs);
    assert.deepEqual(s.planet, c.planet);
    assert.equal(s.texture, undefined);
  }
});

test('a world is listed once only when its columns are identical', () => {
  assert.deepEqual(data.aliases, { moon_no_ozone: 'moon' });
  const changed = structuredClone(product);
  changed.columns.find((c) => c.world === 'moon_no_ozone').rows[5].T += 1e-9;
  assert.deepEqual(prepare(changed, 'test').aliases, {});
});

test('the readout samples by the product rule, bit for bit', () => {
  for (const c of product.columns) {
    const zs = c.rows.map((r) => r.z),
      top = zs.at(-1);
    for (const z of [-5, 0, zs[1], zs[7] * 0.3 + zs[8] * 0.7, 437.25, top / 3, top * 0.999, top + 1])
      assert.deepEqual(sampleAt(shown(c), z), model.sampleRows(c.rows, z), `${c.world} ${c.key} ${z}`);
  }
});

test('the built page embeds the prepared data and the product hash', () => {
  const { outputPath, sha256 } = build(productPath, path.join(dir, 'page.html'));
  const page = fs.readFileSync(outputPath, 'utf8'),
    json = page.split('<script type="application/json" id="column-data">')[1].split('</script>')[0];
  assert.deepEqual(JSON.parse(json), prepare(product, sha256));
  assert.equal(JSON.parse(json).source.sha256, sha256);
  assert.ok(!page.includes('__COLUMN_DATA__') && !page.includes('__VIEWER_SCRIPT__'));
  assert.ok(page.includes('function sampleAt(') && !/^export /m.test(page));
  assert.equal(build(productPath, path.join(dir, 'again.html')).bytes, page.length);
});
