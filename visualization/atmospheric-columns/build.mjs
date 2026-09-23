/* Build the atmospheric column viewer: one self-contained page that displays the
 * atmosphere domain's column-set product.
 *
 * Usage: node build.mjs [--product columns.json] [--output page.html]
 * Without --product it reads atmosphere/column/products/columns.json and runs the
 * domain's exporter first if that file is missing. The page is deterministic for a
 * given product: it records the product's sha256 and carries no build time.
 */
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { prepare } from './src/column-data.mjs';

const here = path.dirname(fileURLToPath(import.meta.url));
const repo = path.resolve(here, '../..');
export const DEFAULT_PRODUCT = path.join(repo, 'atmosphere/column/products/columns.json');
export const DEFAULT_OUTPUT = path.join(here, 'Open_Moon_Atmospheric_Columns.html');

export function build(productPath = DEFAULT_PRODUCT, outputPath = DEFAULT_OUTPUT) {
  if (!fs.existsSync(productPath)) {
    if (path.resolve(productPath) !== DEFAULT_PRODUCT)
      throw new Error(`Column product not found: ${productPath}`);
    execFileSync(process.execPath, [path.join(repo, 'atmosphere/column/export_columns.cjs')], {
      stdio: 'inherit',
    });
  }
  const bytes = fs.readFileSync(productPath),
    sha256 = crypto.createHash('sha256').update(bytes).digest('hex'),
    data = prepare(JSON.parse(bytes), sha256);
  // The data module is inlined ahead of the viewer, so its exports become plain declarations.
  const script =
    fs.readFileSync(path.join(here, 'src/column-data.mjs'), 'utf8').replace(/^export /gm, '') +
    '\n' +
    fs.readFileSync(path.join(here, 'src/column-viewer.js'), 'utf8');
  const json = JSON.stringify(data).replace(/</g, '\\u003c');
  const page = fs
    .readFileSync(path.join(here, 'index.template.html'), 'utf8')
    .replace('__COLUMN_DATA__', () => json)
    .replace('__VIEWER_SCRIPT__', () => script);
  fs.mkdirSync(path.dirname(outputPath), { recursive: true });
  fs.writeFileSync(outputPath, page);
  return { outputPath, sha256, bytes: page.length, data };
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const arg = (name) => {
    const i = process.argv.indexOf(name);
    return i > 0 ? path.resolve(process.argv[i + 1]) : undefined;
  };
  const { outputPath, sha256, bytes, data } = build(arg('--product'), arg('--output'));
  const worlds = Object.keys(data.worlds).join(', '),
    aliases = Object.entries(data.aliases).map(([a, b]) => `${a} = ${b}`);
  console.log(
    `${outputPath} (${(bytes / 1e6).toFixed(1)} MB) from product ${sha256.slice(0, 12)}; ` +
      `worlds ${worlds}${aliases.length ? '; ' + aliases.join(', ') : ''}`,
  );
}
