// Run from the repository with NODE_PATH set to the installed Playwright modules.
// Optional arguments: [existing-app-url] [output-directory]. Without a URL, start the local app.
const {chromium} = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const {spawn, execFileSync} = require('node:child_process');
const root = path.resolve(__dirname, '..');
const python = path.join(root, '.venv/bin/python');
const args = process.argv.slice(2);
const suppliedURL = args[0]?.startsWith('http') ? args.shift() : null;
const out = path.resolve(args[0] || path.join(root, 'paper/output/reliability_v1/t1_rework/ui'));
fs.mkdirSync(out, {recursive: true});
let server, browser, serverLog = '';
const report = {status: 'running', checks: {}, page_errors: []};
const sha256 = value => crypto.createHash('sha256').update(value).digest('hex');
const mapFigure = /^map(?:_(?:activity|class_support|neighbor_recall|structure_group))?\.(svg|png|pdf)$/;

async function startApp() {
  server = spawn(python, ['-m', 'chemgridmap.app', '--no-browser', '--output-dir', path.join(out, 'runs')],
    {cwd: root, env: {...process.env, PYTHONPATH: path.join(root, 'src')}});
  return new Promise((resolve, reject) => {
    server.stdout.on('data', chunk => {
      serverLog += chunk.toString();
      const match = serverLog.match(/http:\/\/127\.0\.0\.1:\d+\/#[^\s]+/);
      if (match) resolve(match[0]);
    });
    server.stderr.on('data', chunk => { serverLog += chunk.toString(); });
    server.on('error', reject);
    server.on('exit', code => reject(new Error(`Local app exited (${code}): ${serverLog}`)));
  });
}

function manifest() {
  const inputs = ['src/chemgridmap/data/chembl205_ic50.csv', 'src/chemgridmap/app.py',
    'src/chemgridmap/plotting.py', 'src/chemgridmap/web/app.js', 'src/chemgridmap/web/index.html',
    'paper/check_reliability_ui.cjs'];
  const outputs = {};
  const scan = directory => {
    for (const entry of fs.readdirSync(directory, {withFileTypes: true})) {
      const file = path.join(directory, entry.name);
      if (entry.isDirectory()) scan(file);
      else if (file !== path.join(out, 'manifest.json')) outputs[path.relative(out, file)] = sha256(fs.readFileSync(file));
    }
  };
  scan(out);
  const version = fs.readFileSync(path.join(root, 'src/chemgridmap/__init__.py'), 'utf8').match(/__version__\s*=\s*["']([^"']+)/)[1];
  fs.writeFileSync(path.join(out, 'manifest.json'), JSON.stringify({software_version: version,
    parameters: {projection: 'umap', random_state: 42, browser: 'chrome', grouping_checks: ['butina', 'none']},
    code_diff_sha256: sha256(execFileSync('git', ['diff', '--binary'], {cwd: root})),
    inputs: Object.fromEntries(inputs.map(file => [file, sha256(fs.readFileSync(path.join(root, file)))])),
    outputs}, null, 2));
}

(async () => {
  const url = suppliedURL || await startApp();
  browser = await chromium.launch({channel: 'chrome', headless: true});
  const page = await browser.newPage({viewport: {width: 1440, height: 1050}, deviceScaleFactor: 2});
  page.on('pageerror', error => report.page_errors.push(error.message));
  await page.goto(url);
  await page.locator('#version').filter({hasText: '0.5.0'}).waitFor();
  await page.locator('#demo').click();
  await page.locator('#page-configure').waitFor({state: 'visible'});
  await page.locator('#projection').selectOption('umap');
  await page.locator('#advanced').evaluate(el => { el.open = true; });
  assert.equal(await page.locator('#grid-refinement').inputValue(), 'neighbors');
  assert.equal(await page.locator('#core-rule').inputValue(), 'coverage');
  assert.equal(await page.locator('#activity-confidence').isChecked(), true);
  assert.equal(await page.locator('#merge-identical-documents').isChecked(), true);
  await page.screenshot({path: path.join(out, '01_settings.png')});
  await page.locator('#build').click();
  await page.locator('#group-controls').waitFor({state: 'visible', timeout: 360000});
  const cells = await page.locator('#map-stage [data-grid-row]').count();
  assert.equal(cells, 889);
  const fills = () => page.locator('#map-stage').evaluate(el =>
    [...el.querySelectorAll('[data-grid-row] > rect, circle[data-grid-row]')].map(n => n.getAttribute('fill')));
  const activity = await fills();
  await page.locator('#color-mode').selectOption('support');
  const support = await fills();
  const faded = support.filter((fill, i) => fill !== activity[i]).length;
  assert.ok(faded > 100 && faded < 889, `uncertain cells recolored: ${faded}`);
  assert.match(await page.locator('#legend').innerText(), /uncertain/i);
  await page.screenshot({path: path.join(out, '02_class_support.png')});
  await page.locator('#color-mode').selectOption('reliability');
  const reliability = new Set(await fills());
  assert.ok(reliability.size >= 4);
  assert.match(await page.locator('#legend').innerText(), /representation neighbors/i);
  await page.screenshot({path: path.join(out, '03_neighbor_reliability.png')});
  await page.locator('#color-mode').selectOption('activity');
  await page.locator('#molecule-id').fill('CHEMBL20');
  await page.locator('#find-molecule button').click();
  await page.locator('#selected-id').filter({hasText: 'CHEMBL20'}).waitFor();
  const evidence = await page.locator('#evidence-summary').innerText();
  assert.match(evidence, /Class support\s+61% active · uncertain/);
  assert.match(evidence, /Evidence units\s+27 \(27 documents\)/);
  assert.match(evidence, /Neighbor recall/);
  await page.screenshot({path: path.join(out, '04_acetazolamide_inspector.png')});
  await page.locator('#group-select').selectOption('G002');
  await page.locator('#compare-group').click();
  await page.locator('#comparison-dialog[open]').waitFor();
  const coverage = await page.locator('#core-coverage').innerText();
  assert.match(coverage, /matches 27\/28 members/);
  assert.match(coverage, /at least 80% of members \(target met\)/);
  await page.screenshot({path: path.join(out, '05_G002_core_rule.png')});
  await page.locator('#comparison-dialog [data-close-dialog]').click();
  const state = await page.evaluate(async () => (await fetch('/api/job', {
    headers: {'X-ChemGridMap-Token': sessionStorage.getItem('chemgridmap-token')}})).json());
  const layers = ['class_support', 'neighbor_recall', 'structure_group'].map(mode => `map_${mode}.svg`);
  const mainMaps = ['svg', 'png', 'pdf'].map(extension => `map.${extension}`);
  assert.deepEqual(state.files.filter(name => mapFigure.test(name)).sort(), [...mainMaps, ...layers].sort());
  assert.ok(state.files.includes('map_document_offsets.csv'));
  await page.locator('#open-export').click();
  assert.equal(await page.locator('#figure-downloads button').count(), 6);
  for (const extension of ['SVG', 'PNG', 'PDF']) {
    assert.equal(await page.locator('#figure-downloads button').filter({hasText: `Observed activity · ${extension}`}).count(), 1);
  }
  for (const title of ['Activity class support', 'Neighbor reliability', 'Structure group']) {
    assert.equal(await page.locator('#figure-downloads button').filter({hasText: `${title} · SVG`}).count(), 1);
  }
  await page.screenshot({path: path.join(out, '06_exports.png')});
  const groupedDownload = page.waitForEvent('download');
  await page.locator('#download-all').click();
  const groupedZip = path.join(out, 'grouped_results.zip');
  await (await groupedDownload).saveAs(groupedZip);
  const groupedNames = JSON.parse(execFileSync(python, ['-c',
    'import json,sys,zipfile; print(json.dumps(zipfile.ZipFile(sys.argv[1]).namelist()))', groupedZip], {encoding: 'utf8'}));
  assert.deepEqual(groupedNames.filter(name => mapFigure.test(name)).sort(), [...mainMaps, ...layers].sort());
  await page.locator('#export-dialog [data-close-dialog]').click();
  report.checks.caii = {cells, uncertain_cells_recolored: faded, reliability_colors: [...reliability],
    acetazolamide_evidence: evidence, g002_core_text: coverage, exported_files: state.files};

  // A small source-record fixture gives a known offset and also checks ungrouped maps.
  const lines = ['Smiles,pChEMBL Value,Molecule ChEMBL ID,Document ChEMBL ID'];
  ['CCO', 'CCN', 'CCC'].forEach((smiles, i) => {
    lines.push(`${smiles},${5 + i * .2},CHEMBL${i + 1},DOC_A`);
    lines.push(`${smiles},${8.5 + i * .2},CHEMBL${i + 1},DOC_B`);
  });
  lines.push('CCO,5,CHEMBL1,DOC_COPY');
  const fixture = path.join(out, 'document_offset_fixture.csv');
  fs.writeFileSync(fixture, lines.join('\n') + '\n');
  await page.locator('#adjust-settings').click();
  await page.locator('#replace-input').click();
  await page.locator('#csv').setInputFiles(fixture);
  await page.locator('#page-configure').waitFor({state: 'visible'});
  await page.locator('#input-format').selectOption('chembl');
  await page.locator('#projection').selectOption('pca');
  await page.locator('#grouping').selectOption('none');
  assert.equal(await page.locator('[data-core-setting]').isVisible(), false);
  await page.locator('#merge-identical-documents').uncheck();
  await page.locator('#build').click();
  await page.locator('#group-controls').waitFor({state: 'visible', timeout: 90000});
  assert.equal(await page.locator('#group-select').isVisible(), false);
  assert.equal(await page.locator('#compare-group').isVisible(), false);
  assert.equal(await page.locator('#color-mode').isEnabled(), true);
  await page.locator('#color-mode').selectOption('reliability');
  assert.match(await page.locator('#legend').innerText(), /representation neighbors/i);
  await page.locator('#molecule-id').fill('CHEMBL1');
  await page.locator('#find-molecule button').click();
  await page.locator('[data-document-offset-warning]').waitFor();
  const warning = await page.locator('[data-document-offset-warning]').innerText();
  assert.match(warning, /DOC_A vs DOC_B \(median offset 3.5 log units, 3 shared entries\)/);
  assert.match(await page.locator('#evidence-summary').innerText(), /Evidence units\s+3 \(3 documents\)/);
  await page.screenshot({path: path.join(out, '07_offset_warning_no_grouping.png')});
  await page.locator('#open-export').click();
  const downloaded = page.waitForEvent('download');
  await page.locator('#download-all').click();
  const zipPath = path.join(out, 'ungrouped_results.zip');
  await (await downloaded).saveAs(zipPath);
  const contents = JSON.parse(execFileSync(python, ['-c',
    'import json,sys,zipfile; z=zipfile.ZipFile(sys.argv[1]); print(json.dumps({"names":z.namelist(),"offsets":z.read("map_document_offsets.csv").decode()}))', zipPath], {encoding: 'utf8'}));
  assert.deepEqual(contents.names.filter(name => mapFigure.test(name)).sort(),
    [...mainMaps, 'map_class_support.svg', 'map_neighbor_recall.svg'].sort());
  assert.match(contents.offsets, /DOC_A,DOC_B,3,3.5/);
  assert.deepEqual(report.page_errors, []);
  report.checks.ungrouped = {warning, export_contents: contents};
  report.status = 'passed';
})().catch(error => {
  const text = String(error.stack || error);
  report.status = /Operation not permitted|PermissionError|EPERM|sandbox.*denied/i.test(text) ? 'sandbox-blocked' : 'failed';
  report.error = text.replace(/http:\/\/127\.0\.0\.1:\d+\/#[^\s]+/g, '[local session URL]');
  process.exitCode = report.status === 'sandbox-blocked' ? 2 : 1;
}).finally(async () => {
  if (browser) await browser.close();
  if (server) server.kill();
  fs.writeFileSync(path.join(out, 'ui_check.json'), JSON.stringify(report, null, 2));
  manifest();
  console.log(JSON.stringify({status: report.status, report: path.join(out, 'ui_check.json')}));
});
