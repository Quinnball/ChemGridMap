#!/usr/bin/env node
// T1 UI capture. This script starts its own app and is intended for an
// unrestricted desktop run by Opus; it is not executed in the timing sandbox.
const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {spawn} = require('node:child_process');

const PAPER = path.resolve(__dirname, '../..');
const CGM = path.resolve(PAPER, '..');
const python = path.join(CGM, '.venv/bin/python');
const defaultOutput = path.join(PAPER, 'strengthening/output/T1/ui');
const appOutput = '/tmp/chemgridmap-demo';

function usage() {
  console.log([
    'Usage: node t1_ui_capture.cjs [--output-dir DIR] [--url URL]',
    '',
    'Captures a 2x CA II class-support map and CHEMBL18 inspector evidence.',
  ].join('\n'));
}

const args = process.argv.slice(2);
if (args.includes('--help') || args.includes('-h')) {
  usage();
  process.exit(0);
}
let output = defaultOutput;
let suppliedURL = null;
for (let i = 0; i < args.length; i += 1) {
  if (args[i] === '--output-dir' && args[i + 1]) output = path.resolve(args[++i]);
  else if (args[i] === '--url' && args[i + 1]) suppliedURL = args[++i];
  else throw new Error(`Unknown argument: ${args[i]}`);
}
fs.mkdirSync(output, {recursive: true});

let server;
let serverLog = '';
const report = {
  status: 'running',
  device_scale_factor: 2,
  selected_entry: 'CHEMBL18',
  selected_group: 'G002',
  screenshots: [],
  checks: {},
  page_errors: [],
};

function digest(value) {
  return crypto.createHash('sha256').update(value).digest('hex');
}

async function startApp() {
  server = spawn(
    python,
    ['-m', 'chemgridmap.app', '--no-browser', '--output-dir', appOutput],
    {cwd: CGM, env: {...process.env, PYTHONPATH: path.join(CGM, 'src')}},
  );
  return new Promise((resolve, reject) => {
    server.stdout.on('data', chunk => {
      serverLog += chunk.toString();
      const match = serverLog.match(/http:\/\/127\.0\.0\.1:\d+\/#[^\s]+/);
      if (match) resolve(match[0]);
    });
    server.stderr.on('data', chunk => { serverLog += chunk.toString(); });
    server.on('error', reject);
    server.on('exit', code => {
      if (!serverLog.includes('/#')) reject(new Error(`Local app exited (${code})`));
    });
  });
}

async function capture(page, name) {
  const file = path.join(output, name);
  await page.screenshot({path: file, type: 'png'});
  report.screenshots.push({file: name, sha256: digest(fs.readFileSync(file))});
}

async function run() {
  const url = suppliedURL || await startApp();
  const browser = await chromium.launch({channel: 'chrome', headless: true});
  try {
    const page = await browser.newPage({
      viewport: {width: 1440, height: 1200},
      deviceScaleFactor: 2,
    });
    page.on('pageerror', error => report.page_errors.push(error.message));
    await page.goto(url);
    await page.locator('#version').filter({hasText: /v0\.5\.0/}).waitFor();
    await page.locator('#demo').click();
    await page.locator('#page-configure').waitFor({state: 'visible'});
    await page.locator('#projection').selectOption('umap');
    await page.locator('#advanced').evaluate(element => { element.open = true; });
    assert.equal(await page.locator('#grid-refinement').inputValue(), 'neighbors');
    assert.equal(await page.locator('#core-rule').inputValue(), 'coverage');
    assert.equal(await page.locator('#activity-confidence').isChecked(), true);
    assert.equal(await page.locator('#merge-identical-documents').isChecked(), true);

    await page.locator('#build').click();
    await page.locator('#group-controls').waitFor({state: 'visible', timeout: 600000});
    const fills = () => page.locator('#map-stage').evaluate(element =>
      [...element.querySelectorAll('[data-grid-row] > rect, circle[data-grid-row]')]
        .map(node => node.getAttribute('fill')));
    const activity = await fills();
    await page.locator('#color-mode').selectOption('support');
    const support = await fills();
    const recolored = support.filter((fill, index) => fill !== activity[index]).length;
    assert.ok(recolored > 100 && recolored < 889, `support recolored ${recolored} cells`);
    assert.match(await page.locator('#legend').innerText(), /uncertain/i);
    await capture(page, '01_caii_class_support.png');

    // CHEMBL18 is a non-CHEMBL20 multi-unit G002 member with four documents.
    await page.locator('#molecule-id').fill('CHEMBL18');
    await page.locator('#find-molecule button').click();
    await page.locator('#selected-id').filter({hasText: 'CHEMBL18'}).waitFor();
    const evidence = await page.locator('#evidence-summary').innerText();
    assert.match(evidence, /Evidence units\s+4 \(4 documents\)/);
    assert.match(evidence, /Source records\s+\d+/);
    assert.match(evidence, /Class support\s+\d+% active · supported/);
    assert.match(evidence, /Neighbor recall\s+\d+\/\d+/);
    const groupSummary = await page.locator('#selected-group-summary').innerText();
    assert.match(groupSummary, /G002 \| 28 members \| core covers 27\/28/);
    await page.locator('#inspector details').evaluateAll(elements => {
      elements.forEach(element => { element.open = true; });
    });
    const neighborSummary = await page.locator('#neighbor-summary').innerText();
    assert.match(neighborSummary, /shared neighbors out of/);
    const bodyText = await page.locator('body').innerText();
    assert.doesNotMatch(bodyText, new RegExp('\\x2fUsers\\x2f|chemical map\\x2fchemgridmap', 'i'));
    await capture(page, '02_CHEMBL18_inspector.png');

    await page.locator('#view-records').click();
    await page.locator('#records-dialog[open]').waitFor();
    const recordText = await page.locator('#record-table').innerText();
    assert.match(recordText, /document_chembl_id/);
    assert.match(recordText, /CHEMBL\d+/);
    await capture(page, '03_CHEMBL18_source_records.png');
    await page.locator('#records-dialog [data-close-dialog]').click();

    await page.locator('#group-select').selectOption('G002');
    await page.locator('#compare-group').click();
    await page.locator('#comparison-dialog[open]').waitFor();
    const core = await page.locator('#core-coverage').innerText();
    assert.match(core, /matches 27\/28 members/);
    assert.match(core, /at least 80% of members \(target met\)/);
    const dialogText = await page.locator('#comparison-dialog').innerText();
    assert.doesNotMatch(dialogText, new RegExp('\\x2fUsers\\x2f|chemical map\\x2fchemgridmap', 'i'));
    await capture(page, '04_G002_core_coverage.png');
    report.checks.caii = {
      cells: await page.locator('#map-stage [data-grid-row]').count(),
      support_cells_recolored: recolored,
      evidence: evidence.replace(/\s+/g, ' ').trim(),
      group_summary: groupSummary,
      neighbor_summary: neighborSummary.replace(/\s+/g, ' ').trim(),
      core_coverage: core.replace(/\s+/g, ' ').trim(),
    };
    assert.deepEqual(report.page_errors, []);
    report.status = 'passed';
  } finally {
    await browser.close();
  }
}

(async () => {
  try {
    await run();
  } catch (error) {
    const text = String(error.stack || error);
    report.status = /Operation not permitted|PermissionError|EPERM|sandbox.*denied/i.test(text)
      ? 'sandbox-blocked' : 'failed';
    report.error = text.replace(/http:\/\/127\.0\.0\.1:\d+\/#[^\s]+/g, '[local session URL]');
    process.exitCode = report.status === 'sandbox-blocked' ? 2 : 1;
  } finally {
    if (server) server.kill();
    fs.writeFileSync(
      path.join(output, 'capture_report.json'),
      JSON.stringify(report, null, 2) + '\n',
    );
    console.log(JSON.stringify({status: report.status, report: path.join(output, 'capture_report.json')}));
  }
})();
