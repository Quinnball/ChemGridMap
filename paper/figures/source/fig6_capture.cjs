#!/usr/bin/env node
// F2 图 6 截图。仅由 Opus 在沙箱外运行；本任务只检查语法。
const {chromium} = require('playwright');
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const {spawn} = require('node:child_process');

const CGM = path.resolve(__dirname, '../../..');
const python = path.join(CGM, '.venv/bin/python');
const defaultOutput = path.resolve(__dirname, '../data_fig6');
const appOutput = '/tmp/chemgridmap-f2-demo';

function usage() {
  console.log([
    'Usage: node fig6_capture.cjs [--output-dir DIR] [--url URL]',
    '',
    'Captures the CA II class-support map with its legend, the CHEMBL18 inspector, its source records and a neighbor view.',
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
  device_scale_factor: 3,
  viewport: {width: 1440, height: 1600},
  selected_entry: 'CHEMBL18',
  selected_group: 'G002',
  screenshots: [],
  checks: {},
  page_errors: [],
};

function digest(value) {
  return crypto.createHash('sha256').update(value).digest('hex');
}

function screenshotRecord(file) {
  return {file: path.basename(file), sha256: digest(fs.readFileSync(file))};
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

async function capture(page, file, options) {
  const {locator, ...screenshotOptions} = options || {};
  if (locator) await locator.screenshot({path: file, type: 'png', ...screenshotOptions});
  else await page.screenshot({path: file, type: 'png', ...screenshotOptions});
  report.screenshots.push(screenshotRecord(file));
}

async function tileBox(tile) {
  return tile.evaluate(element => {
    const rect = element.getBoundingClientRect();
    return {x: rect.x, y: rect.y, width: rect.width, height: rect.height};
  });
}

async function zoomTileTo(page, tile, targetCSSPx) {
  let box = await tileBox(tile);
  const center = {x: box.x + box.width / 2, y: box.y + box.height / 2};
  await page.mouse.move(center.x, center.y);
  // Chrome may deliver a large wheel delta only in part, so step until the tile reaches the target.
  for (let i = 0; i < 15; i += 1) {
    box = await tileBox(tile);
    if (Math.abs(box.width - targetCSSPx) <= 0.2) return box;
    const factor = targetCSSPx / box.width;
    const deltaY = -Math.log(factor) / 0.002;
    await page.mouse.wheel(0, deltaY);
    await page.waitForTimeout(150);
  }
  return tileBox(tile);
}

async function run() {
  const url = suppliedURL || await startApp();
  const browser = await chromium.launch({channel: 'chrome', headless: true});
  try {
    const page = await browser.newPage({
      viewport: report.viewport,
      deviceScaleFactor: report.device_scale_factor,
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
    await page.locator('#color-mode').selectOption('support');
    // The class-support legend lists the uncertain shades; it has no "Activity class support" title.
    const legendText = await page.locator('#legend').innerText();
    assert.match(legendText, /uncertain/i);
    assert.match(legendText, /class support below 0\.8/);
    // CA II has no flagged document pair, so the offset marker item must not be listed
    // (the explanatory sentence still mentions offset pairs).
    assert.doesNotMatch(legendText, /◣/);
    report.checks.legend = legendText.replace(/\s+/g, ' ').trim();

    const controlsBox = await page.locator('#group-controls').boundingBox();
    const cardBox = await page.locator('.map-card').boundingBox();
    const mapBottomBox = await page.locator('.map-card .map-bottom').boundingBox();
    assert.ok(controlsBox && cardBox && mapBottomBox);
    // Include the legend row below the viewport.
    const mapClip = {
      x: cardBox.x,
      y: controlsBox.y,
      width: cardBox.width,
      height: mapBottomBox.y + mapBottomBox.height - controlsBox.y,
    };
    await capture(page, path.join(output, '6a_map.png'), {clip: mapClip});
    report.map_clip_css = mapClip;
    // Where the map viewport and the legend row lie inside 6a_map.png, in CSS px from its top left.
    const mapViewportBox = await page.locator('#map-viewport').boundingBox();
    assert.ok(mapViewportBox);
    report.map_viewport_in_clip_css = {
      x: mapViewportBox.x - mapClip.x,
      y: mapViewportBox.y - mapClip.y,
      width: mapViewportBox.width,
      height: mapViewportBox.height,
    };
    report.map_legend_in_clip_css = {
      x: mapBottomBox.x - mapClip.x,
      y: mapBottomBox.y - mapClip.y,
      width: mapBottomBox.width,
      height: mapBottomBox.height,
    };

    await page.locator('#molecule-id').fill('CHEMBL18');
    await page.locator('#find-molecule button').click();
    await page.locator('#selected-id').filter({hasText: 'CHEMBL18'}).waitFor();
    const evidence = await page.locator('#evidence-summary').innerText();
    assert.match(evidence, /Evidence units\s+4 \(4 documents\)/);
    assert.match(evidence, /Source records\s+\d+/);
    assert.match(evidence, /Class support\s+(> 99%|\d+%) active · supported/);
    assert.match(evidence, /Support status\s+supported/);
    assert.match(evidence, /σ_doc\s+0\.796/);
    assert.match(evidence, /Neighbor recall \(tie-inclusive\)\s+\d+\/\d+/);
    assert.match(evidence, /Neighbor recall \(fixed k\)\s+\d+%/);
    assert.match(evidence, /10th fingerprint neighbor\s+Tanimoto 0\.\d\d/);
    // CHEMBL18 has a wide record range; the inspector explains it neutrally.
    assert.match(evidence, /Records span more than one log unit/);
    assert.doesNotMatch(evidence, /Conflict flagged/);
    const groupSummary = await page.locator('#selected-group-summary').innerText();
    assert.match(groupSummary, /G002 \| 28 members \| core covers 27\/28/);
    // The neighbor summary sits in a collapsed <details>; innerText is empty until it is opened.
    await page.locator('#inspector details').evaluateAll(elements => {
      elements.forEach(element => { element.open = true; });
    });
    const neighborSummary = await page.locator('#neighbor-summary').innerText();
    assert.match(neighborSummary, /shared neighbors out of \d+ fingerprint neighbors/);
    assert.equal((await page.locator('#neighbor-title').innerText()).trim(), 'Fingerprint vs grid neighbors');
    assert.equal(await page.locator('#show-neighbor-links').isChecked(), true);
    const selectedViewportBox = await page.locator('#map-viewport').boundingBox();
    assert.ok(selectedViewportBox);

    const inspector = page.locator('#inspector');
    await inspector.evaluate(element => {
      element.style.maxHeight = 'none';
      element.style.overflow = 'visible';
    });
    const crop = await page.evaluate(() => {
      const inspectorElement = document.querySelector('#inspector');
      const summary = document.querySelector('#evidence-summary');
      const label = [...document.querySelectorAll('#evidence-summary dt')]
        .find(element => element.textContent.trim() === 'Neighbor recall (fixed k)');
      if (!inspectorElement || !summary || !label) throw new Error('Inspector rows not found');
      const top = inspectorElement.getBoundingClientRect().top;
      const labelRect = label.getBoundingClientRect();
      const value = label.nextElementSibling?.getBoundingClientRect();
      const recallBottom = Math.max(labelRect.bottom, value ? value.bottom : labelRect.bottom);
      const heading = document.querySelector('#selected-id');
      return {
        // Top of the entry heading (CHEMBL18), so that the group header above it can be cropped away.
        heading: heading ? heading.getBoundingClientRect().top - top : 0,
        recall: recallBottom - top,
        // The bottom of the evidence summary includes the wide-record-range note.
        evidence: summary.getBoundingClientRect().bottom - top,
      };
    });
    report.neighbor_recall_bottom_css_px = crop.recall;
    report.evidence_bottom_css_px = crop.evidence;
    report.heading_top_css_px = crop.heading;
    await capture(page, path.join(output, '6b_inspector.png'), {locator: inspector});
    // The element screenshot scrolls the page; restore the layout before measuring the map again.
    await inspector.evaluate(element => {
      element.style.maxHeight = '';
      element.style.overflow = '';
    });
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.waitForTimeout(150);

    // The source records behind CHEMBL18, as the "View source records" dialog shows them.
    await page.locator('#view-records').click();
    await page.locator('#record-section').waitFor({state: 'visible'});
    const recordsDialog = page.locator('#records-dialog');
    const recordText = await recordsDialog.innerText();
    assert.match(recordText, /Source activity records/);
    assert.match(recordText, /CHEMBL\d+/);
    assert.doesNotMatch(recordText, new RegExp('\\x2fUsers\\x2f|\\x2ftmp\\x2f|\\x2fprivate\\x2f', 'i'));
    report.checks.records = {
      count_line: (await page.locator('#record-count').innerText()).replace(/\s+/g, ' ').trim(),
      rows: await page.locator('#record-table tbody tr').count(),
    };
    await capture(page, path.join(output, '6d_records.png'), {locator: recordsDialog});
    // A crop of the first five columns, legible at print size; the assay description is omitted.
    const tableClip = await page.evaluate(() => {
      const table = document.querySelector('#record-table');
      const headers = [...document.querySelectorAll('#record-table thead th')];
      if (!table || headers.length < 5) throw new Error('Record table columns not found');
      const tableRect = table.getBoundingClientRect();
      const fifth = headers[4].getBoundingClientRect();
      return {
        x: tableRect.x,
        y: tableRect.y,
        width: fifth.right - tableRect.x,
        height: tableRect.height,
        columns: headers.slice(0, 5).map(element => element.textContent.trim()),
      };
    });
    assert.deepEqual(tableClip.columns, [
      'source_row', 'molecule_identity_key', 'activity_pchembl_raw', 'assay_chembl_id', 'document_chembl_id',
    ]);
    report.records_table_clip_css = tableClip;
    await capture(page, path.join(output, '6d_records_table.png'), {
      clip: {x: tableClip.x, y: tableClip.y, width: tableClip.width, height: tableClip.height},
    });
    await page.locator('#records-dialog [data-close-dialog]').click();
    await recordsDialog.waitFor({state: 'hidden'});
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.waitForTimeout(150);

    const savedLocation = page.locator('#saved-location');
    if (await savedLocation.count()) await savedLocation.evaluate(element => { element.textContent = ''; });
    const bodyText = await page.locator('body').innerText();
    assert.doesNotMatch(bodyText, new RegExp('\\x2fUsers\\x2f|\\x2ftmp\\x2f|\\x2fprivate\\x2f|chemical map\\x2fchemgridmap', 'i'));

    const tile = page.locator('#map-stage .selected-cell').first();
    await tile.waitFor({state: 'visible'});
    let box = await zoomTileTo(page, tile, 23.0);
    assert.ok(box.width >= 22 && box.width <= 23.9, `zoomed tile width ${box.width}`);
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    box = await tileBox(tile);
    // Size the square window to the farthest outlined neighbor (CHEMBL18 has one
    // fingerprint-only neighbor more than three cells away), with one cell of margin.
    const reach = await page.locator('#map-stage [data-grid-row]').evaluateAll(elements => {
      const selectedCell = elements.find(element => element.classList.contains('selected-cell'));
      const row = Number(selectedCell.dataset.gridRow), col = Number(selectedCell.dataset.gridCol);
      let farthest = 0;
      for (const element of elements) {
        if (!['neighbor-shared', 'neighbor-grid_only', 'neighbor-representation_only']
          .some(name => element.classList.contains(name))) continue;
        farthest = Math.max(farthest, Math.abs(Number(element.dataset.gridRow) - row),
          Math.abs(Number(element.dataset.gridCol) - col));
      }
      return farthest;
    });
    const gridCells = 2 * (reach + 1) + 1;
    const side = gridCells * box.width;
    const center = {x: box.x + box.width / 2, y: box.y + box.height / 2};
    const neighborClip = {
      x: center.x - side / 2,
      y: center.y - side / 2,
      width: side,
      height: side,
    };
    const zoomedViewportBox = await page.locator('#map-viewport').boundingBox();
    assert.ok(zoomedViewportBox);
    assert.ok(
      neighborClip.x >= zoomedViewportBox.x
        && neighborClip.y >= zoomedViewportBox.y
        && neighborClip.x + side <= zoomedViewportBox.x + zoomedViewportBox.width
        && neighborClip.y + side <= zoomedViewportBox.y + zoomedViewportBox.height,
      `neighbor clip ${JSON.stringify(neighborClip)} extends outside map viewport ${JSON.stringify(zoomedViewportBox)}`,
    );
    const selected = await page.locator('#map-stage [data-grid-row].selected-cell').evaluateAll(
      elements => elements.map(element => {
        const rect = element.getBoundingClientRect();
        return {x: rect.x, y: rect.y, width: rect.width, height: rect.height};
      }),
    );
    const neighborCounts = {
      selected_cell: selected.length,
      'neighbor-shared': 0,
      'neighbor-grid_only': 0,
      'neighbor-representation_only': 0,
    };
    const countMatch = neighborSummary.replace(/\s+/g, ' ').match(
      /(\d+) shared neighbors out of (\d+) fingerprint neighbors and (\d+) grid neighbors/,
    );
    assert.ok(countMatch, `neighbor summary counts not found: ${neighborSummary}`);
    const expectedCounts = {
      selected_cell: 1,
      'neighbor-shared': Number(countMatch[1]),
      'neighbor-grid_only': Number(countMatch[3]) - Number(countMatch[1]),
      'neighbor-representation_only': Number(countMatch[2]) - Number(countMatch[1]),
    };
    assert.equal(selected.length, expectedCounts.selected_cell);
    assert.ok(
      selected.every(rect => (
        rect.x >= neighborClip.x
        && rect.y >= neighborClip.y
        && rect.x + rect.width <= neighborClip.x + neighborClip.width
        && rect.y + rect.height <= neighborClip.y + neighborClip.height
      )),
      'selected cell is not fully contained in the neighbor capture',
    );
    for (const className of [
      'neighbor-shared',
      'neighbor-grid_only',
      'neighbor-representation_only',
    ]) {
      const boxes = await page.locator(`#map-stage [data-grid-row].${className}`).evaluateAll(
        elements => elements.map(element => {
          const rect = element.getBoundingClientRect();
          return {x: rect.x, y: rect.y, width: rect.width, height: rect.height};
        }),
      );
      neighborCounts[className] = boxes.length;
      // Cell positions inside 6c_neighbors.png, in CSS px, for the print labels of Fig. 6c.
      report.neighbor_cells = (report.neighbor_cells || []).concat(boxes.map(rect => ({
        class: className,
        x: rect.x - neighborClip.x,
        y: rect.y - neighborClip.y,
        width: rect.width,
        height: rect.height,
      })));
      assert.equal(
        boxes.length,
        expectedCounts[className],
        `${className} count ${boxes.length} does not match summary ${expectedCounts[className]}`,
      );
      assert.ok(
        boxes.every(rect => (
          rect.x >= neighborClip.x
          && rect.y >= neighborClip.y
          && rect.x + rect.width <= neighborClip.x + neighborClip.width
          && rect.y + rect.height <= neighborClip.y + neighborClip.height
        )),
        `${className} has cells outside the neighbor capture`,
      );
    }
    report.selected_cell_in_clip_css = selected.map(rect => ({
      x: rect.x - neighborClip.x, y: rect.y - neighborClip.y, width: rect.width, height: rect.height,
    }));
    await capture(page, path.join(output, '6c_neighbors.png'), {clip: neighborClip});
    report.neighbor_clip_css = neighborClip;
    report.neighbor_cell_width_css = box.width;
    report.neighbor_clip_grid_cells = gridCells;
    report.neighbor_reach_cells = reach;
    report.neighbor_counts = neighborCounts;
    report.checks.caii = {
      cells: await page.locator('#map-stage [data-grid-row]').count(),
      evidence: evidence.replace(/\s+/g, ' ').trim(),
      group_summary: groupSummary.replace(/\s+/g, ' ').trim(),
      neighbor_summary: neighborSummary.replace(/\s+/g, ' ').trim(),
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
    report.error = text.replace(/http:\/\/127\.0\.0\.1:\d+\/#[^\s#]+/g, '[local session URL]');
    process.exitCode = report.status === 'sandbox-blocked' ? 2 : 1;
  } finally {
    if (server) server.kill();
    fs.writeFileSync(
      path.join(output, 'capture_report.json'),
      `${JSON.stringify(report, null, 2)}\n`,
    );
    console.log(JSON.stringify({status: report.status, report: path.join(output, 'capture_report.json')}));
  }
})();
