// Exercise the real browser workflow against a running local test server.
const {chromium} = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');

(async () => {
  const out = path.resolve(process.argv[3]);
  fs.mkdirSync(out, {recursive: true});
  const browser = await chromium.launch({channel: 'chrome', headless: true});
  const page = await browser.newPage({viewport: {width: 1440, height: 1050}});
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(process.argv[2]);
  await page.locator('#version').filter({hasText:'0.4.0'}).waitFor();
  await page.locator('#step-import').click();
  await page.locator('#demo').click();
  await page.locator('#page-configure').waitFor({state: 'visible'});
  await page.locator('#projection').selectOption('umap');
  assert.equal(await page.locator('#options').evaluate(el => el.checkValidity()), true);
  await page.locator('#build').click();
  await page.locator('#group-controls').waitFor({state: 'visible', timeout: 180000});
  const positions = () => page.locator('#map-stage').evaluate(el =>
    [...el.querySelectorAll('[data-grid-row]')].map(n => [n.getAttribute('data-grid-row'), n.getAttribute('data-grid-col'), n.getAttribute('transform')]));
  const before = await positions();
  assert.equal(before.length, 889);
  await page.locator('#color-mode').selectOption('structure');
  assert.deepEqual(await positions(), before);
  assert.equal(await page.locator('#group-select option').count(), 314);
  await page.screenshot({path: path.join(out, 'desktop_groups.png')});
  await page.locator('#group-select').selectOption('G002');
  await page.locator('#compare-group').click();
  await page.locator('#comparison-dialog').waitFor({state: 'visible'});
  await page.locator('#comparison-svg svg').waitFor();
  assert.match(await page.locator('#core-coverage').innerText(), /27\s*\/\s*28/);
  await page.screenshot({path: path.join(out, 'desktop_comparison.png')});
  const downloadPromise = page.waitForEvent('download');
  await page.locator('#download-comparison').click();
  const download = await downloadPromise;
  await download.saveAs(path.join(out, 'comparison.zip'));
  await page.getByRole('button', {name: 'Close structural comparison', exact: true}).click();
  await page.locator('#molecule-id').fill('CHEMBL20');
  await page.getByRole('button', {name: 'Find molecule', exact: true}).click();
  await page.locator('#selected-id').filter({hasText:'CHEMBL20'}).waitFor();
  await page.locator('#neighbor-key').waitFor({state:'visible'});
  await page.screenshot({path:path.join(out,'desktop_neighbors.png')});
  await page.setViewportSize({width: 390, height: 844});
  await page.screenshot({path: path.join(out, 'mobile_map.png'), fullPage:true});
  await page.locator('#close-inspector').click();
  await page.locator('#group-select').selectOption('G002');
  await page.locator('#compare-group').click();
  await page.locator('#comparison-dialog').waitFor({state:'visible'});
  await page.screenshot({path: path.join(out, 'mobile_comparison.png')});
  const overflow = await page.evaluate(() => document.documentElement.scrollWidth > innerWidth);
  assert.equal(overflow, false);
  assert.deepEqual(errors, []);
  fs.writeFileSync(path.join(out,'browser_checks.json'), JSON.stringify({
    passed:true, groups:314, core_coverage:'27/28', coordinate_transforms_unchanged: true,
    tested_molecules:before.length, comparison_download:true, mobile_overflow:false,
    javascript_errors:errors, browser:await browser.version()
  },null,2));
  await browser.close();
})().catch(error => { console.error(error); process.exit(1); });
