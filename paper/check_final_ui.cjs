// Check the installed application's native path; never modify page content for screenshots.
const {chromium} = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');

(async () => {
  const [url, output, input] = process.argv.slice(2);
  const out = path.resolve(output);
  fs.mkdirSync(out, {recursive: true});
  const browser = await chromium.launch({channel: 'chrome', headless: true});
  const page = await browser.newPage({viewport: {width: 1440, height: 1050}, deviceScaleFactor: 2});
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(url);
  await page.locator('#version').filter({hasText:'0.4.0'}).waitFor();
  await page.locator('#step-import').click();
  if (input === 'demo') await page.locator('#demo').click();
  else await page.locator('#csv').setInputFiles(path.resolve(input));
  await page.locator('#page-configure').waitFor({state: 'visible'});
  await page.locator('#projection').selectOption('umap');
  await page.locator('#build').click();
  await page.locator('#group-controls').waitFor({state:'visible', timeout:180000});
  const positions = () => page.locator('#map-stage').evaluate(el =>
    [...el.querySelectorAll('[data-grid-row]')].map(n => [n.dataset.gridRow,n.dataset.gridCol,n.getAttribute('transform')]));
  const before = await positions();
  assert.equal(before.length,889);
  await page.locator('#color-mode').selectOption('structure');
  await page.locator('#group-select').selectOption('G002');
  const comparisonResponse = page.waitForResponse(response => new URL(response.url()).pathname === '/api/compare');
  await page.locator('#compare-group').click();
  const comparison = await (await comparisonResponse).json();
  await page.locator('#comparison-dialog').waitFor({state:'visible'});
  await page.locator('#comparison-svg svg').waitFor();
  const coverage = await page.locator('#core-coverage').innerText();
  const preview = await page.locator('#comparison-rule').innerText();
  assert.match(coverage,/27\/28/);
  assert.match(preview,/12\/28 members shown/);
  const previewIds = await page.locator('#comparison-svg [data-comparison-id]').evaluateAll(nodes => nodes.map(n => n.dataset.comparisonId));
  assert.equal(previewIds.length,12);
  await page.locator('#comparison-dialog').screenshot({path:path.join(out,'native_group_preview.png')});
  const downloadPromise = page.waitForEvent('download');
  await page.locator('#download-comparison').click();
  await (await downloadPromise).saveAs(path.join(out,'native_comparison.zip'));
  await page.getByText('All group members and scaffold coverage',{exact:true}).click();
  const rows = page.locator('#group-members tbody tr');
  assert.equal(await rows.count(),28);
  const outside = comparison.members.find(member => !previewIds.includes(member.comparison_id));
  assert.ok(outside);
  await rows.filter({hasText:outside.comparison_id}).click();
  await page.locator('#selected-id').filter({hasText:outside.molecule_id}).waitFor();
  await page.locator('#selected-group-button').click();
  await page.locator('#comparison-dialog').waitFor({state:'visible'});
  const memberDetails = page.locator('#group-members').locator('..');
  if (!(await memberDetails.evaluate(el => el.open)))
    await page.getByText('All group members and scaffold coverage',{exact:true}).click();
  const unmatched = comparison.members.find(member => String(member.molecule_id).split(';').includes('CHEMBL85636'));
  assert.ok(unmatched);
  const exception = rows.filter({hasText:unmatched.comparison_id});
  await exception.scrollIntoViewIfNeeded();
  assert.match(await exception.innerText(),/false/);
  await page.locator('#comparison-dialog').screenshot({path:path.join(out,'native_full_members.png')});
  await exception.click();
  await page.locator('#selected-id').filter({hasText:'CHEMBL85636'}).waitFor();
  await page.locator('#molecule-svg svg').waitFor();
  assert.equal(await page.locator('#comparison-dialog').isVisible(),false);
  await page.screenshot({path:path.join(out,'native_unmatched_inspector.png')});
  const unmatchedNeighbors = await page.locator('#neighbor-summary').textContent();
  assert.match(unmatchedNeighbors,/Recall/);
  await page.locator('#download-audit').scrollIntoViewIfNeeded();
  const evidencePromise = page.waitForEvent('download');
  await page.locator('#download-audit').click();
  await (await evidencePromise).saveAs(path.join(out,'native_unmatched_evidence.zip'));
  await page.locator('#molecule-id').fill('CHEMBL20');
  await page.getByRole('button',{name:'Find molecule',exact:true}).click();
  await page.locator('#selected-id').filter({hasText:'CHEMBL20'}).waitFor();
  await page.locator('#view-records').click();
  await page.locator('#records-dialog').waitFor({state:'visible'});
  assert.equal(await page.locator('#record-table tbody tr').count(),54);
  await page.locator('#records-dialog').screenshot({path:path.join(out,'native_source_records.png')});
  await page.getByRole('button',{name:'Close source records',exact:true}).click();
  const recordsPromise = page.waitForEvent('download');
  await page.locator('#download-evidence').click();
  await (await recordsPromise).saveAs(path.join(out,'native_acetazolamide_records.csv'));
  assert.deepEqual(await positions(),before);
  assert.deepEqual(errors,[]);
  fs.writeFileSync(path.join(out,'native_ui_checks.json'),JSON.stringify({
    passed:true, software:'0.4.0', input_csv:input === 'demo' ? 'bundled chembl205_ic50.csv' : path.basename(input), molecular_entries:889,
    preview_count:12, total_members:28, core_coverage:coverage, preview_rule:preview,
    preview_ids:previewIds, unmatched_in_default_preview:previewIds.includes(unmatched.comparison_id),
    unmatched_comparison_key:unmatched.comparison_id,
    outside_preview_member:outside.molecule_id, outside_preview_access_verified:true,
    unmatched_id:'CHEMBL85636', full_members_table_count:28, unmatched_access:'All group members table -> select row -> molecule inspector',
    unmatched_neighbors:unmatchedNeighbors, source_records_displayed:54,
    downloads:['native_comparison.zip','native_unmatched_evidence.zip','native_acetazolamide_records.csv'],
    unchanged_grid_positions:true, screenshots:'Unmodified native UI; element screenshots and viewport captures only',
    browser:await browser.version(),javascript_errors:errors
  },null,2));
  await browser.close();
})().catch(error => {console.error(error);process.exit(1);});
