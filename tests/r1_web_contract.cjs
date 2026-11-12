'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(process.argv[2], 'utf8');
const elements = new Map(), calls = [];
const dictionary = JSON.parse(fs.readFileSync(require('node:path').join(require('node:path').dirname(process.argv[2]), 'i18n.json'), 'utf8')).en;
const i18n = {getLanguage: () => 'en', t: (key,values={}) => (dictionary[key] || key.slice(2)).replace(/\{(\w+)\}/g,(_,name)=>values[name] ?? `{${name}}`)};
// Only the nodes touched by the tested functions are needed; this is not a DOM emulator.
function node() {
  return {children: [], dataset: {}, style: {setProperty(key, value) { this[key] = value; }},
    classList: {toggle() {}, add() {}, remove() {}}, textContent: '', hidden: false, value: '',
    replaceChildren(...children) { this.children = children; }, append(...children) { this.children.push(...children); },
    setAttribute() {}, get firstElementChild() { return this.children[0]; },
    createTHead() { const head = node(); this.append(head); return head; },
    createTBody() { const body = node(); this.append(body); return body; },
    insertRow() { const row = node(); this.append(row); return row; },
    insertCell() { const cell = node(); this.append(cell); return cell; },
    querySelector() { return null; }, querySelectorAll() { return this.children; }};
}
const $ = id => { if (!elements.has(id)) elements.set(id, node()); return elements.get(id); };
const document = {createElement: node, createTextNode: text => ({textContent: text})};
function realFunction(name) {
  const start = source.search(new RegExp(`^(?:async )?function ${name}\\(`, 'm'));
  assert.ok(start >= 0, `Missing function ${name}`);
  const rest = source.slice(start), next = rest.slice(1).search(/^(?:async )?function \w+\(/m);
  return next < 0 ? rest : rest.slice(0, next + 1);
}
function fragment(functionName, start, end) {
  const body = realFunction(functionName), begin = body.indexOf(start), finish = body.indexOf(end, begin);
  assert.ok(begin >= 0 && finish > begin, `Missing ${functionName} display fragment`);
  return body.slice(begin, finish);
}
const measured = {molecule_id: 'CHEMBL538536', molecule_identity_key: 'CHEMBL1189346', structure_group: 'G161'};
const query = {query_id: 'Q1', molecule_id: 'Q1', structure_group: 'G161', measurement_status: 'not_measured',
  nearest_measured_id: 'CHEMBL1189346', nearest_measured_similarity: .8};
const group = {structure_group: 'G161', size: 1, core_smiles: 'c1ccccc1', core_coverage_count: 1, core_exact_scaffold_count: 1};
const context = vm.createContext({$, document, t:i18n.t, ChemGridMapI18n:i18n, newMolecules:null, points: [measured, query], structureGroups: [group],
  highlightedGroup: '', groupRequestNumber: 0, job: {options: {input_format: 'molecule'}, molecules: 889,
    report: {input_rows: 1100, retained_activity_records: 1100, repeated_conflicted_molecules: 17}},
  transform: {s: 1, x: 0, y: 0}, view: 'grid', tileWidth: 100, hasStructures: false,
  busy: false, requestNumber: 0, requestAnimationFrame() {}, status() {}, fit() {},
  colors: {inactive: '#c44e52', medium: '#e6a23c', active: '#2a9d8f'},
  fadedColors: {inactive: '#ebc0c1', medium: '#f6dfb8', active: '#b8e0db'},
  highlight() {}, updateNavigation() {}, coreRule() { return ''; },
  renderColorLegend() { $('legend').replaceChildren(); }, getJSON: async () => ({core_svg: ''}),
  fail(error) { throw error; }, inspect: query => calls.push(query),
  table(rows, container, limit, columns) { container.columns = columns; container.replaceChildren(...rows.map(data => ({...node(), data}))); },
});
$('color-mode').value = 'activity'; $('measurement-filter').value = 'all';
$('measurement-filter').selectedOptions = [{textContent: 'All molecules'}];
const names = ['renderMeasuredResult', 'resultSource', 'resetMeasurementForm', 'updateMeasurementForm', 'readMeasurement', 'saveMeasuredResult', 'addedFeedbackCounts', 'resultCounts', 'stat', 'countLabel', 'measuredLabel', 'showQuerySelection', 'selectGroup',
  'format', 'percentage', 'legendItem', 'renderLegend', 'renderMismatchOverlay', 'applyTransform', 'showDataAudit',
  'showSourceInput', 'loadInput', 'showInput', 'selectColumns', 'showSources', 'renderNewBatchSources', 'clearSourceHighlight', 'highlightSources', 'loadCoreComparison'];
const summaryBody = fragment('showResult', '  const report = job.report;', "  $('summary').hidden = false;");
const dropdownBody = fragment('renderGroupOptions', "  $('group-select').replaceChildren();", "  $('group-select').value = highlightedGroup;");
const inspectorBody = fragment('inspect', '    const group = structureGroups.find', '    const wasClosed');
const inspectorTitle = realFunction('inspect').split('\n').find(line => line.includes("$('selected-id').textContent ="));
const outlineBody = fragment('colorMap', "    const unmeasured = point.measurement_status", '    let fill');
vm.runInContext(names.map(realFunction).join('\n') + `\n
function showSummary() {${summaryBody}}
function showDropdown() {${dropdownBody}}
function inspectorGroup(row) {${inspectorBody}}
function inspectorTitle(row) {${inspectorTitle}}
function outlineClass(point, mode, el) {${outlineBody}}
${realFunction('renderColorLegend').replace('function renderColorLegend(', 'function actualColorLegend(')}
${realFunction('inspect').replace('async function inspect(', 'async function runInspection(')}
${realFunction('table').replace('function table(', 'function renderTable(')}
`, context);
const summary = () => $('summary').children.map(el => [el.children[1].textContent, el.children[0].textContent]);
const legend = () => $('legend').children.find(el => el.className === 'group-legend').textContent;
const scenarios = {
  async neighbors() {
    context.showQuerySelection({molecule: query, query_neighbors: [{rank: 1, measured_id: measured.molecule_identity_key, tanimoto: .8}]});
    const row = $('query-neighbor-table').children[0];
    assert.equal(row.data.ID, 'CHEMBL538536 (parent CHEMBL1189346)');
    assert.equal($('query-summary').children[0].children[1].textContent, 'CHEMBL538536 (parent CHEMBL1189346) · Tanimoto 0.8');
    row.onclick(); assert.equal(calls.at(-1), 'id=CHEMBL1189346');
    row.onkeydown({key: 'Enter'}); assert.equal(calls.at(-1), 'id=CHEMBL1189346');
    context.inspectorTitle(measured);
    assert.equal(row.data.ID.split(' (parent ')[0], $('selected-id').textContent);
    assert.equal(context.measuredLabel('Q1'), 'Q1'); assert.equal(context.measuredLabel('unknown'), 'unknown');
  },
  async summary() {
    context.showSummary();
    assert.deepEqual(summary(), [['Input records', '1,100'], ['Mapped molecules', '889'], ['Retained records', '1,100'], ['Wide record ranges', '17']]);
    context.job.report.query_molecules = {mapped: 45}; context.job.molecules = 934; context.showSummary();
    assert.deepEqual(summary(), [['Input records', '1,100'], ['Reference molecules', '889'], ['New molecules', '45'], ['Retained records', '1,100'], ['Wide record ranges', '17']]);
    context.job.report.query_molecules.mapped = 0; context.showSummary();
    assert.ok(!summary().some(([label]) => label === 'New molecules'));
    delete context.job.report.query_molecules; context.job.options.input_format = 'chembl'; context.showSummary();
    assert.equal(summary().length, 4);
  },
  async groups() {
    context.points.push({...measured, structure_group: 'G002'}); await context.selectGroup('G161');
    assert.equal($('group-core-title').textContent, 'G161 · 1 reference member + 1 new');
    assert.equal($('group-core-coverage').textContent, 'Core in 1 of 1 reference member · exact scaffold in 1');
    assert.equal(legend(), 'Outlined: G161, 2 molecules (1 new). Others faded.');
    context.inspectorGroup(query);
    assert.equal($('selected-group-summary').textContent, 'G161 | 1 member | core covers 1/1');
    context.structureGroups = [{...group, size: 2, core_coverage_count: 2, core_exact_scaffold_count: 2}];
    context.points = [measured, {...measured}, ...Array.from({length: 8}, () => ({...query}))];
    await context.selectGroup('G161');
    assert.equal($('group-core-title').textContent, 'G161 · 2 reference members + 8 new');
    assert.equal($('group-core-coverage').textContent, 'Core in 2 of 2 reference members · exact scaffold in 2');
    context.structureGroups[0].core_smiles = ''; await context.selectGroup('G161');
    assert.equal($('group-core-coverage').textContent, '2 reference members · no ring scaffold');
  },
  async plural() {
    for (const count of [0, 1, 3]) assert.equal(context.countLabel(count, 'member'), `${count} member${count === 1 ? '' : 's'}`);
    for (const size of [1, 3]) {
      context.structureGroups = [{...group, size, core_smiles: ''}];
      context.points = Array.from({length: size}, () => ({...measured})); context.showDropdown();
      const members = `${size} member${size === 1 ? '' : 's'}`;
      assert.equal($('group-select').children[1].textContent, `G161 | ${members} | core 1/${size}`);
      await context.selectGroup('G161');
      assert.equal($('group-core-title').textContent, `G161 · ${members}`);
      assert.equal($('group-core-coverage').textContent, `${members} · no ring scaffold`);
      assert.equal(legend(), `Outlined: G161, ${size} molecule${size === 1 ? '' : 's'}. Others faded.`);
      context.inspectorGroup(measured);
      assert.equal($('selected-group-summary').textContent, `G161 | ${members} | core covers 1/${size}`);
    }
  },
  async outlines() {
    for (const point of [measured, query]) for (const mode of ['activity', 'support', 'reliability', 'structure']) {
      const el = node(), applied = {};
      el.classList.toggle = (name, active) => { applied[name] = !!active; };
      context.outlineClass(point, mode, el);
      assert.equal(applied['unmeasured-outline'], point === query && mode === 'activity');
      assert.equal(applied['new-molecule-outline'], point === query && mode !== 'activity');
      assert.equal(applied['reference-evidence-faded'], false);
    }
    for (const scale of [.005, .05, .2, 1, 4, 15]) {
      context.transform = {s: scale, x: 0, y: 0}; context.applyTransform();
      assert.ok(Math.abs(1.5 * Number($('map-stage').style['--zoom-inverse']) * scale - 1.5) < 1e-12);
    }
  },
  async support_legend() {
    context.job.options = {input_format: 'molecule', compound_column: 'compound_id', lower_threshold: 6, upper_threshold: 7, offset_display: 'none'};
    context.points = [{...measured, offset_pair_touched: true}];
    for (const mode of ['activity', 'support']) {
      context.actualColorLegend(mode);
      assert.ok($('legend').children.some(node => node.textContent === '◥ Wide record range'));
    }
    assert.ok($('legend').children.some(node => node.textContent === '◣ offset document pair (display unchanged)'));
  },
  async input_reload() {
    const file = {name: 'same.csv'}, input = {filename: file.name, rows: 2,
      columns: ['compound_id', 'smiles', 'activity'], targets: [], preview: [],
      detected: {molecule_id: 'compound_id', smiles: 'smiles', pchembl: 'activity'}};
    context.getJSON = async () => input;
    context.setBusy = value => { context.busy = value; };
    context.resetResult = () => { context.job = null; };
    context.showQueryInput = () => {};
    context.updateInputType = () => {};
    context.showPage = page => { context.page = page; };
    context.fail = error => { calls.push(error.message); };
    $('csv').files = [file]; $('query-csv').value = 'query.csv';
    for (let load = 0; load < 2; load++) {
      $('csv').value = file.name; $('value-col').value = 'previous-choice'; context.page = 'import';
      await context.loadInput();
      assert.equal($('csv').value, ''); assert.equal(context.page, 'configure');
      assert.equal($('value-col').value, 'activity'); assert.equal($('compound-column').value, 'compound_id');
    }
    context.getJSON = async () => { throw Error('invalid CSV'); };
    $('csv').value = file.name; await context.loadInput();
    assert.equal($('csv').value, ''); assert.equal(calls.at(-1), 'invalid CSV');
    assert.equal($('query-csv').value, 'query.csv');
    context.getJSON = async (path, options) => { assert.equal(options.body, file); return input; };
    $('csv').files = []; await context.loadInput(false, file);
    assert.equal(context.page, 'configure');
  },
  async audit() {
    const pairs = [
      {document_a: 'R1', document_b: 'R3', median_offset: 1.18, shared_entries: 5, sign_agreement: 1, systematic_offset: true, touched_entry_count: 24},
      {document_a: 'R2', document_b: 'R3', median_offset: 1.1, shared_entries: 3, sign_agreement: 1, systematic_offset: true, touched_entry_count: 28},
      {document_a: 'R1', document_b: 'R2', median_offset: -.03, shared_entries: 13, sign_agreement: .62, systematic_offset: false, touched_entry_count: 33}
    ];
    const report = {generic_activity: {input_rows: 134, excluded_activity_records: 30, document_column: 'run_id',
      excluded_reasons: {non_exact_relation: 28, invalid_activity_value: 2, invalid_smiles: 0}},
      activity_confidence: {document_pair_count: 3, observed_document_pair_count: 8, document_pair_possible_count: 10,
        discordant_document_pair_count: 2, offset_min_shared: 3, offset_threshold: 1, offset_sign_agreement: .8,
        offset_display: 'none', document_offsets_top10: pairs}};
    context.job.options = {input_format: 'molecule', compound_column: 'compound_id', document_column: 'run_id'};
    context.showDataAudit(report);
    assert.equal($('audit-body').children.find(node => node.dataset.excludedRecords).textContent,
      'Excluded records: 30 of 134. Censored or inexact values (relation other than =): 28; empty or non-numeric activity values: 2.');
    const offsetSection = () => $('audit-body').children.find(node => node.children[0]?.textContent === 'Evidence units and offsets');
    const section = offsetSection();
    assert.match(section.children[1].textContent, /run_id column; 3 pairs with at least 3 shared entries were evaluated, 2 pairs were flagged/);
    assert.match(section.children[1].textContent, /\|median difference\| ≥ 1 log units and sign agreement ≥ 80%/);
    const rows = section.children[2].children;
    assert.deepEqual(Object.keys(rows[0].data), ['Pair', 'Shared entries', 'Median difference (B − A)', 'Sign agreement', 'Flagged', 'Entries touched']);
    assert.deepEqual(rows.map(node => Object.values(node.data)), [
      ['R1 vs R3', 5, '+1.18', '100%', 'Yes', 24], ['R2 vs R3', 3, '+1.10', '100%', 'Yes', 28],
      ['R1 vs R2', 13, '−0.03', '62%', 'No', 33]]);
    assert.deepEqual(rows.map(node => node.dataset.offsetFlagged), ['true', 'true', 'false']);
    assert.equal(section.children[3].textContent, 'Display unchanged: entries with records from a flagged pair are marked ◣ on Activity class support.');
    assert.equal(section.children[4].textContent, 'Full table: map_document_offsets.csv (Export).');
    report.activity_confidence.offset_display = 'uncertain'; context.showDataAudit(report);
    assert.equal(offsetSection().children[3].textContent, 'Entries with records from a flagged pair are shown as uncertain.');
    delete report.generic_activity.excluded_reasons; report.generic_activity.excluded_non_exact_relation_records = 28;
    context.showDataAudit(report);
    assert.equal($('audit-body').children.find(node => node.dataset.excludedRecords).textContent,
      'Excluded records: 30 of 134. Censored or inexact values (relation other than =): 28.');
    context.job.options.document_column = null; context.job.options.assumed_sigma = .64;
    report.activity_confidence.document_offset_detection = 'not applicable'; context.showDataAudit(report);
    assert.equal($('audit-body').children.find(node => node.dataset.assumedSigma).textContent,
      'σ (assumed) 0.64 log units. Evidence units: one per compound. Offset detection needs an evidence-unit column, so it was not run.');
    assert.equal(offsetSection(), undefined);
    context.job.options = {input_format: 'chembl'}; delete report.generic_activity;
    delete report.activity_confidence.document_offset_detection; report.exclusion_counts = {invalid_smiles: 2};
    context.showDataAudit(report);
    assert.match(offsetSection().children[1].textContent, /Evidence units come from document;/);
    assert.ok($('audit-body').children.some(node => node.textContent === 'Exclusion reasons (may overlap): invalid_smiles: 2'));
  },
  async generic_inspector() {
    const explanation = 'Generic activity records: the listed rows are the retained records for this compound; values were converted to pActivity.';
    const recordColumns = ['source_row', 'compound_id', 'run_id', 'ic50', 'unit', 'relation', 'activity_pchembl_raw', '2', '1', 'smiles'];
    const record = {source_row: 0, compound_id: 'LAB-001', run_id: 'R1', ic50: 10, unit: 'nM', relation: '=', activity_pchembl_raw: 8, '1': 'a', '2': 'b', smiles: 'CCO'};
    const detail = {molecule: {molecule_id: 'LAB-001', activity_pchembl: 8, pchembl_min: 8, pchembl_max: 8, n_records: 101, activity_class: 'active'},
      summary: {input_format: 'generic_activity', n_records: 101, median_verified: true, record_count_verified: true, interpretation: explanation,
        warnings: ['Offset document pair: R1 vs R3; same text as the ChEMBL branch.'], neighborhood: {}},
      neighbors: [], records: Array.from({length: 101}, () => record), record_columns: recordColumns, total_records: 101};
    context.job.options = {input_format: 'molecule', compound_column: 'compound_id', value_col: 'ic50'};
    context.getJSON = async () => detail; await context.runInspection('id=LAB-001');
    const fields = $('evidence-summary').children[0].children.map(node => node.textContent);
    assert.deepEqual(fields, ['Median pActivity', '8', 'pActivity range', '8–8', 'Source records', '101',
      'Activity class', 'active', 'Source-record audit', 'Median recomputed from 101 listed records']);
    assert.ok($('evidence-summary').children.some(node => node.textContent === explanation));
    assert.ok(!$('evidence-summary').children.some(node => node.textContent.includes('Quality flags not verified')));
    assert.ok($('evidence-summary').children.some(node => node.dataset.documentOffsetWarning));
    assert.equal($('record-section').hidden, false); assert.equal($('view-records').hidden, false); assert.equal($('download-evidence').hidden, false);
    assert.equal($('record-count').textContent, '101 records'); assert.equal($('record-table').children.length, 101);
    assert.deepEqual($('record-table').columns, recordColumns);
    const container = node(); context.renderTable([record], container, 999, recordColumns);
    assert.deepEqual(container.children[0].children[0].children[0].children.map(node => node.textContent), recordColumns);
    detail.summary = {median_verified: false, neighborhood: {}, interpretation: 'Molecule-table input: no assay-record audit is available.'};
    detail.records = []; detail.total_records = 0; context.job.options.compound_column = null;
    await context.runInspection('id=LAB-001');
    assert.equal($('record-section').hidden, true); assert.equal($('record-table').children.length, 0);
    assert.ok($('evidence-summary').children.some(node => node.textContent === detail.summary.interpretation));
  },
};
assert.ok(process.argv[3] in scenarios, 'Unknown scenario');
scenarios[process.argv[3]]().catch(error => { console.error(error); process.exitCode = 1; });
