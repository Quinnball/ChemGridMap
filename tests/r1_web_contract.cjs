'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(process.argv[2], 'utf8');
const elements = new Map(), calls = [];
// Only the nodes touched by the tested functions are needed; this is not a DOM emulator.
function node() {
  return {children: [], dataset: {}, style: {setProperty(key, value) { this[key] = value; }},
    classList: {toggle() {}}, textContent: '', hidden: false, value: '',
    replaceChildren(...children) { this.children = children; }, append(...children) { this.children.push(...children); },
    setAttribute() {}, get firstElementChild() { return this.children[0]; },
    querySelectorAll() { return this.children; }};
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
const context = vm.createContext({$, document, points: [measured, query], structureGroups: [group],
  highlightedGroup: '', groupRequestNumber: 0, job: {options: {input_format: 'molecule'}, molecules: 889,
    report: {input_rows: 1100, retained_activity_records: 1100, repeated_conflicted_molecules: 17}},
  transform: {s: 1, x: 0, y: 0}, view: 'grid', tileWidth: 100, hasStructures: false,
  highlight() {}, updateNavigation() {}, coreRule() { return ''; },
  renderColorLegend() { $('legend').replaceChildren(); }, getJSON: async () => ({core_svg: ''}),
  fail(error) { throw error; }, inspect: query => calls.push(query),
  table(rows, container) { container.replaceChildren(...rows.map(data => ({...node(), data}))); },
});
$('color-mode').value = 'activity'; $('measurement-filter').value = 'all';
$('measurement-filter').selectedOptions = [{textContent: 'All molecules'}];
const names = ['stat', 'countLabel', 'measuredLabel', 'showQuerySelection', 'selectGroup',
  'format', 'percentage', 'legendItem', 'renderLegend', 'applyTransform'];
const summaryBody = fragment('showResult', '  const report = job.report;', "  $('summary').hidden = false;");
const dropdownBody = fragment('showResult', "  $('group-select').replaceChildren();", "  $('group-select').value = highlightedGroup;");
const inspectorBody = fragment('inspect', '    const group = structureGroups.find', '    const wasClosed');
const inspectorTitle = realFunction('inspect').split('\n').find(line => line.includes("$('selected-id').textContent ="));
const outlineBody = fragment('colorMap', "    const unmeasured = point.measurement_status", '    let fill');
vm.runInContext(names.map(realFunction).join('\n') + `\n
function showSummary() {${summaryBody}}
function showDropdown() {${dropdownBody}}
function inspectorGroup(row) {${inspectorBody}}
function inspectorTitle(row) {${inspectorTitle}}
function outlineClass(point, mode, el) {${outlineBody}}
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
    assert.deepEqual(summary(), [['Input records', '1,100'], ['Mapped molecules', '934'], ['Unmeasured molecules', '45'], ['Retained records', '1,100'], ['Wide record ranges', '17']]);
    context.job.report.query_molecules.mapped = 0; context.showSummary();
    assert.ok(!summary().some(([label]) => label === 'Unmeasured molecules'));
    delete context.job.report.query_molecules; context.job.options.input_format = 'chembl'; context.showSummary();
    assert.equal(summary().length, 4);
  },
  async groups() {
    context.points.push({...measured, structure_group: 'G002'}); await context.selectGroup('G161');
    assert.equal($('group-core-title').textContent, 'G161 · 1 member · 1 unmeasured');
    assert.equal($('group-core-coverage').textContent, 'Core in 1 of 1 member · exact scaffold in 1');
    assert.equal(legend(), 'Outlined: G161, 2 molecules (1 unmeasured). Others faded.');
    context.inspectorGroup(query);
    assert.equal($('selected-group-summary').textContent, 'G161 | 1 member | core covers 1/1');
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
      const el = node(); let applied;
      el.classList.toggle = (name, active) => { assert.equal(name, 'unmeasured-outline'); applied = active; };
      context.outlineClass(point, mode, el);
      assert.equal(applied, point === query && ['activity', 'support'].includes(mode));
    }
    for (const scale of [.005, .05, .2, 1, 4, 15]) {
      context.transform = {s: scale, x: 0, y: 0}; context.applyTransform();
      assert.ok(Math.abs(1.5 * Number($('map-stage').style['--zoom-inverse']) * scale - 1.5) < 1e-12);
    }
  },
};
assert.ok(process.argv[3] in scenarios, 'Unknown scenario');
scenarios[process.argv[3]]().catch(error => { console.error(error); process.exitCode = 1; });
