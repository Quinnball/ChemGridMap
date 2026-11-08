'use strict';
const assert = require('node:assert/strict');
const vm = require('node:vm');
const {runtime, realFunction} = require('./u4_web_contract.cjs');

const measured = {molecule_id: 'CHEMBL538536', molecule_identity_key: 'CHEMBL1189346',
  parent_molecule_id: 'CHEMBL1189346', grid_row: 3, grid_col: 4, activity_class: 'active'};
const match = {query_id: 'GEN-0005', status: 'in_reference', reference_id: 'CHEMBL1189346', first_query_id: ''};
const source = {id: 'CHEMBL-DOC', doc_type: 'PUBLICATION', n_records: 2,
  links: [{kind: 'doi', label: 'DOI 10.1000/test', url: 'https://doi.org/10.1000/test'}], references: []};
const auditRow = {query_id: match.query_id, reference_id: match.reference_id,
  molecule_id: measured.molecule_id, parent_molecule_id: measured.parent_molecule_id,
  source_count: 1, linked_source_count: 1};

function detail(molecule = measured, queryMatch = match) {
  return {molecule, query_match: queryMatch,
    summary: {median_verified: true, quality_annotations_verified: true, median_pchembl: 8,
      minimum_pchembl: 8, maximum_pchembl: 8, n_records: 2, n_documents: 1, n_assays: 1,
      neighborhood: {representation: 'morgan', representation_metric: 'Tanimoto', effective_k: 1,
        shared_count: 1, representation_count: 1, grid_count: 1, recall: 1}},
    neighbors: [], query_neighbors: [], records: [], total_records: 0, svg: '', neighborhood_svg: '', sources: [source]};
}

function report(rows = [auditRow]) {
  return {query_molecules: {mapped: 1, invalid_smiles: 0, duplicate_query: 0, in_reference: rows.length},
    query_reference_matches: rows, warnings: []};
}

async function setup() {
  const env = await runtime(), requests = [];
  const $ = id => { const item = env.document.getElementById(id); assert.ok(item, `Missing element: ${id}`); return item; };
  Object.assign(env.context, {$, URL, t: env.api.t, translateServer: env.api.translateServer,
    job: {options: {input_format: 'chembl', representation: 'morgan'}, report: report()},
    points: [measured], structureGroups: [], highlightedGroup: '', currentCell: null,
    requestNumber: 0, busy: false, inputInfo: null, queryInputInfo: null, sourceInputInfo: null, umapAvailable: null,
    highlight() {}, fit() {}, updateInputType() {}, applyTransform() {}, translateSVGText() {}, refreshValidation() {},
    showResult: () => env.context.showDataAudit(env.context.job.report),
    getJSON: async url => { requests.push(url); return detail(); },
    fetch() { throw Error('U6 must not make network requests'); },
  });
  const names = ['format', 'ordinal', 'percentage', 'countLabel', 'table', 'measuredLabel',
    'showQuerySelection', 'showSources', 'showDataAudit', 'inspect', 'status', 'fail', 'refreshLanguage'];
  vm.runInContext(names.map(realFunction).join('\n'), env.context);
  env.context.showResult.grouping = {};
  return {...env, $, requests};
}

const scenarios = {
  async matched() {
    const {context, $, api} = await setup();
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      await context.inspect('id=GEN-0005');
      assert.equal($('query-match-note').hidden, false);
      assert.equal($('selected-id').textContent, measured.molecule_id);
      assert.equal($('query-match-note').textContent, language === 'en'
        ? 'GEN-0005 is not shown separately: it has the same structure as the reference molecule CHEMBL538536 (parent CHEMBL1189346). Its evidence and sources are shown below.'
        : 'GEN-0005 未单独入图：它与参照分子 CHEMBL538536（母体 CHEMBL1189346） 结构相同。下面显示该参照分子的证据与来源文献。');
      assert.equal($('source-list').children.length, 1);
      assert.equal($('sources').hidden, false);
      const noParent = {...measured, molecule_id: match.reference_id};
      await context.inspect(null, detail(noParent));
      assert.doesNotMatch($('query-match-note').textContent, /\(parent|（母体/);
      // 普通 ID、格子和已映射未测分子不保留上一次别名查找的说明。
      for (const queryMatch of [undefined, {...match, status: 'mapped'}]) {
        const ordinary = detail();
        if (queryMatch === undefined) delete ordinary.query_match;
        else ordinary.query_match = queryMatch;
        await context.inspect(null, ordinary);
        assert.equal($('query-match-note').hidden, true);
        assert.equal($('query-match-note').textContent, '');
      }
    }
    const hostile = {...match, query_id: '<img src=x onerror=alert(1)>'};
    await context.inspect(null, detail(measured, hostile));
    assert.match($('query-match-note').textContent, /<img src=x/);
    assert.equal($('query-match-note').querySelectorAll('img').length, 0);
  },
  async duplicates() {
    const {context, $, api} = await setup();
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      const repeated = {query_id: 'GEN-0123', status: 'duplicate_query', reference_id: '', first_query_id: 'GEN-0045'};
      const unmeasured = {...measured, query_id: 'GEN-0045', molecule_id: 'GEN-0045', measurement_status: 'not_measured',
        nearest_measured_id: measured.molecule_identity_key, nearest_measured_similarity: .8};
      await context.inspect(null, detail(unmeasured, repeated));
      assert.equal($('selected-id').textContent, 'GEN-0045');
      assert.equal($('query-match-note').textContent, language === 'en'
        ? 'GEN-0123 repeats GEN-0045 (same structure).'
        : 'GEN-0123 与 GEN-0045 结构相同（重复）。');
      assert.equal($('query-selection').hidden, false);
      assert.equal($('sources').hidden, true);
      repeated.reference_id = measured.molecule_identity_key;
      await context.inspect(null, detail(measured, repeated));
      assert.equal($('selected-id').textContent, measured.molecule_id);
      assert.match($('query-match-note').textContent, language === 'en' ? /GEN-0123 repeats GEN-0045/ : /GEN-0123 与 GEN-0045/);
      assert.match($('query-match-note').textContent, language === 'en'
        ? /GEN-0045 is not shown separately:.*CHEMBL538536 \(parent CHEMBL1189346\)/
        : /GEN-0045 未单独入图：.*CHEMBL538536（母体 CHEMBL1189346）/);
      assert.equal($('query-selection').hidden, true);
      assert.equal($('sources').hidden, false);
    }
  },
  async audit() {
    const {context, $, api, requests} = await setup();
    const rows = Array.from({length: 5}, (_, index) => ({...auditRow, query_id: 'GEN-000' + (index + 1)}));
    context.getJSON = async url => {
      requests.push(url);
      return detail(measured, {...match, query_id: decodeURIComponent(url.split('id=')[1])});
    };
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); context.showDataAudit(report(rows));
      const container = $('audit-body').querySelector('.query-reference-matches');
      const section = container.parentNode, heading = section.children[5];
      assert.equal(section.children.slice(1, 5).filter(item => item.dataset.queryStatus).length, 4);
      assert.equal(heading.textContent, language === 'en' ? 'Same structure as a reference molecule' : '与参照分子结构相同');
      assert.deepEqual(container.querySelectorAll('th').map(item => item.textContent), language === 'en'
        ? ['New molecule ID', 'Matched reference molecule ID', 'Sources', 'Sources with links']
        : ['新分子 ID', '对应的参照分子 ID', '来源文献数', '带链接的来源文献数']);
      const renderedRows = container.querySelectorAll('tbody tr');
      assert.equal(renderedRows.length, 5);
      assert.deepEqual(renderedRows[0].children.map(item => item.textContent),
        [rows[0].query_id, language === 'en' ? 'CHEMBL538536 (parent CHEMBL1189346)' : 'CHEMBL538536（母体 CHEMBL1189346）', '1', '1']);
      for (const [index, row] of renderedRows.entries()) {
        $('audit-dialog').showModal(); row.click(); await new Promise(setImmediate);
        assert.equal(requests.at(-1), '/api/inspect?id=' + rows[index].query_id);
        assert.equal($('audit-dialog').open, false);
        assert.equal($('selected-id').textContent, measured.molecule_id);
        assert.match($('query-match-note').textContent, new RegExp(rows[index].query_id));
        assert.equal($('sources').hidden, false);
      }
      for (const key of ['Enter', ' ']) {
        let prevented = false;
        renderedRows[0].onkeydown({key, preventDefault() { prevented = true; }}); await new Promise(setImmediate);
        assert.equal(prevented, true);
        assert.equal(requests.at(-1), '/api/inspect?id=' + rows[0].query_id);
      }
    }
  },
  async limit() {
    const {context, $, api} = await setup();
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      const rows = Array.from({length: 53}, (_, index) => ({...auditRow, query_id: `GEN-${index}`}));
      context.showDataAudit(report(rows));
      const container = $('audit-body').querySelector('.query-reference-matches');
      assert.equal(container.querySelectorAll('tbody tr').length, 50);
      assert.equal(container.querySelectorAll('tbody tr').at(-1).dataset.queryId, 'GEN-49');
      assert.match(container.parentNode.textContent, language === 'en'
        ? /3 more rows\. Full list: new_molecules_check\.csv \(Export\)/
        : /另有 3 行。完整列表：new_molecules_check\.csv（导出）/);
      context.showDataAudit(report(rows.slice(0, 50)));
      assert.doesNotMatch($('audit-body').textContent, /new_molecules_check\.csv/);
    }
  },
  async no_sources() {
    const {context, $, api} = await setup();
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      context.showDataAudit(report([{...auditRow, source_count: 0, linked_source_count: 0}]));
      const rows = $('audit-body').querySelector('.query-reference-matches').querySelectorAll('tbody tr');
      assert.equal(rows.length, 1);
      assert.deepEqual(rows[0].children.slice(2).map(item => item.textContent), ['0', '0']);
      context.showDataAudit(report([]));
      assert.equal($('audit-body').querySelector('.query-reference-matches'), null);
    }
  },
  async language() {
    const {context, $, api, requests, document} = await setup();
    await context.inspect('id=GEN-0005'); context.showDataAudit(context.job.report);
    assert.equal(requests.length, 1);
    context.getJSON = () => { throw Error('Language changes must reuse the inspected match'); };
    let languageWork;
    document.addEventListener('chemgridmap:languagechange', () => { languageWork = context.refreshLanguage(); });
    for (const language of ['zh', 'en']) {
      api.setLanguage(language); await languageWork;
      assert.match($('query-match-note').textContent, language === 'en' ? /is not shown separately/ : /未单独入图/);
      assert.match($('audit-body').textContent, language === 'en' ? /Same structure as a reference molecule/ : /与参照分子结构相同/);
      assert.equal($('sources-title').textContent, language === 'en' ? 'Sources' : '来源文献');
      assert.equal($('selected-id').textContent, measured.molecule_id);
      assert.equal(context.inspect.lastDetail.query_match.query_id, 'GEN-0005');
      assert.equal(requests.length, 1);
    }
  },
  async invalid() {
    const {context, $, api} = await setup();
    context.getJSON = async () => { throw Error('New molecule GEN-BAD has invalid SMILES.'); };
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); await context.inspect('id=GEN-BAD');
      assert.equal($('status').textContent, language === 'en'
        ? 'New molecule GEN-BAD has invalid SMILES.' : '新分子 GEN-BAD 的 SMILES 无效。');
      assert.equal($('status').className, 'error');
      assert.equal($('status').hidden, false);
    }
  },
};

(async () => {
  const scenario = process.argv[3]; assert.ok(scenario in scenarios);
  await scenarios[scenario](); console.log(`U6 ${scenario}: passed`);
})().catch(error => { console.error(error); process.exitCode = 1; });
