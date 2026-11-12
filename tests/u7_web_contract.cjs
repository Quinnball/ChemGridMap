'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {runtime, realFunction} = require('./u4_web_contract.cjs');

const lookupFile = 'CA II export.csv';
const source = {id: 'CHEMBL1124628', doc_type: 'PUBLICATION', journal: 'J. Med. Chem.', year: 1989,
  n_records: 3, same_values: ['7.20', '8.00'], references: [], links: [
    {kind: 'doi', label: 'DOI 10.1021/jm00131a011', url: 'https://doi.org/10.1021/jm00131a011'},
    {kind: 'pubmed', label: 'PubMed 1', url: 'https://pubmed.ncbi.nlm.nih.gov/1/'},
    {kind: 'chembl', label: 'ChEMBL', url: 'https://www.ebi.ac.uk/chembl/document_report_card/CHEMBL1124628/'},
    {kind: 'patent', label: 'Patent US-1', url: 'https://patents.google.com/patent/US1'},
    {label: 'Bad script', url: 'javascript:alert(1)'}, {label: 'Bad data', url: 'data:text/html,test'}]};

function detail(extra = {}) {
  return {molecule: {molecule_id: 'CMP-0001', grid_row: 0, grid_col: 1}, sources: [],
    source_lookup_file: lookupFile, possible_sources: [source],
    summary: {input_format: 'generic_activity', median_verified: true, n_records: 1,
      interpretation: '', neighborhood: {representation: 'morgan', representation_metric: 'Tanimoto',
        effective_k: 1, shared_count: 1, representation_count: 1, grid_count: 1, recall: 1}},
    neighbors: [], records: [], total_records: 0, svg: '', neighborhood_svg: '', ...extra};
}

async function setup() {
  const env = await runtime(), requests = [];
  const $ = id => { const item = env.document.getElementById(id); assert.ok(item, `Missing element: ${id}`); return item; };
  Object.assign(env.context, {$, URL, t: env.api.t, translateServer: env.api.translateServer,
    job: {options: {input_format: 'molecule', compound_column: 'compound', representation: 'morgan'}, report: {}},
    points: [], structureGroups: [], highlightedGroup: '', currentCell: null, requestNumber: 0,
    busy: false, inputInfo: null, queryInputInfo: null, sourceInputInfo: null, umapAvailable: null,
    highlight() {}, fit() {}, updateNavigation() {}, applyTransform() {}, translateSVGText() {}, refreshValidation() {},
    showResult: () => env.context.showDataAudit(env.context.job.report),
    getJSON: async (url, options) => { requests.push({url, options}); return detail(); },
    fetch() { throw Error('U7 must not make network requests'); },
  });
  const names = ['format', 'ordinal', 'percentage', 'countLabel', 'table', 'stat', 'measuredLabel',
    'showQuerySelection', 'showSources', 'referenceMapDetails', 'renderReferenceMap', 'newMoleculeWarning', 'newMoleculeInsertedMessage', 'showDataAudit', 'inspect', 'status', 'fail', 'refreshLanguage',
    'updateInputType', 'updateEvidenceOptions', 'showSourceInput', 'loadSourceInput', 'removeSourceInput',
    'loadInput', 'showInput', 'selectColumns', 'showQueryInput', 'options'];
  vm.runInContext(names.map(realFunction).join('\n'), env.context);
  env.context.showResult.grouping = {};
  return {...env, $, requests};
}

const scenarios = {
  async sources() {
    const {context, $, api} = await setup();
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); await context.inspect(null, detail());
      assert.equal($('possible-sources').hidden, false);
      assert.equal($('possible-sources-title').textContent, language === 'en'
        ? 'Possible sources (same structure in the ChEMBL file)' : '可能的来源文献（ChEMBL 文件中结构相同）');
      assert.equal($('possible-sources-note').textContent, api.t('u7.sourceNote', {file: lookupFile}));
      assert.match($('possible-source-list').textContent, language === 'en' ? /Publication.*J\. Med\. Chem\..*1989.*ChEMBL records: 3/ : /论文.*J\. Med\. Chem\..*1989.*ChEMBL 记录：3 条/);
      assert.equal($('possible-source-list').querySelector('.source-same-value').textContent,
        language === 'en' ? 'Same value (pChEMBL 7.20, 8.00)' : '数值相同（pChEMBL 7.20、8.00）');
      assert.equal($('possible-sources-empty').hidden, true);
      const anchors = $('possible-source-list').querySelectorAll('a'); assert.equal(anchors.length, 4);
      for (const anchor of anchors) {
        assert.equal(anchor.getAttribute('target'), '_blank');
        assert.equal(anchor.getAttribute('rel'), 'noopener noreferrer');
        assert.match(anchor.getAttribute('href'), /^https:\/\//);
      }
      assert.doesNotMatch($('possible-source-list').textContent, /Bad script|Bad data/);
      await context.inspect(null, detail({possible_sources: [{...source, same_values: []}]}));
      assert.equal($('possible-source-list').querySelector('.source-same-value'), null);
      await context.inspect(null, detail({possible_sources: []}));
      assert.equal($('possible-sources').hidden, false); assert.equal($('possible-sources-empty').hidden, false);
      assert.equal($('possible-sources-empty').textContent, language === 'en'
        ? `No record with the same structure in ${lookupFile}.` : `${lookupFile} 中没有结构相同的记录。`);
      const own = detail({sources: [source]}); delete own.source_lookup_file; delete own.possible_sources;
      await context.inspect(null, own);
      assert.equal($('possible-sources').hidden, true); assert.equal($('sources').hidden, false);
      await context.inspect(null, detail({molecule: {molecule_id: 'Q1', measurement_status: 'not_measured'}}));
      assert.equal($('possible-sources').hidden, true);
    }
    const hostile = '<img src=x onerror=alert(1)>.csv';
    context.showSources([], {source_lookup_file: hostile, possible_sources: []});
    assert.match($('possible-sources-note').textContent, /<img src=x/);
    assert.equal($('possible-sources').querySelectorAll('img').length, 0);
    context.showSources([], {source_lookup_file: lookupFile,
      possible_sources: Array.from({length: 27}, (_, i) => ({...source, id: 'CHEMBL' + i}))});
    assert.equal($('possible-source-list').children.length, 27);
  },
  async audit() {
    const {context, $, api} = await setup();
    const lookup = {file: lookupFile, found: 889, total: 889, same: 881, target_count: 2};
    const report = {document_coverage: {kind: 'reference', total: 2, linked: 0}, source_lookup: lookup};
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); context.showDataAudit(report);
      const children = $('audit-body').children;
      const coverageIndex = children.findIndex(row => row.dataset.documentCoverage);
      assert.equal(children[coverageIndex + 1].dataset.sourceLookup, 'true');
      assert.equal(children[coverageIndex + 1].textContent, language === 'en'
        ? `Source lookup (${lookupFile}): 889 of 889 molecules have the same structure in the file; 881 of them match a reported pChEMBL.`
        : `来源文献查找（${lookupFile}）：889 个分子中，889 个在该文件中有结构相同的记录，其中 881 个与某篇文献报告的 pChEMBL 相同。`);
      assert.equal(children[coverageIndex + 2].textContent, language === 'en'
        ? 'The source lookup file contains 2 targets; possible sources may come from other targets.'
        : '来源文献查找文件包含 2 个靶标，可能的来源文献也可能来自其他靶标。');
      context.showDataAudit({...report, source_lookup: {...lookup, target_count: 1}});
      assert.equal($('audit-body').children.some(row => row.dataset.sourceLookupTargets), false);
      context.showDataAudit({}); assert.equal($('audit-body').children.some(row => row.dataset.sourceLookup), false);
    }
  },
  async controls() {
    const {context, $, api, requests, document} = await setup();
    const input = {filename: lookupFile, rows: 2281, columns: ['canonical_smiles', 'document_chembl_id']};
    const measured = {filename: 'records.csv', rows: 1, targets: [], columns: ['compound', 'smiles', 'activity'],
      detected: {molecule_id: 'compound', smiles: 'smiles', pchembl: 'activity'}, preview: []};
    const file = {name: lookupFile};
    let resets = 0;
    context.setBusy = value => { context.busy = value; };
    context.resetResult = () => { resets++; context.job = null; };
    context.showPage = value => { context.page = value; };
    context.getJSON = async (url, options) => { requests.push({url, options}); return url === '/api/input' ? measured : input; };
    assert.equal(document.querySelectorAll('[data-help]').length, 77);
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      $('input-format').value = 'generic'; context.updateInputType();
      assert.equal($('source-input-control').hidden, false); assert.equal($('source-csv').disabled, false);
      const help = $('source-input-control').querySelector('[data-help]');
      help.dispatchEvent({type: 'pointerenter', pointerType: 'mouse'});
      assert.match($('parameter-tooltip').textContent, language === 'en' ? /default is not to use.*map and exported files do not change.*not uploaded.*few seconds/ : /默认不使用.*地图与导出文件不变.*不会上传.*几秒/);
      assert.match(help.getAttribute('aria-label'), language === 'en' ? /Source lookup: ChEMBL export/ : /来源文献查找：ChEMBL 导出文件/);
      $('source-csv').files = [file]; $('source-csv').value = lookupFile;
      await context.loadSourceInput();
      assert.equal(requests.at(-1).url, '/api/source-input'); assert.equal(requests.at(-1).options.body, file);
      assert.equal(requests.at(-1).options.headers['X-Filename'], encodeURIComponent(lookupFile));
      assert.equal(context.job, null); assert.equal(context.page, 'configure');
      assert.equal($('source-csv').value, ''); assert.equal($('source-input-name').textContent, lookupFile);
      assert.equal($('source-input-rows').textContent, api.t('d.{count} rows', {count: '2,281'}));
      assert.equal($('source-input-info').hidden, false);
      assert.equal(Object.keys(context.options()).some(key => key.includes('source_lookup') || key === 'source_input'), false);
      for (const type of ['chembl', 'molecule']) {
        $('input-format').value = type; context.updateInputType();
        assert.equal($('source-input-control').hidden, true); assert.equal($('source-csv').disabled, true);
      }
      $('input-format').value = 'generic'; context.updateInputType(); await context.removeSourceInput();
      assert.equal(requests.at(-1).url, '/api/source-input/clear');
      assert.equal($('source-input-info').hidden, true); assert.equal(context.sourceInputInfo, null);
      context.showSourceInput(input); $('csv').files = [{name: 'records.csv'}];
      await context.loadInput(); assert.equal(context.sourceInputInfo, null);
      assert.equal($('source-input-info').hidden, true);
    }
    assert.equal(resets, 6);
  },
  async errors() {
    const {context, $, api} = await setup();
    const message = 'The source lookup file needs SMILES and ChEMBL document ID columns, as in a ChEMBL activity export.';
    context.setBusy = value => { context.busy = value; };
    context.getJSON = async () => { throw Error(message); };
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); $('source-csv').files = [{name: 'bad.csv'}];
      await context.loadSourceInput();
      assert.equal($('status').textContent, language === 'en' ? message
        : '来源文献查找文件需要 SMILES 列与 ChEMBL 文献编号列，即 ChEMBL 活性导出文件的格式。');
      assert.equal($('status').className, 'error'); assert.equal(context.busy, false);
      context.status('Matching structures to the source lookup file...', 'running');
      assert.equal($('build-progress').textContent, language === 'en'
        ? 'Matching structures to the source lookup file...' : '正在按结构匹配来源文献查找文件……');
    }
  },
  async language() {
    const {context, $, api} = await setup();
    context.showSourceInput({filename: lookupFile, rows: 2281});
    await context.inspect(null, detail());
    context.job.report.source_lookup = {file: lookupFile, total: 2, found: 1, same: 1, target_count: 1};
    context.getJSON = () => { throw Error('Language changes must reuse cached lookup results'); };
    for (const language of ['zh', 'en']) {
      api.setLanguage(language); await context.refreshLanguage();
      assert.equal($('source-input-name').textContent, lookupFile);
      assert.equal($('source-input-rows').textContent, api.t('d.{count} rows', {count: '2,281'}));
      assert.equal($('possible-sources-note').textContent, api.t('u7.sourceNote', {file: lookupFile}));
      assert.match($('audit-body').textContent, language === 'en' ? /Source lookup/ : /来源文献查找/);
    }
  },
  async summary_columns() {
    await setup();
    const css = fs.readFileSync(path.join(process.argv[2], 'app.css'), 'utf8').replace(/\/\*[\s\S]*?\*\//g, '');
    const declarations = selector => [...css.matchAll(/([^{}]+)\{([^{}]*)\}/g)]
      .filter(([, selectors]) => selectors.split(',').some(item => item.trim() === selector))
      .map(([, , body]) => body).join('\n');
    const grid = declarations('#evidence-summary dl');
    assert.match(grid, /display:\s*grid\s*;/);
    assert.match(grid, /grid-template-columns:\s*minmax\(0,\s*1\.25fr\)\s+minmax\(0,\s*1fr\)\s*;/);
    assert.match(grid, /gap:\s*8px\s*;/);
    assert.match(grid, /font-size:\s*11px\s*;/);
    for (const selector of ['#evidence-summary dt', '#evidence-summary dd']) {
      const rule = declarations(selector);
      assert.match(rule, /min-width:\s*0\s*;/);
      assert.match(rule, /overflow-wrap:\s*anywhere\s*;/);
      assert.doesNotMatch(rule, /white-space:\s*nowrap/);
    }
    // VM 没有浏览器排版引擎；这里只约束列宽与换行规则，不证明实际行数或无溢出。
  },
  async median_plural() {
    const {context, $, api} = await setup();
    context.getJSON = () => { throw Error('Language changes must reuse cached records'); };
    for (const count of [0, 1, 2]) {
      const result = detail(); result.summary.n_records = count;
      api.setLanguage('en'); await context.inspect(null, result);
      for (const language of ['en', 'zh', 'en']) {
        api.setLanguage(language); await context.refreshLanguage();
        const fields = $('evidence-summary').children[0].children;
        const auditIndex = fields.findIndex(node => node.tagName === 'DT'
          && node.textContent === (language === 'en' ? 'Source-record audit' : '来源记录审计'));
        assert.ok(auditIndex >= 0);
        assert.equal(fields[auditIndex + 1].textContent, language === 'en'
          ? `Median recomputed from ${count} listed ${count === 1 ? 'record' : 'records'}`
          : `已根据列出的 ${count} 条记录重新计算中位数`);
      }
    }
  },
  async same_value_separator() {
    const {context, $, api} = await setup();
    context.getJSON = () => { throw Error('Language changes must reuse cached sources'); };
    for (const values of [[], ['8.05'], ['8.05', '8.10']]) {
      const originalValues = [...values];
      const result = detail({possible_sources: [{...source, same_values: values}]});
      api.setLanguage('en'); await context.inspect(null, result);
      for (const language of ['en', 'zh', 'en']) {
        api.setLanguage(language); await context.refreshLanguage();
        const same = $('possible-source-list').querySelector('.source-same-value');
        if (!values.length) assert.equal(same, null);
        else assert.equal(same.textContent, language === 'en'
          ? `Same value (pChEMBL ${values.join(', ')})` : `数值相同（pChEMBL ${values.join('、')}）`);
      }
      assert.deepEqual(result.possible_sources[0].same_values, originalValues);
    }
  },
  async control_order() {
    const {context, $, api, document} = await setup();
    const sourceControl = $('source-input-control'), parent = sourceControl.parentElement;
    const relationRow = $('relation-column').closest('.form-row');
    assert.equal($('label-col').closest('.form-row'), relationRow);
    assert.equal(relationRow.parentElement, parent);
    assert.equal($('reference-column-control').parentElement, parent);
    const siblings = parent.children, relationIndex = siblings.indexOf(relationRow);
    assert.equal(siblings[relationIndex - 1], $('reference-column-control'));
    assert.equal(siblings[relationIndex + 1], sourceControl);
    assert.equal(siblings[relationIndex + 2].dataset.i18nText, 's.060');
    assert.equal(document.querySelectorAll('[data-help]').length, 77);
    assert.equal(sourceControl.querySelectorAll('[data-help]').length, 1);
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      for (const type of ['generic', 'chembl', 'molecule']) {
        $('input-format').value = type; context.updateInputType();
        assert.equal(sourceControl.hidden, type !== 'generic');
        assert.equal($('source-csv').disabled, type !== 'generic');
      }
    }
  },
  async u6_layout() {
    const {context, $, api} = await setup();
    const report = {document_coverage: {kind: 'chembl', total: 2, doi_or_pubmed: 2, patent: 0, chembl_only: 0},
      generic_activity: {input_rows: 4, excluded_activity_records: 1, document_column: 'run_id', excluded_reasons: {non_exact_relation: 1}},
      exclusion_counts: {invalid_smiles: 2}, query_molecules: {mapped: 1, in_reference: 0, duplicate_query: 0, invalid_smiles: 0},
      activity_confidence: {document_pair_count: 1, offset_min_shared: 3, discordant_document_pair_count: 0,
        offset_threshold: 1, offset_sign_agreement: .8, document_offsets_top10: []}};
    context.job.options.assumed_sigma = .64; context.job.options.document_column = 'run_id';
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      for (const type of ['chembl', 'molecule']) {
        context.job.options.input_format = type; context.showDataAudit(report);
        const rows = $('audit-body').children;
        const coverage = rows.findIndex(row => row.dataset.documentCoverage);
        const sigma = rows.findIndex(row => row.dataset.assumedSigma);
        const excluded = rows.findIndex(row => row.dataset.excludedRecords);
        const chemblExcluded = rows.findIndex(row => row.textContent.startsWith(api.t('d.Exclusion reasons (may overlap): ')));
        const queries = rows.findIndex(row => row.querySelector('h3')?.textContent === api.t('u8.newMolecules'));
        assert.ok(coverage < sigma && sigma < excluded && excluded < chemblExcluded && chemblExcluded < queries);
        assert.equal(rows[excluded].dataset.excludedRecords, 'true');
        assert.match($('audit-body').textContent, language === 'en'
          ? type === 'chembl' ? /Evidence units come from document;/ : /Evidence units come from run_id column;/
          : type === 'chembl' ? /证据单元来自文献；/ : /证据单元来自 run_id 列；/);
        assert.doesNotMatch($('audit-body').textContent, /证据单元来自 文献/);
      }
    }
  },
};

(async () => {
  const scenario = process.argv[3]; assert.ok(scenario in scenarios);
  await scenarios[scenario](); console.log(`U7 ${scenario}: passed`);
})().catch(error => { console.error(error); process.exitCode = 1; });
