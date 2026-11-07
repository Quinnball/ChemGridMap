'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {runtime, realFunction} = require('./u4_web_contract.cjs');

async function setup(names) {
  const env = await runtime();
  const $ = id => env.document.getElementById(id);
  Object.assign(env.context, {$, URL, busy: false, queryInputInfo: null,
    job: {options: {input_format: 'chembl'}}, updateNavigation() {}, densityPreview() {},
    fetch() { throw Error('Sources must not make network requests'); }});
  vm.runInContext(names.map(realFunction).join('\n'), env.context);
  return {...env, $};
}

const scenarios = {
  async sources() {
    const {context, $, api} = await setup(['showSources']);
    const samples = Array.from({length: 27}, (_, i) => ({
      id: 'CHEMBL' + (1124628 + i), doc_type: 'PUBLICATION', journal: 'J. Med. Chem.',
      year: 1989 + i, n_records: i + 1,
      links: [{kind: 'doi', label: 'DOI 10.1021/jm00131a011', url: 'https://doi.org/10.1021/jm00131a011'}],
      references: [],
    }));
    samples.push({id: 'R1', doc_type: '', year: null, n_records: 1, links: [
      {label: 'javascript:alert(1)', url: 'javascript:alert(1)'},
      {label: 'data:text/html,<b>bad</b>', url: 'data:text/html,<b>bad</b>'},
      {label: 'file:///secret', url: 'file:///secret'}], references: [
      {value: 'javascript:alert(1)', url: ''}, {value: 'data:text/html,<b>bad</b>', url: ''}]
    });
    for (const language of ['en', 'zh']) {
      api.setLanguage(language); context.showSources(samples);
      assert.equal($('sources').hidden, false);
      assert.equal($('sources-title').textContent, language === 'en' ? 'Sources' : '来源文献');
      assert.equal($('source-list').children.length, 28, 'Long lists retain all sources');
      assert.deepEqual($('source-list').children.map(row => row.dataset.sourceId), samples.map(row => row.id));
      const anchors = $('sources').querySelectorAll('a'); assert.equal(anchors.length, 27);
      for (const anchor of anchors) {
        assert.equal(anchor.getAttribute('target'), '_blank');
        assert.equal(anchor.getAttribute('rel'), 'noopener noreferrer');
        assert.match(anchor.getAttribute('href'), /^https:\/\/doi\.org\//);
        assert.match(anchor.getAttribute('aria-label'), /DOI 10\.1021\/jm00131a011/);
      }
      const last = $('source-list').children.at(-1);
      assert.equal(last.querySelectorAll('a').length, 0);
      assert.match(last.textContent, /javascript:alert\(1\)/);
      assert.match(last.textContent, /data:text\/html,<b>bad<\/b>/);
      assert.equal(last.querySelectorAll('b').length, 0, 'Reference values remain text');
      assert.match(last.textContent, language === 'en' ? /No available links/ : /无可用链接/);
      assert.match($('source-list').children[0].textContent, language === 'en' ? /Publication.*Retained records: 1/ : /论文.*保留记录：1 条/);
    }
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      context.showSources([{id: 'CHEMBL-PATENT', doc_type: 'PATENT', n_records: 2, links: [
        {kind: 'patent', label: 'Patent US-10174040-B2', url: 'https://patents.google.com/patent/US10174040B2'}]}]);
      assert.equal($('source-list').querySelector('a').textContent, language === 'en' ? 'Patent US-10174040-B2' : '专利 US-10174040-B2');
    }
    const siblings = $('measured-selection').children;
    assert.ok(siblings.indexOf($('evidence-summary')) < siblings.indexOf($('sources')));
    assert.ok(siblings.indexOf($('sources')) < siblings.indexOf($('view-records')));
    const css = fs.readFileSync(path.join(process.argv[2], 'app.css'), 'utf8');
    assert.match(css, /\.source-list\s*\{[^}]*max-height:\s*300px;[^}]*overflow-y:\s*auto;/);
    assert.equal($('source-list').getAttribute('tabindex'), '0');
    context.showSources(); assert.equal($('sources').hidden, true); assert.equal($('source-list').children.length, 0);
  },
  async references() {
    const {context, $, api, document} = await setup(['selectColumns', 'showInput', 'table', 'stat', 'updateInputType', 'updateEvidenceOptions', 'options', 'restoreOptions']);
    context.updateGrouping = () => {};
    context.showInput({filename: 'records.csv', rows: 2, targets: [], columns: ['compound', 'smiles', 'activity', 'reference'],
      detected: {molecule_id: 'compound', smiles: 'smiles', pchembl: 'activity'}, preview: []});
    assert.equal($('input-format').value, 'generic');
    assert.equal($('reference-column').value, '');
    assert.equal(Object.hasOwn(context.options(), 'reference_column'), false);
    $('reference-column').value = 'reference';
    assert.equal(context.options().reference_column, 'reference');
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      const control = $('reference-column-control');
      assert.match(control.textContent, language === 'en' ? /Reference column/ : /参考文献列/);
      const help = control.querySelector('[data-help]');
      help.dispatchEvent({type: 'pointerenter', pointerType: 'mouse'});
      assert.match($('parameter-tooltip').textContent, language === 'en' ? /do not change evidence units, σ or support/ : /不改变证据单元、σ 或支持度/);
      assert.match(help.getAttribute('aria-label'), language === 'en' ? /^Help: Reference column/ : /^说明：参考文献列/);
      assert.equal($('reference-column').value, 'reference');
    }
    const saved = context.options();
    $('reference-column').value = ''; context.restoreOptions(saved);
    assert.equal($('reference-column').value, 'reference');
    delete saved.reference_column; context.restoreOptions(saved);
    assert.equal($('reference-column').value, '');
    for (const inputType of ['chembl', 'molecule']) {
      $('input-format').value = inputType; $('reference-column').value = 'reference'; context.updateInputType();
      assert.equal($('reference-column-control').hidden, true);
      assert.equal($('reference-column').disabled, true);
      assert.equal(Object.hasOwn(context.options(), 'reference_column'), false);
    }
  },
  async coverage() {
    const {context, $, api} = await setup(['showDataAudit']);
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      context.showDataAudit({document_coverage: {kind: 'chembl', total: 100, doi_or_pubmed: 93, patent: 6, chembl_only: 1}, warnings: []});
      const text = $('audit-body').children.find(item => item.dataset.documentCoverage).textContent;
      assert.match(text, language === 'en' ? /Of 100 documents, 93 have a DOI or PubMed ID, 6 are patents, and 1 have only a ChEMBL document page/ : /100 份文献中，93 份有 DOI 或 PubMed ID，6 份为专利，1 份只有 ChEMBL 文献页/);
      context.showDataAudit({document_coverage: {kind: 'reference', total: 4, linked: 3}, warnings: []});
      assert.match($('audit-body').children.find(item => item.dataset.documentCoverage).textContent,
        language === 'en' ? /3 of 4 evidence units have links/ : /4 个证据单元有 3 个带有链接/);
    }
  },
  async punctuation() {
    const {context, $, api} = await setup(['showDataAudit', 'table', 'format']);
    api.setLanguage('zh');
    assert.equal(api.t('s.032'), '准备好探索分子空间了吗？');
    assert.equal(api.t('d.Medium ({lower}, {upper})', {lower: 6, upper: 7}), '中等活性 (6, 7)');
    assert.match(api.t('d.Core in {covered} of {members} · exact scaffold in {exact}', {covered: 1, members: '1 个已测成员', exact: 1}), /核心骨架覆盖 1 个已测成员中的 1 个/);
    context.job.options = {input_format: 'molecule', compound_column: 'compound', document_column: 'run_id'};
    const report = {query_molecules: {mapped: 45, invalid_smiles: 0, duplicate_query: 0, in_reference: 0},
      generic_activity: {input_rows: 134, excluded_activity_records: 30, document_column: 'run_id', excluded_reasons: {non_exact_relation: 28, invalid_activity_value: 2}},
      activity_confidence: {document_pair_count: 1, offset_min_shared: 3, discordant_document_pair_count: 0, offset_threshold: 1, offset_sign_agreement: .8, document_offsets_top10: []}, warnings: []};
    context.showDataAudit(report);
    const text = $('audit-body').textContent;
    assert.match(text, /已映射：45/);
    assert.match(text, /共 134 条。删失值或非精确值（关系符不是 =）：28；空值或非数值的活性值：2。/);
    assert.match(text, /证据单元来自 run_id 列；/);
    const message = 'This entry has records from R1, R3; the display is unchanged because the offset rule is off.';
    assert.match(api.translateServer(message), /包含来自 R1、R3 的记录/);
    api.setLanguage('en'); assert.equal(api.translateServer(message), message);
    context.showDataAudit(report);
    assert.match($('audit-body').textContent, /mapped: 45/);
    assert.match($('audit-body').textContent, /of 134\. Censored or inexact values \(relation other than =\): 28; empty or non-numeric activity values: 2\./);
  },
};

(async () => {
  const scenario = process.argv[3]; assert.ok(scenario in scenarios);
  await scenarios[scenario](); console.log(`U5 ${scenario}: passed`);
})().catch(error => { console.error(error); process.exitCode = 1; });
