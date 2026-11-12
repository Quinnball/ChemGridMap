'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assets = process.argv[2], scenario = process.argv[3];
const html = fs.readFileSync(path.join(assets, 'index.html'), 'utf8');
const source = fs.readFileSync(path.join(assets, 'app.js'), 'utf8');
const dictionaries = JSON.parse(fs.readFileSync(path.join(assets, 'i18n.json'), 'utf8'));
const controls = ['open-previous-results','previous-run-select','open-selected-run','single-smiles','single-id','single-label','single-mode','single-value','single-unit','single-class','single-batch','single-note','measured-result','result-mode','result-value','result-unit','result-class','result-batch','result-note','new-mismatch-filter', 'input-format', 'target', 'activity-type', 'compound-column', 'smiles-col',
  'value-col', 'value-unit', 'unit-column', 'document-column', 'reference-column', 'relation-column', 'label-col',
  'new-activity-col', 'new-unit', 'new-unit-col', 'new-relation-col', 'new-document-col', 'new-label-col', 'new-class-filter', 'new-label-filter', 'source-csv', 'query-csv', 'query-smiles-col', 'query-id-col', 'representation', 'projection', 'occupancy',
  'grouping', 'group-similarity', 'cluster-col', 'cluster-source', 'core-rule', 'core-coverage-target',
  'fmcs-timeout', 'grid-refinement', 'activity-confidence', 'merge-identical-documents',
  'use-assumed-sigma', 'assumed-sigma', 'evidence-model', 'support-level', 'offset-display',
  'offset-min-shared', 'offset-min-difference', 'offset-sign-agreement', 'identity-conflicts',
  'lower', 'upper', 'conflict-range', 'seed', 'identity', 'structure-policy', 'detail', 'raster',
  'color-mode', 'measurement-filter', 'group-select', 'map-core-rule', 'group-core-comparison', 'show-neighbor-links'];

function decoded(text) {
  return text.replace(/&(?:amp|lt|gt|quot|apos|#39|nbsp);/g, value =>
    ({'&amp;': '&', '&lt;': '<', '&gt;': '>', '&quot;': '"', '&apos;': "'", '&#39;': "'", '&nbsp;': '\u00a0'})[value]);
}

// 只实现这些契约实际访问的 DOM 表面；不代替浏览器排版验收。
function environment({saved = null, browserLanguage = 'en-US', blockedStorage = false} = {}) {
  const nodes = [], storage = new Map(saved ? [['chemgridmap.language', saved]] : []), listeners = new Map();
  let document;
  function textNode(text) { return {nodeType: 3, textContent: text, parentNode: null}; }
  function node(tag = 'div') {
    const result = {nodeType: 1, tagName: tag.toUpperCase(), childNodes: [], dataset: {}, attributes: {},
      parentNode: null, style: {setProperty(key, value) { this[key] = value; }}, hidden: false,
      disabled: false, checked: false, open: false, className: '', events: new Map(),
      validity: {valid: true}, setCustomValidity(message) { this.validationMessage = message; },
      get children() { return this.childNodes.filter(child => child.nodeType === 1); },
      get parentElement() { return this.parentNode; },
      get firstElementChild() { return this.children[0]; },
      get textContent() { return this.childNodes.map(child => child.textContent).join(''); },
      set textContent(value) { this.replaceChildren(textNode(String(value))); },
      get innerHTML() { return this.textContent; }, set innerHTML(value) { this.textContent = value; },
      get value() {
        if (this._value != null) return this._value;
        if (this.tagName === 'SELECT') return (this.children.find(item => item.attributes.selected != null) || this.children[0])?.value || '';
        return this.attributes.value || '';
      },
      set value(value) { this._value = String(value); },
      get options() { return this.children; },
      get selectedOptions() { return this.children.filter(item => item.value === this.value); },
      get offsetWidth() { return this.id === 'parameter-tooltip' ? 300 : 16; },
      get offsetHeight() { return this.id === 'parameter-tooltip' ? 140 : 16; },
      append(...items) { for (let item of items) { if (typeof item === 'string') item = textNode(item); item.remove?.(); item.parentNode = this; this.childNodes.push(item); } },
      appendChild(item) { this.append(item); return item; },
      cloneNode(deep=false) { const copy=node(this.tagName.toLowerCase()); for(const [key,value] of Object.entries(this.attributes))copy.setAttribute(key,value); Object.assign(copy.dataset,this.dataset); if(deep)for(const child of this.childNodes)copy.append(child.nodeType===3?textNode(child.textContent):child.cloneNode(true)); return copy; },
      prepend(...items) { for (const item of items.reverse()) { item.parentNode = this; this.childNodes.unshift(item); } },
      replaceChildren(...items) { this.childNodes = []; this.append(...items); },
      remove() { if (this.parentNode) this.parentNode.childNodes = this.parentNode.childNodes.filter(item => item !== this); this.parentNode = null; },
      contains(item) { return item === this || this.childNodes.some(child => child.contains?.(item)); },
      setAttribute(key, value) {
        this.attributes[key] = String(value);
        if (key === 'id' || key === 'class' || key === 'type') this[key === 'class' ? 'className' : key] = String(value);
        if (key.startsWith('data-')) this.dataset[key.slice(5).replace(/-([a-z])/g, (_, char) => char.toUpperCase())] = String(value);
        if (key === 'checked') this.checked = true;
        if (key === 'hidden') this.hidden = true;
      },
      getAttribute(key) { return this.attributes[key] ?? null; },
      hasAttribute(key) { return Object.hasOwn(this.attributes, key); },
      removeAttribute(key) { delete this.attributes[key]; },
      matches(selector) {
        const tagName = selector.match(/^[a-z][a-z0-9-]*/i)?.[0];
        if (tagName && this.tagName !== tagName.toUpperCase()) return false;
        const id = selector.match(/#([\w-]+)/)?.[1]; if (id && this.id !== id) return false;
        for (const [, name] of selector.matchAll(/\.([\w-]+)/g)) if (!this.className.split(/\s+/).includes(name)) return false;
        for (const [, name, expected] of selector.matchAll(/\[([^=\]\s]+)(?:=["']?([^"'\]]*)["']?)?\]/g)) {
          if (!this.hasAttribute(name)) return false;
          if (expected != null && this.getAttribute(name) !== expected) return false;
        }
        return true;
      },
      closest(selector) { for (let item = this; item; item = item.parentNode) if (item.matches?.(selector)) return item; return null; },
      querySelectorAll(selector) {
        const candidates = [];
        function visit(item) { for (const child of item.children) { candidates.push(child); visit(child); } }
        visit(this);
        return candidates.filter(item => selector.split(',').some(part => {
          const sequence = part.trim().split(/\s+/);
          if (!item.matches(sequence.pop())) return false;
          let ancestor = item.parentNode;
          while (sequence.length) {
            const wanted = sequence.pop();
            while (ancestor && !ancestor.matches(wanted)) ancestor = ancestor.parentNode;
            if (!ancestor) return false;
            ancestor = ancestor.parentNode;
          }
          return true;
        }));
      },
      querySelector(selector) { return this.querySelectorAll(selector)[0] || null; },
      addEventListener(type, callback) { if (!this.events.has(type)) this.events.set(type, []); this.events.get(type).push(callback); },
      dispatchEvent(event) {
        event.target ||= this; event.currentTarget = this;
        event.preventDefault ||= function () { this.defaultPrevented = true; };
        event.stopPropagation ||= function () { this.stopped = true; };
        for (const callback of this.events.get(event.type) || []) callback(event);
        if (event.bubbles && !event.stopped) this.parentNode?.dispatchEvent?.(event);
        return !event.defaultPrevented;
      },
      click() { const event = {type: 'click', bubbles: true}; this.dispatchEvent(event); this.onclick?.(event); return event; },
      focus() { document.activeElement = this; this.dispatchEvent({type: 'focus'}); this.dispatchEvent({type: 'focusin', bubbles: true}); },
      getBoundingClientRect() { return this.rect || {left: 20, top: 20, right: 20 + this.offsetWidth,
        bottom: 20 + this.offsetHeight, width: this.offsetWidth, height: this.offsetHeight}; },
      createTHead() { const item = node('thead'); this.append(item); return item; },
      createTBody() { const item = node('tbody'); this.append(item); return item; },
      insertRow() { const item = node('tr'); this.append(item); return item; },
      insertCell() { const item = node('td'); this.append(item); return item; },
      showModal() { this.open = true; }, close() { this.open = false; this.dispatchEvent({type: 'close'}); },
    };
    result.classList = {add(...names) { result.className = [...new Set(result.className.split(/\s+/).concat(names))].join(' ').trim(); },
      remove(...names) { result.className = result.className.split(/\s+/).filter(name => !names.includes(name)).join(' '); },
      toggle(name, force) { const enabled = force ?? !this.contains(name); this[enabled ? 'add' : 'remove'](name); return enabled; },
      contains(name) { return result.className.split(/\s+/).includes(name); }};
    nodes.push(result); return result;
  }
  const root = node('document'), stack = [root];
  const voidTags = new Set(['area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'source', 'wbr']);
  for (const token of html.match(/<!--[\s\S]*?-->|<![^>]*>|<[^>]+>|[^<]+/g) || []) {
    if (token.startsWith('<!')) continue;
    if (token.startsWith('</')) { if (stack.length > 1) stack.pop(); continue; }
    if (!token.startsWith('<')) { stack.at(-1).append(textNode(decoded(token))); continue; }
    const match = token.match(/^<([\w-]+)([\s\S]*?)\/?\s*>$/); if (!match) continue;
    const item = node(match[1]);
    for (const [, key, doubleValue, singleValue, bareValue] of match[2].matchAll(/([^\s=]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+)))?/g))
      item.setAttribute(key, decoded(doubleValue ?? singleValue ?? bareValue ?? ''));
    stack.at(-1).append(item); if (!voidTags.has(match[1].toLowerCase())) stack.push(item);
  }
  document = root;
  document.documentElement = nodes.find(item => item.tagName === 'HTML');
  document.body = nodes.find(item => item.tagName === 'BODY');
  document.getElementById = id => nodes.find(item => item.id === id) || null;
  document.createElement = node; document.createElementNS = (_, tag) => node(tag);
  document.createTextNode = textNode; document.readyState = 'complete';
  Object.defineProperty(document, 'title', {get: () => nodes.find(item => item.tagName === 'TITLE').textContent,
    set: value => { nodes.find(item => item.tagName === 'TITLE').textContent = value; }});
  const window = {innerWidth: 1280, innerHeight: 800, scrollX: 0, scrollY: 0, document,
    addEventListener(type, callback) { if (!listeners.has(type)) listeners.set(type, []); listeners.get(type).push(callback); },
    dispatchEvent(event) { for (const callback of listeners.get(event.type) || []) callback(event); }, scrollTo() {}};
  const localStorage = {getItem(key) { if (blockedStorage) throw Error('Storage disabled'); return storage.get(key) ?? null; },
    setItem(key, value) { if (blockedStorage) throw Error('Storage disabled'); storage.set(key, value); }};
  Object.assign(window, {localStorage, navigator: {language: browserLanguage, languages: [browserLanguage]}});
  const context = vm.createContext({window, document, localStorage, navigator: window.navigator, console,
    Node: {TEXT_NODE: 3, ELEMENT_NODE: 1}, HTMLElement: Object,
    CustomEvent: class { constructor(type, options = {}) { this.type = type; Object.assign(this, options); } },
    requestAnimationFrame(callback) { callback(); return 1; }, cancelAnimationFrame() {}, setTimeout, clearTimeout,
    fetch: async url => { assert.equal(url, '/i18n.json'); return {ok: true, json: async () => dictionaries}; },
    getComputedStyle: () => ({direction: 'ltr'}),
  });
  window.requestAnimationFrame = context.requestAnimationFrame;
  return {context, document, window, storage, node, nodes};
}

async function runtime(options) {
  const env = environment(options);
  env.context.mapCells = [];
  const start = source.indexOf('// BEGIN CHEMGRIDMAP I18N RUNTIME');
  const end = source.indexOf('// END CHEMGRIDMAP I18N RUNTIME', start);
  assert.ok(start >= 0 && end > start, 'Production i18n runtime markers');
  vm.runInContext(source.slice(start, end), env.context);
  env.api = env.window.ChemGridMapI18n || env.context.ChemGridMapI18n;
  assert.ok(env.api, 'Production i18n API');
  env.context.ChemGridMapI18n = env.api;
  await env.api.ready;
  // 新增控件的依赖函数与旧契约共享同一生产实现。
  const coreAndSourceFunctions = ['previousRunLabel','renderPreviousRuns','loadPreviousRuns','restoredSessionMessage','applySessionSnapshot','openSelectedRun','pollCoreProgress','mismatchMarkerSize','mismatchMarkerTitle','updateMismatchMarker','renderMismatchOverlay','setNewMoleculeTab','updateMeasurementForm','readMeasurement','resetMeasurementForm','renderSinglePreview','previewSingleMolecule','scheduleSinglePreview','addedFeedbackCounts','resultCounts','resultSource','confirmNewMoleculeReplacement','patchNewMoleculeCell','refreshNewMoleculeMutation','addSingleMolecule','renderMeasuredResult','saveMeasuredResult','bindMeasuredMolecules', 'clearSourceHighlight','renderSourceHighlight','highlightSources','renderNewBatchSources','coreRuleName','renderCoreRuleOptions','renderGroupOptions','renderCoreRuleStatus','addCoreDownloads','switchCoreRule','renderCoreComparison','loadCoreComparison'];
  vm.runInContext(coreAndSourceFunctions.map(realFunction).join('\n'), env.context);
  return env;
}

function realFunction(name) {
  const start = source.search(new RegExp(`^(?:async )?function ${name}\\(`, 'm'));
  assert.ok(start >= 0, `Missing function ${name}`);
  const rest = source.slice(start), firstLine = rest.split('\n')[0];
  if (firstLine.trim().endsWith('}')) return firstLine;
  const end = rest.search(/^\}/m);
  assert.ok(end > 0, `Missing closing brace for ${name}`);
  return rest.slice(0, end + 1);
}

const scenarios = {
  async dictionary() {
    assert.deepEqual(Object.keys(dictionaries.en).sort(), Object.keys(dictionaries.zh).sort());
    for (const language of ['en', 'zh']) for (const [key, value] of Object.entries(dictionaries[language])) {
      assert.equal(typeof value, 'string', `${language}:${key}`);
      assert.ok(value.trim(), `${language}:${key} must not be empty`);
    }
    for (const [key, english] of Object.entries(dictionaries.en)) {
      const placeholders = text => [...text.matchAll(/\{(\w+)\}/g)].map(match => match[1]).sort();
      assert.deepEqual(placeholders(english), placeholders(dictionaries.zh[key]), `Placeholder parity: ${key}`);
      if (key.startsWith('d.')) assert.equal(english, key.slice(2), `English text preserved: ${key}`);
    }
  },
  async static() {
    const {nodes} = environment(); let count = 0;
    for (const item of nodes) for (const [attribute, value] of Object.entries(item.attributes)) {
      if (!attribute.startsWith('data-i18n')) continue;
      let keys = attribute === 'data-i18n-text' ? value.split('|') : [value];
      if (attribute === 'data-i18n-option' || attribute === 'data-i18n-child') keys = [value.slice(value.indexOf(':') + 1)];
      for (const key of keys.filter(key => key !== '-')) { assert.ok(Object.hasOwn(dictionaries.en, key), `${attribute}: ${key}`); count++; }
    }
    assert.ok(count > 100, 'Static visible text and attributes carry dictionary keys');
  },
  async help() {
    const {api, document, window} = await runtime();
    const buttons = document.querySelectorAll('[data-help]');
    assert.deepEqual(buttons.map(button => button.dataset.helpControl).sort(), [...controls].sort());
    const tooltip = document.getElementById('parameter-tooltip');
    assert.equal(tooltip.getAttribute('role'), 'tooltip');
    assert.equal(tooltip.parentNode, document.body, 'Tooltip escapes scrolling containers');
    for (const language of ['en', 'zh']) {
      api.setLanguage(language);
      for (const button of buttons) {
        assert.equal(button.tagName, 'BUTTON'); assert.equal(button.getAttribute('type'), 'button');
        assert.equal(button.getAttribute('aria-describedby'), tooltip.id);
        assert.match(button.getAttribute('aria-label'), language === 'en' ? /^Help: .+/ : /^说明：.+/);
        assert.ok(dictionaries[language][button.dataset.help]?.trim(), `Help translation: ${button.dataset.help}`);
        button.dispatchEvent({type: 'pointerenter', pointerType: 'mouse'});
        assert.equal(tooltip.hidden, false);
        assert.equal(tooltip.textContent, api.t(button.dataset.help, {activity: 'pChEMBL'}));
        document.dispatchEvent({type: 'keydown', key: 'Escape'});
        assert.equal(tooltip.hidden, true);
      }
    }
    const checkbox = document.getElementById('activity-confidence');
    const button = buttons.find(item => item.dataset.helpControl === checkbox.id);
    checkbox.checked = true;
    const first = button.click();
    assert.equal(first.defaultPrevented, true); assert.equal(first.stopped, true);
    assert.equal(checkbox.checked, true); assert.equal(tooltip.hidden, false);
    button.click(); assert.equal(tooltip.hidden, true); assert.equal(checkbox.checked, true);
    button.focus(); assert.equal(tooltip.hidden, false);
    document.dispatchEvent({type: 'pointerdown', target: document.body}); assert.equal(tooltip.hidden, true);
    button.click(); buttons[0].click();
    assert.equal(document.querySelectorAll('[role=tooltip]').length, 1);
    assert.equal(tooltip.textContent, api.t(buttons[0].dataset.help));
    for (const width of [360, 390, 1280, 1600]) {
      window.innerWidth = width;
      button.rect = {left: width - 20, right: width - 4, top: 780, bottom: 796, width: 16, height: 16};
      api.showHelp(button);
      const left = parseFloat(tooltip.style.left), top = parseFloat(tooltip.style.top);
      assert.ok(left >= 0 && left + tooltip.offsetWidth <= width, `Tooltip fits ${width}px viewport`);
      assert.ok(top >= 0 && top + tooltip.offsetHeight <= window.innerHeight);
      assert.ok(parseFloat(tooltip.style.maxWidth) <= 320);
    }
    document.getElementById('input-format').value = 'generic'; api.updateHelpLabels();
    const lower = buttons.find(item => item.dataset.helpControl === 'lower');
    api.showHelp(lower); assert.match(lower.getAttribute('aria-label'), /pActivity/); assert.match(tooltip.textContent, /pActivity/);
    assert.doesNotMatch(tooltip.textContent, /pChEMBL/);
  },
  async language() {
    const env = await runtime({saved: 'en', browserLanguage: 'zh-CN'});
    const {api, document, storage} = env;
    assert.equal(api.getLanguage(), 'en', 'Saved language takes priority');
    const build = document.getElementById('build'), title = document.title;
    const previous = build.textContent;
    document.getElementById('target').value = 'CHEMBL205';
    document.getElementById('occupancy').value = '67';
    document.getElementById('audit-dialog').open = true;
    document.querySelector('[data-language=zh]').click();
    assert.equal(api.getLanguage(), 'zh'); assert.equal(document.documentElement.lang, 'zh-CN');
    assert.equal(storage.get('chemgridmap.language'), 'zh');
    assert.notEqual(build.textContent, previous); assert.match(build.textContent, /生成地图/);
    assert.notEqual(document.title, title); assert.match(document.title, /分子工作台/);
    assert.equal(document.querySelector('[data-language=zh]').getAttribute('aria-pressed'), 'true');
    assert.equal(document.querySelector('[data-language=en]').getAttribute('aria-pressed'), 'false');
    assert.equal(document.getElementById('target').value, 'CHEMBL205');
    assert.equal(document.getElementById('occupancy').value, '67');
    assert.equal(document.getElementById('audit-dialog').open, true);
    api.setLanguage('en'); assert.equal(build.textContent, previous); assert.equal(document.title, title);
    assert.equal((await runtime({browserLanguage: 'zh-TW'})).api.getLanguage(), 'zh');
    assert.equal((await runtime({browserLanguage: 'fr-FR'})).api.getLanguage(), 'en');
    const blocked = await runtime({browserLanguage: 'zh-CN', blockedStorage: true});
    assert.equal(blocked.api.getLanguage(), 'zh'); assert.doesNotThrow(() => blocked.api.setLanguage('en'));
    assert.equal(blocked.document.documentElement.lang, 'en');
  },
  async dynamic() {
    const {api, context, document, node} = await runtime();
    const $ = id => { const item = document.getElementById(id); assert.ok(item, `Missing fixture element: ${id}`); return item; };
    const measured = {molecule_id: 'LAB-001', grid_row: 3, grid_col: 4, structure_group: 'G001',
      activity_class: 'active', class_status: 'supported', activity_pchembl: 7.4, pchembl_min: 7.2,
      pchembl_max: 7.6, n_records: 2, neighbor_recall: .75};
    const group = {structure_group: 'G001', size: 1, core_coverage_count: 1, core_smiles: 'c1ccccc1'};
    const report = {input_rows: 3, retained_activity_records: 2, repeated_conflicted_molecules: 0,
      generic_activity: {input_rows: 3, excluded_activity_records: 1, document_column: 'run_id', excluded_reasons: {non_exact_relation: 1}},
      activity_confidence: {document_pair_count: 1, offset_min_shared: 3, discordant_document_pair_count: 1,
        offset_threshold: 1, offset_sign_agreement: .8, offset_display: 'none', document_offsets_top10: [
          {document_a: 'R1', document_b: 'R3', median_offset: 1.18, shared_entries: 5, sign_agreement: 1, systematic_offset: true, touched_entry_count: 24}]},
      warnings: []};
    const record = {document: 'document', activity_class: 'active', Class: 'Class', source_row: 1};
    const recordColumns = Object.keys(record);
    const detail = {molecule: measured, summary: {input_format: 'generic_activity', n_records: 2,
      median_verified: true, record_count_verified: true,
      interpretation: 'Generic activity records: the listed rows are the retained records for this compound; values were converted to pActivity.',
      warnings: [], neighborhood: {representation: 'morgan', representation_metric: 'Tanimoto', effective_k: 5,
        shared_count: 3, representation_count: 5, grid_count: 5, recall: .6}},
      neighbors: [], records: [record], record_columns: recordColumns, total_records: 1, svg: '', neighborhood_svg: ''};
    const transform = {s: 2.5, x: 17, y: -24};
    Object.assign(context, {newMolecules:null,$, t: api.t, translateServer: api.translateServer, job: {report, molecules: 1,
      options: {input_format: 'molecule', compound_column: 'compound_id', document_column: 'run_id',
        representation: 'morgan', projection: 'umap', grid_occupancy: .4, target_id: 'CHEMBL205',
        lower_threshold: 6, upper_threshold: 7, support_level: .8, offset_display: 'none'},
      files: ['map.svg', 'map_class_support.svg', 'map_metrics.csv'], metrics: {neighbor_recall: .75}, output_dir: '/saved/results'},
      inputInfo: null, queryInputInfo: null, sourceInputInfo: null, umapAvailable: null, busy: false, points: [measured],
      structureGroups: [group], highlightedGroup: '', transform, view: 'grid', tileWidth: 100, hasStructures: true,
      currentCell: {row: 3, col: 4}, neighborRelations: new Map(), requestNumber: 1, groupRequestNumber: 0,
      colors: {inactive: '#c44e52', medium: '#e6a23c', active: '#2a9d8f'},
      fadedColors: {inactive: '#ebc0c1', medium: '#f6dfb8', active: '#b8e0db'},
      getJSON() { throw Error('Language changes must use cached results'); },
      api() { throw Error('Language changes must not request exported figures'); },
      restoreOptions() { throw Error('Language changes must preserve current inputs'); },
      selectGroup() {}, compareGroup() {}, highlight() {}, updateNavigation() {}, download() {},
      fit() { throw Error('Language changes must preserve map zoom'); },
      setView() { throw Error('Language changes must preserve the rendered map'); },
      fail(error) { throw error; },
    });
    const names = ['stat', 'format', 'ordinal', 'percentage', 'countLabel', 'legendItem', 'renderColorLegend',
      'renderLegend', 'coreRule', 'selectGroup', 'table', 'selectColumns', 'referenceMapDetails', 'renderReferenceMap', 'newMoleculeWarning', 'newMoleculeInsertedMessage', 'showDataAudit', 'showResult', 'measuredLabel', 'showQuerySelection', 'showSources', 'inspect',
      'updateInputType', 'updateEvidenceOptions', 'applyTransform', 'translateSVGText', 'refreshValidation', 'status', 'refreshLanguage'];
    vm.runInContext(names.map(realFunction).join('\n'), context);
    context.showResult.grouping = {groups: [group], protocol: {method: 'Butina'}};
    context.inspect.lastDetail = detail;
    $('input-format').value = 'generic'; $('occupancy').value = '67'; $('target').value = 'USER-SETTING';
    context.selectColumns('document-column', ['document', 'active', 'Class'], 'document', true);
    context.table([record], $('preview-table'), 12, recordColumns);
    $('inspector').hidden = false; $('color-mode').value = 'activity';
    for (const id of ['audit-dialog', 'export-dialog', 'records-dialog']) $(id).open = true;
    const originalMap = node('svg'); originalMap.dataset.originalMap = 'true'; $('map-stage').append(originalMap);
    await context.refreshLanguage();
    assert.match($('summary').textContent, /Input records/);
    assert.match($('legend').textContent, /Inactive ≤6/);
    assert.match($('evidence-summary').textContent, /Source records/);
    assert.match($('audit-body').textContent, /Evidence units and offsets/);
    assert.doesNotMatch($('projection').querySelector('[value=umap]').textContent, /not installed/,
      'A pending capability response must not claim UMAP is unavailable');
    let languageWork;
    document.addEventListener('chemgridmap:languagechange', () => { languageWork = context.refreshLanguage(); });
    api.setLanguage('zh'); await languageWork;
    assert.match($('summary').textContent, /输入记录/);
    assert.match($('legend').textContent, /非活性/);
    assert.match($('evidence-summary').textContent, /来源记录/);
    assert.match($('audit-body').textContent, /证据单元与偏移/);
    assert.match($('figure-downloads').textContent, /实测活性/);
    assert.doesNotMatch($('summary').textContent, /Input records/);
    assert.equal($('occupancy').value, '67'); assert.equal($('target').value, 'USER-SETTING');
    assert.equal($('document-column').value, 'document');
    assert.deepEqual($('document-column').children.slice(1).map(item => item.textContent), ['document', 'active', 'Class']);
    for (const id of ['record-table', 'preview-table']) {
      assert.deepEqual($(id).querySelectorAll('th').map(item => item.textContent), recordColumns);
      assert.deepEqual($(id).querySelectorAll('td').map(item => item.textContent), Object.values(record).map(String));
    }
    assert.match($('downloads').textContent, /map_metrics\.csv/);
    assert.equal($('map-stage').children[0], originalMap);
    assert.deepEqual(JSON.parse(JSON.stringify(context.transform)), transform);
    assert.deepEqual(JSON.parse(JSON.stringify(context.currentCell)), {row: 3, col: 4});
    for (const id of ['audit-dialog', 'export-dialog', 'records-dialog']) assert.equal($(id).open, true);
    const unmeasured = {...measured, molecule_id: 'QUERY-1', query_id: 'QUERY-1', measurement_status: 'not_measured',
      nearest_measured_id: 'LAB-001', nearest_measured_similarity: .8};
    context.inspect.lastDetail = {...detail, molecule: unmeasured, query_neighbors: []};
    api.setLanguage('en'); await languageWork;
    assert.match($('query-summary').textContent, /Most similar reference ID/);
    api.setLanguage('zh'); await languageWork;
    assert.match($('query-summary').textContent, /最相似的参照分子 ID/);
    assert.equal($('selected-id').textContent, 'QUERY-1');
    assert.equal($('record-section').hidden, true);
    context.umapAvailable = false; await context.refreshLanguage();
    assert.match($('projection').querySelector('[value=umap]').textContent, /未安装/);
    const secondGroup = {...group, structure_group: 'G002'};
    context.showResult.grouping.groups.push(secondGroup);
    context.selectGroup.lastResult = {summary: {structure_group: 'G001'}, core_svg: 'CORE-A'};
    let resolveGroup, requests = 0;
    context.getJSON = () => { requests++; return new Promise(resolve => { resolveGroup = resolve; }); };
    const pending = context.selectGroup('G002');
    api.setLanguage('en'); await languageWork;
    assert.equal(context.highlightedGroup, 'G002');
    assert.match($('group-core-title').textContent, /^G002/);
    assert.doesNotMatch($('group-core-svg').textContent, /CORE-A/);
    assert.equal(requests, 1, 'Changing language does not issue another group request');
    resolveGroup({summary: {structure_group: 'G002'}, core_svg: 'CORE-B'}); await pending;
    assert.equal($('group-core-svg').textContent, 'CORE-B', 'Pending correct group remains eligible to render');
  },
  async server() {
    const {api} = await runtime();
    const examples = [
      ['Preparing input records...', []],
      ['Validating structures and auditing repeated measurements...', []],
      ['Estimating how well source documents support each activity class...', []],
      ["Calculating the projection, assigning grid cells and rendering. UMAP's first run may take longer.", ['UMAP']],
      ['Map ready. Select a cell to inspect its evidence.', []],
      ['Class support uses an assumed σ of 0.64 log units; it was not estimated from these data.', ['0.64', 'σ']],
      ['σ_doc was not estimated: fewer than 5 entries have two or more evidence units; the literature value 0.64 is used', ['σ_doc', '5', '0.64']],
      ['Entry CHEMBL20 supplies 40.6% of σ_doc degrees of freedom and 63.1% of the pooled sum of squares.', ['CHEMBL20', '40.6%', '63.1%']],
      ['Activity class support was not estimated: the input has no document identifiers. Use --assumed-sigma to calculate class support with an assumed uncertainty.', []],
      ['Offset document pair: R1 vs R3, median difference 1.18 log units, n shared 5; the pair touches 24 entries. This entry has records from R1, R3; the display is unchanged because the offset rule is off.', ['R1', 'R3', '1.18', '5', '24']],
    ];
    for (const [text] of examples) assert.equal(api.translateServer(text), text, 'English server information stays exact');
    api.setLanguage('zh');
    for (const [text, preserved] of examples) {
      const translated = api.translateServer(text);
      assert.notEqual(translated, text, `Known server template: ${text}`);
      assert.match(translated, /[\u3400-\u9fff]/);
      for (const value of preserved) assert.ok(translated.includes(value), `Server placeholder preserved: ${value}`);
    }
    const unknown = 'An unknown server diagnostic 42: USER-COLUMN';
    assert.equal(api.translateServer(unknown), unknown);
  },
};

module.exports = {environment, runtime, realFunction, dictionaries};
if (require.main === module) (async () => {
  assert.ok(scenario in scenarios, `Unknown scenario: ${scenario}`);
  await scenarios[scenario]();
  console.log(`U4 ${scenario}: passed`);
})().catch(error => { console.error(error); process.exitCode = 1; });
