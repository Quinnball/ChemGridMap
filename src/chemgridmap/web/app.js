// BEGIN CHEMGRIDMAP I18N RUNTIME
// 仅匹配软件已知的英文消息模板；未命中的服务器信息保持原文。
function installServerTranslations(i18n) {
  let templates;
  const extraTemplates = [
    'd.{term} vs grid neighbors',
    'd.{representation} ({metric}), k={k}: {shared} shared neighbors out of {representationCount} {term} neighbors and {gridCount} grid neighbors. Recall {recall}%.',
    'd.{ordinal} fingerprint neighbor', 'd.{ordinal} representation neighbor', 'd.distance {distance}'
  ];
  const translatedParameters = new Set(['reason', 'touches', 'display', 'term', 'representation', 'metric', 'activityClass']);
  const escape = text => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  i18n.translateServer = function translateServer(message) {
    if (message == null || i18n.getLanguage() === 'en') return message;
    const original = String(message), dictionary = i18n.dictionaries.en;
    if (!templates) {
      templates = Object.keys(dictionary).filter(key => key.startsWith('server.') || extraTemplates.includes(key)).map(key => {
        const names = [], pattern = dictionary[key].split(/(\{\w+\})/g).map(part => {
          if (!/^\{\w+\}$/.test(part)) return escape(part);
          names.push(part.slice(1, -1)); return '(.+?)';
        }).join('');
        return {key, names, pattern: new RegExp('^' + pattern + '$')};
      });
    }
    // 保留状态信息的秒数，并且仅在主消息已识别时翻译时间后缀。
    const elapsed = original.match(/^(.*) \((\d+) s\)$/s);
    if (elapsed) {
      const body = translateServer(elapsed[1]);
      if (body !== elapsed[1]) return body + i18n.t('server.elapsed', {seconds: elapsed[2]});
    }
    // SVG 仅在屏幕上替换说明，保留分子 ID 与源文件字节。
    const projection = original.match(/^(.*) · (active|medium|inactive|unlabelled|not_measured)$/);
    if (projection) return i18n.t('svg.projection_title', {id: projection[1], activityClass: translateServer(projection[2])});
    const card = original.match(/^(pChEMBL [-\d.]+|activity not supplied) \| (.*?)(?: \| (core absent|unaligned: stereo safety))?$/);
    if (card) return i18n.t('svg.card', {activity: translateServer(card[1]), activityClass: translateServer(card[2]), alignment: card[3] ? ' | ' + translateServer(card[3]) : ''});
    const gridTitle = original.match(/^(.*) \| (molecule_id=.*)$/);
    if (gridTitle) {
      const details = gridTitle[2].split('; '), fields = new Set(['molecule_id', 'activity_pchembl', 'n_records', 'pchembl_min', 'pchembl_max', 'is_conflicted', 'has_class_boundary_crossing', 'source_rows']);
      if (details.every(part => fields.has(part.split('=')[0]))) return i18n.t('svg.grid_title', {id: gridTitle[1], details: details.map(part => {
        const at = part.indexOf('='), field = part.slice(0, at);
        return i18n.t('d.' + field, {activityTerm: document.getElementById('input-format')?.value === 'chembl' ? 'pChEMBL' : 'pActivity'}) + part.slice(at);
      }).join('; ')});
    }
    for (const template of templates) {
      const match = original.match(template.pattern);
      if (!match) continue;
      const values = Object.fromEntries(template.names.map((name, index) => {
        let value = match[index + 1];
        if (name === 'ordinal') value = value.replace(/(st|nd|rd|th)$/, '');
        else if (name === 'documents') value = value.split(', ').join(i18n.t('u.idSeparator'));
        else if (translatedParameters.has(name)) {
          const nested = translateServer(value);
          value = nested !== value ? nested : i18n.t('d.' + value);
        }
        return [name, value];
      }));
      return i18n.t(template.key, values);
    }
    // 已知类别和方法名可单独出现；ID、文件名和来源列名不经过此函数。
    if (Object.hasOwn(dictionary, 'd.' + original)) return i18n.t('d.' + original);
    return message;
  };
}

(() => {
  'use strict';
  const dictionaries = {en: {}, zh: {}};
  let language = 'en';
  let activeHelp = null;
  let helpPinned = false;
  let tooltip = null;
  const storageKey = 'chemgridmap.language';
  function initialLanguage() {
    let saved;
    try { saved = window.localStorage.getItem(storageKey); } catch (_) { /* Private storage may be disabled. */ }
    if (saved === 'en' || saved === 'zh') return saved;
    return String(window.navigator.language || '').toLowerCase().startsWith('zh') ? 'zh' : 'en';
  }
  function t(key, values = {}) {
    const template = dictionaries[language][key] ?? dictionaries.en[key] ?? (key.startsWith('d.') ? key.slice(2) : key);
    return String(template).replace(/\{([\w]+)\}/g, (match, name) => Object.prototype.hasOwnProperty.call(values, name) ? String(values[name]) : match);
  }
  function getLanguage() { return language; }
  function activityName() {
    return document.getElementById('input-format')?.value === 'chembl' ? 'pChEMBL' : 'pActivity';
  }
  function helpText(button) { return t(button.dataset.help, {activity: activityName()}); }
  function updateHelpLabels() {
    document.querySelectorAll('[data-help]').forEach(button => {
      const id = button.dataset.helpControl;
      const name = id === 'lower' || id === 'upper' ? t(`u.${id}`, {activity: activityName()}) : t(`c.${id}`);
      button.setAttribute('aria-label', t('u.help', {name}));
    });
    if (activeHelp) {
      tooltip.textContent = helpText(activeHelp);
      positionHelp();
    }
  }
  function applyTranslations(root = document) {
    root.querySelectorAll('[data-i18n]').forEach(element => { element.textContent = t(element.dataset.i18n); });
    // Translate direct text nodes only: icons, controls and live values remain intact.
    root.querySelectorAll('[data-i18n-text]').forEach(element => {
      const keys = element.dataset.i18nText.split('|');
      const nodes = Array.from(element.childNodes).filter(node => node.nodeType === 3 && node.textContent.trim());
      nodes.forEach((node, index) => {
        if (keys[index] && keys[index] !== '-') {
          const leading = node.textContent.match(/^\s*/)[0], trailing = node.textContent.match(/\s*$/)[0];
          node.textContent = leading + t(keys[index]) + trailing;
        }
      });
    });
    for (const attribute of ['aria-label', 'placeholder', 'alt', 'title']) {
      root.querySelectorAll(`[data-i18n-${attribute}]`).forEach(element => {
        element.setAttribute(attribute, t(element.getAttribute(`data-i18n-${attribute}`)));
      });
    }
    root.querySelectorAll('[data-i18n-child]').forEach(element => {
      const [id, key] = element.dataset.i18nChild.split(':');
      const child = document.getElementById(id);
      if (child) child.textContent = t(key);
    });
    root.querySelectorAll('[data-i18n-option]').forEach(element => {
      const [id, key] = element.dataset.i18nOption.split(':');
      const select = document.getElementById(id);
      if (select?.options.length) select.options[0].textContent = t(key);
    });
    updateHelpLabels();
  }
  function setLanguage(next) {
    language = next === 'zh' || next === 'zh-CN' ? 'zh' : 'en';
    try { window.localStorage.setItem(storageKey, language); } catch (_) { /* Language switching works without storage. */ }
    document.documentElement.lang = language === 'zh' ? 'zh-CN' : 'en';
    applyTranslations();
    document.querySelectorAll('[data-language]').forEach(button => {
      button.setAttribute('aria-pressed', String(button.dataset.language === language));
    });
    document.dispatchEvent(new CustomEvent('chemgridmap:languagechange', {detail: {language}}));
    return language;
  }
  function positionHelp() {
    if (!activeHelp || !tooltip || tooltip.hidden) return;
    const viewport = window.visualViewport;
    const left = viewport?.offsetLeft || 0, top = viewport?.offsetTop || 0;
    const width = viewport?.width || window.innerWidth, height = viewport?.height || window.innerHeight;
    const margin = 8, gap = 7;
    tooltip.style.maxWidth = `${Math.max(0, Math.min(320, width - margin * 2))}px`;
    tooltip.style.maxHeight = `${Math.max(0, height - margin * 2)}px`;
    const anchor = activeHelp.getBoundingClientRect();
    const box = tooltip.getBoundingClientRect();
    const x = Math.min(Math.max(anchor.left, left + margin), left + width - box.width - margin);
    const below = anchor.bottom + gap;
    const preferred = below + box.height <= top + height - margin ? below : anchor.top - box.height - gap;
    const y = Math.min(Math.max(preferred, top + margin), top + height - box.height - margin);
    tooltip.style.left = `${x}px`;
    tooltip.style.top = `${y}px`;
  }
  function closeHelp() {
    activeHelp = null; helpPinned = false;
    if (tooltip) tooltip.hidden = true;
  }
  function showHelp(button) {
    if (!tooltip) return;
    // A modal dialog is in the top layer: its help must share that layer.
    const dialog = button.closest('dialog');
    const host = dialog?.open ? dialog : document.body;
    if (tooltip.parentNode !== host) host.append(tooltip);
    activeHelp = button; helpPinned = false;
    tooltip.textContent = helpText(button);
    tooltip.hidden = false;
    positionHelp();
  }
  function initializeHelp() {
    tooltip = document.createElement('div');
    tooltip.id = 'parameter-tooltip'; tooltip.className = 'parameter-tooltip';
    tooltip.setAttribute('role', 'tooltip'); tooltip.hidden = true;
    document.body.append(tooltip);
    document.querySelectorAll('[data-help]').forEach(button => {
      button.addEventListener('pointerenter', event => {
        if (event.pointerType !== 'touch') showHelp(button);
      });
      button.addEventListener('pointerleave', event => {
        if (activeHelp === button && !helpPinned && document.activeElement !== button && !tooltip.contains(event.relatedTarget)) closeHelp();
      });
      button.addEventListener('focus', () => { showHelp(button); });
      button.addEventListener('blur', event => { if (!tooltip.contains(event.relatedTarget)) closeHelp(); });
      button.addEventListener('click', event => {
        // Prevent the enclosing label from activating its checkbox/select.
        event.preventDefault(); event.stopPropagation();
        if (activeHelp === button && helpPinned) closeHelp();
        else { showHelp(button); helpPinned = true; }
      });
    });
    tooltip.addEventListener('pointerleave', event => {
      if (!helpPinned && event.relatedTarget !== activeHelp && document.activeElement !== activeHelp) closeHelp();
    });
    document.addEventListener('pointerdown', event => {
      if (activeHelp && !activeHelp.contains(event.target) && !tooltip.contains(event.target)) closeHelp();
    });
    document.addEventListener('keydown', event => {
      if (event.key !== 'Escape') return;
      // 第一次 Esc 只收起对话框内的提示；下一次仍可关闭对话框。
      if (activeHelp?.closest('dialog')?.open) event.preventDefault();
      closeHelp();
    });
    document.querySelectorAll('dialog').forEach(dialog => {
      dialog.addEventListener('close', () => {
        if (activeHelp?.closest('dialog') === dialog) closeHelp();
      });
    });
    window.addEventListener('resize', positionHelp);
    window.addEventListener('scroll', positionHelp, true);
    window.visualViewport?.addEventListener('resize', positionHelp);
    window.visualViewport?.addEventListener('scroll', positionHelp);
    document.getElementById('input-format')?.addEventListener('change', updateHelpLabels);
  }
  function initialize() {
    initializeHelp();
    document.querySelectorAll('[data-language]').forEach(button => {
      button.addEventListener('click', () => setLanguage(button.dataset.language));
    });
    setLanguage(initialLanguage());
  }
  // Server-template translation is supplied alongside the same dictionary.
  function translateServer(value) { return value; }
  const api = {t, setLanguage, getLanguage, applyTranslations, updateHelpLabels, closeHelp, showHelp, positionHelp, dictionaries, translateServer};
  installServerTranslations(api);
  window.ChemGridMapI18n = api;
  api.ready = fetch('/i18n.json').then(response => {
    if (!response.ok) throw new Error('Unable to load interface translations.');
    return response.json();
  }).then(data => {
    Object.assign(dictionaries.en, data.en);
    Object.assign(dictionaries.zh, data.zh);
    if (document.readyState === 'loading') {
      return new Promise(resolve => document.addEventListener('DOMContentLoaded', () => { initialize(); resolve(api); }, {once: true}));
    }
    initialize();
    return api;
  });
})();
// END CHEMGRIDMAP I18N RUNTIME

async function initializeWorkspace() {
await window.ChemGridMapI18n.ready;
'use strict';
const $ = id => document.getElementById(id);
const t = (key, values) => window.ChemGridMapI18n.t(key, values);
const translateServer = text => window.ChemGridMapI18n.translateServer(text);
const token = location.hash.slice(1) || sessionStorage.getItem('chemgridmap-token');
if (token) sessionStorage.setItem('chemgridmap-token', token);
history.replaceState(null, '', location.pathname);
let busy = false, job = null, points = [], gridSVG = '', view = 'grid', currentCell = null, requestNumber = 0;
let transform = {s: 1, x: 0, y: 0}, drawingSize = {w: 1, h: 1};
let tileWidth = 100, hasStructures = false;
let inputInfo = null, queryInputInfo = null, sourceInputInfo = null, page = 'import', umapAvailable = null, demoAvailable = false;
let structureGroups = [], neighborRelations = new Map();
let newMolecules = null, newListSort = {key:'query_id', direction:1};
let highlightedGroup = '', groupRequestNumber = 0, mapCells = [];
const colors = {active: '#2a9d8f', medium: '#e6a23c', inactive: '#c44e52', unlabelled: '#b8c2cc'};
// Uncertain classes keep their hue at reduced saturation; reliability uses one sequential blue scale.
const fadedColors = {active: '#b8e0db', medium: '#f6dfb8', inactive: '#ebc0c1', unlabelled: '#dfe4e9'};
const reliabilityBins = [[0, '#cde1f5'], [.2, '#a0bdda'], [.4, '#7599bf'], [.6, '#4a77a4'], [.8, '#1d5689']];
function reliabilityColor(value) {
  if (value == null || Number.isNaN(value)) return colors.unlabelled;
  return reliabilityBins.filter(([low]) => value >= low - 1e-9).pop()[1];
}

function showPage(name) {
  if (name === 'configure' && !inputInfo || name === 'explore' && !job && !busy) return;
  page = name;
  for (const key of ['import', 'configure', 'explore']) {
    $('page-' + key).hidden = key !== name;
    const button = $('step-' + key);
    if (key === name) button.setAttribute('aria-current', 'step');
    else button.removeAttribute('aria-current');
    button.classList.toggle('completed', key === 'import' && !!inputInfo || key === 'configure' && !!job);
  }
  window.scrollTo({top: 0});
  if (name === 'explore') requestAnimationFrame(fit);
}
function updateNavigation() {
  $('step-import').disabled = busy;
  $('step-configure').disabled = busy || !inputInfo;
  $('step-explore').disabled = !job && !busy;
  $('replace-input').disabled = busy;
  $('adjust-settings').disabled = busy;
  $('open-audit').disabled = $('open-export').disabled = !job || busy;
  updateNewNavigation();
  $('build').disabled = busy || !inputInfo;
  $('demo').disabled = busy || !demoAvailable;
  $('projection').querySelector('[value=umap]').disabled = !umapAvailable;
  document.querySelectorAll('#find-molecule input,#find-molecule button,.tabs button,.zoom-tools button').forEach(el => el.disabled = !job || busy);
  document.querySelectorAll('#group-controls input,#group-controls select,#group-controls button').forEach(el => el.disabled = !job || busy);
  $('compare-group').disabled = !job || busy || !highlightedGroup;
}
function densityPreview() {
  $('occupancy-value').textContent = $('occupancy').value + '%';
  $('density-preview').replaceChildren();
  const count = Math.round(Number($('occupancy').value) / 4);
  for (let i = 0; i < 25; i++) {
    const cell = document.createElement('i');
    cell.className = (i * 7) % 25 < count ? 'filled' : '';
    $('density-preview').append(cell);
  }
}
function closeInspector() {
  $('inspector').hidden = true; $('map-layout').classList.remove('has-selection');
  currentCell = null; neighborRelations.clear(); requestNumber++; highlight(); requestAnimationFrame(fit);
}

async function api(path, options = {}) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const response = await fetch(path, {...options, headers: {...options.headers, 'X-ChemGridMap-Token': token || ''}});
  if (!response.ok) {
    const data = await response.json();
    throw new Error(data.error || 'The request failed.');
  }
  return response;
}
const getJSON = async (path, options) => (await api(path, options)).json();
function status(message, kind = '') {
  status.last = {message, kind};
  message = typeof message === 'function' ? message() : translateServer(message);
  $('status').textContent = message; $('status').className = kind;
  $('status').hidden = !message || kind === 'running';
  if (kind === 'running') $('build-progress').textContent = message;
}
function fail(error) { status(error.message, 'error'); }
function setBusy(value) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  busy = value;
  document.querySelectorAll('#options input,#options select,#options button,#new-molecules-form input,#new-molecules-form select,#new-molecules-form button,#csv,#demo').forEach(el => el.disabled = value);
  $('quit').disabled = value;
  $('build').textContent = value ? t("d.Building map...") : t("d.Generate map →");
  updateNavigation(); updateEvidenceOptions();
}
function table(rows, container, columnLimit = 999, orderedColumns = null) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  container.replaceChildren();
  if (!rows.length) { container.textContent = t("d.No records available."); return; }
  const columns = (orderedColumns || Object.keys(rows[0])).slice(0, columnLimit), el = document.createElement('table');
  const sourceTable = container === $('preview-table') || container === $('record-table');
  const head = el.createTHead().insertRow();
  for (const col of columns) { const th = document.createElement('th'); th.textContent = sourceTable ? col : t('d.' + col, {activityTerm: job?.options.input_format === 'chembl' ? 'pChEMBL' : 'pActivity'}); head.append(th); }
  const body = el.createTBody();
  for (const row of rows) {
    const tr = body.insertRow();
    for (const col of columns) { const td = tr.insertCell(); const value = row[col] ?? '—'; td.textContent = sourceTable ? value : typeof value === 'boolean' ? t('d.' + value) : typeof value !== 'string' ? value : (['Metric','neighbor_relation','activity_class','class_status','measurement_status','Class','Display class','Class status'].includes(col) || col === 'Value' && row.Metric !== 'map_name' ? t('d.' + value) : value); td.title = td.textContent; }
  }
  container.append(el);
}
function selectColumns(id, columns, preferred, optional = false) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $(id).replaceChildren();
  for (const name of (optional ? [''].concat(columns) : columns)) {
    const option = document.createElement('option'); option.value = name; option.textContent = name || t("d.None"); $(id).append(option);
  }
  $(id).value = columns.includes(preferred) ? preferred : (optional ? '' : columns[0]);
}
function resetResult() {
  inspect.lastDetail = null; compareGroup.lastResult = null; selectGroup.lastResult = null;
  job = null; currentCell = null; points = []; gridSVG = ''; requestNumber++;
  structureGroups = []; neighborRelations.clear(); $('group-controls').hidden = true;
  highlightedGroup = ''; groupRequestNumber++; mapCells = [];
  $('group-select').value = ''; $('group-core-card').hidden = true; $('group-core-svg').replaceChildren();
  $('measurement-filter').value = 'all'; $('measurement-filter-control').hidden = true;
  if ($('comparison-dialog').open) $('comparison-dialog').close();
  ['summary','audit','metrics','selection','record-section','download-all'].forEach(id => $(id).hidden = true);
  $('selection-hint').hidden = false; $('downloads').replaceChildren(); $('map-stage').replaceChildren();
  $('empty-map').hidden = false; $('saved-location').textContent = '';
  $('figure-downloads').replaceChildren(); $('legend').replaceChildren();
  $('run-caption').textContent = ''; closeInspector(); updateNavigation();
}
async function loadInput(demo = false, droppedFile = null) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  try {
    if (busy) return;
    const file = droppedFile || $('csv').files[0];
    if (!demo && !file) return;
    setBusy(true); status("Reading the CSV...");
    const input = await getJSON(demo ? '/api/demo' : '/api/input', {
      method: 'POST', body: demo ? '' : file, headers: demo ? {} : {'X-Filename': encodeURIComponent(file.name)}
    });
    resetResult(); newMolecules = null; showQueryInput(null); showSourceInput(null); showInput(input, demo);
    status(input.targets.length > 1 ? 'Multiple targets detected. Select a target before building; they will not be pooled.' : '');
    setBusy(false);
    showPage('configure');
  } catch (error) { setBusy(false); $('build').disabled = !job && $('input-summary').hidden; fail(error); }
  finally { $('csv').value = ''; }
}
function showInput(input, demo = false) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
    inputInfo = input;
    $('input-summary').hidden = false;
    $('input-name').textContent = input.filename;
    $('input-summary').textContent = t("d.CSV loaded · stored only on this computer");
    $('input-stats').replaceChildren(stat(input.rows, t("d.Input rows")),stat(input.columns.length, t("d.Columns")));
    $('target').value = input.targets.length === 1 ? input.targets[0] : '';
    $('target-list').replaceChildren();
    input.targets.forEach(target => { const option = document.createElement('option'); option.value = target; $('target-list').append(option); });
    selectColumns('smiles-col', input.columns, input.detected.smiles || 'smiles');
    selectColumns('value-col', input.columns, input.detected.pchembl || 'activity_pchembl', true);
    selectColumns('compound-column', input.columns, input.detected.molecule_id || '', true);
    selectColumns('document-column', input.columns, input.detected.document_id || '', true);
    selectColumns('reference-column', input.columns, '', true);
    selectColumns('relation-column', input.columns, input.detected.standard_relation || '', true);
    selectColumns('label-col', input.columns, 'activity_class', true);
    selectColumns('unit-column', input.columns, '', true);
    $('value-unit').value = 'p'; $('use-assumed-sigma').checked = false; $('assumed-sigma').value = .64;
    $('input-format').value = demo || (input.detected.pchembl && (input.detected.standard_type || input.detected.target_id)) ? 'chembl' : input.detected.molecule_id ? 'generic' : 'molecule';
    updateInputType();
    table(input.preview, $('preview-table'), 12); $('preview').hidden = false;
    updateNavigation();
}
function updateInputType() {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const chembl = $('input-format').value === 'chembl';
  $('chembl-options').hidden = !chembl; $('molecule-options').hidden = chembl;
  $('lower-label').textContent = t("d.Lower {term}", {term: chembl ? 'pChEMBL' : 'pActivity'});
  $('upper-label').textContent = t("d.Upper {term}", {term: chembl ? 'pChEMBL' : 'pActivity'});
  updateEvidenceOptions();
  window.ChemGridMapI18n.updateHelpLabels?.();
}
function updateEvidenceOptions() {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const moleculeTable = $('input-format').value === 'molecule';
  $('reference-column-control').hidden = $('input-format').value !== 'generic';
  $('reference-column').disabled = busy || $('input-format').value !== 'generic';
  $('source-input-control').hidden = $('input-format').value !== 'generic';
  $('source-csv').disabled = busy || $('input-format').value !== 'generic';
  $('remove-source-input').disabled = busy || $('input-format').value !== 'generic';
  $('compound-column').disabled = busy || moleculeTable;
  $('compound-column').required = $('input-format').value === 'generic';
  $('value-unit').disabled = busy || !!$('unit-column').value;
  if (moleculeTable) $('use-assumed-sigma').checked = false;
  $('use-assumed-sigma').disabled = busy || moleculeTable;
  $('assumed-sigma').disabled = busy || moleculeTable || !$('use-assumed-sigma').checked;
  $('assumed-sigma').setCustomValidity(Number($('assumed-sigma').value) > 0 ? '' : t("d.Assumed σ must be greater than zero."));
  $('assumed-sigma-unavailable').hidden = !moleculeTable;
}
function updateNewNavigation() {
  const available = !!job && job.options.representation === 'morgan' && !job.options.coordinates && !job.options.embedding_columns;
  $('open-new-molecules').disabled = busy || !available;
  $('new-molecules-availability').textContent = !job ? t('u8.buildFirst') : !available ? t('u8.morganOnly') : '';
  $('open-new-list').hidden = !newMolecules?.rows?.length;
  $('open-new-list').disabled = busy;
  $('remove-query-input').hidden = !queryInputInfo && !newMolecules?.rows?.length;
  $('remove-query-input').disabled = busy;
  $('insert-new-molecules').disabled = busy || !queryInputInfo || !available;
  const measured = !!$('new-activity-col').value;
  for (const id of ['new-unit','new-unit-col','new-relation-col','new-document-col']) $(id).disabled = busy || !measured;
  $('new-unit').disabled = busy || !measured || !!$('new-unit-col').value;
}
function referenceMapDetails() {
  return {
    description: job.reference_map?.description || job.options.target_id || inputInfo?.filename || '—',
    count: job.reference_map?.count ?? (job.molecules - (job.report.query_molecules?.mapped || 0))
  };
}
function renderReferenceMap() {
  const details = job ? referenceMapDetails() : null;
  for (const id of ['new-molecules-reference','new-list-reference']) $(id).textContent = details ? t('u8d.referenceMap', details) : '';
}
function newMoleculeWarning(warning) {
  if (typeof warning === 'string') return translateServer(warning);
  if (warning.code === 'structural_mismatch') return t('u8d.structuralMismatch', {...warning, median: format(warning.median)});
  if (warning.code === 'target_mismatch') return t('u8d.targetMismatch', {...warning, targets: warning.targets.join(t('u8d.separator'))});
  return '';
}
function newMoleculeInsertedMessage(report) {
  const counts = ['active','medium','inactive','mixed','insufficient','none']
    .filter(category => report.neighbor_evidence_counts?.[category] > 0)
    .map(category => `${t('u8d.' + category)} ${report.neighbor_evidence_counts[category]}`).join(t('u8d.separator'));
  return [t('u8d.inserted', {placed: report.mapped || 0, counts: counts || '0'}), ...(report.warnings || []).map(newMoleculeWarning)].join(' ');
}
function showQueryInput(input, saved = {}) {
  queryInputInfo = input;
  $('query-input-info').hidden = !input;
  if (!input) { $('query-csv').value = ''; return; }
  $('query-input-name').textContent = input.filename;
  $('query-input-rows').textContent = t('d.{count} rows', {count: input.rows.toLocaleString()});
  selectColumns('query-smiles-col', input.columns, saved.smiles_col || input.detected?.smiles || 'smiles');
  selectColumns('query-id-col', input.columns, saved.id_col || ['query_id','gen_id','compound_id','molecule_id','id'].find(name => input.columns.includes(name)) || '', true);
  for (const [id,key] of [['new-activity-col','activity_col'],['new-unit-col','unit_col'],['new-relation-col','relation_col'],['new-document-col','document_col'],['new-label-col','label_col']]) selectColumns(id, input.columns, saved[key] || saved[({activity_col:'value_col',unit_col:'unit_column',relation_col:'relation_column',document_col:'document_column',label_col:'label_column'})[key]] || '', true);
  $('new-unit').value = saved.unit || saved.value_unit || 'p';
}
async function loadQueryInput() {
  if (busy || !$('query-csv').files[0]) return;
  try {
    const file = $('query-csv').files[0]; setBusy(true);
    const input = await getJSON('/api/new-molecules/input', {method:'POST', body:file, headers:{'X-Filename':encodeURIComponent(file.name)}});
    showQueryInput(input); $('new-molecules-feedback').textContent = t('u8.loaded');
  } catch (error) { fail(error); $('new-molecules-feedback').textContent = translateServer(error.message); }
  finally { setBusy(false); }
}
function newMoleculeOptions() {
  const options = {};
  for (const [id,key] of [['query-smiles-col','smiles_col'],['query-id-col','id_col'],['new-activity-col','activity_col'],['new-unit','unit'],['new-unit-col','unit_col'],['new-relation-col','relation_col'],['new-document-col','document_col'],['new-label-col','label_col']]) options[key] = $(id).value || null;
  if (options.unit_col) options.unit = 'p';
  return options;
}
async function insertNewMolecules(event) {
  event.preventDefault();
  if (busy || !job || !queryInputInfo) return;
  try {
    const options = newMoleculeOptions(); setBusy(true);
    $('new-molecules-feedback').textContent = t('u8.inserting');
    job = await getJSON('/api/new-molecules', {method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(options)});
    // 只由本次插入成功触发；刷新、重新建图与语言切换不会重置用户所选图层。
    $('color-mode').value = 'neighbor_evidence';
    closeInspector(); await showResult(); $('new-molecules-dialog').close();
    const inserted = job.report.query_molecules || {};
    status(() => newMoleculeInsertedMessage(inserted), inserted.warnings?.length ? 'warning' : '');
  } catch (error) { fail(error); $('new-molecules-feedback').textContent = translateServer(error.message); }
  finally { setBusy(false); }
}
async function removeQueryInput() {
  if (busy) return;
  try {
    setBusy(true); const result = await getJSON('/api/new-molecules/clear', {method:'POST'});
    if (result.job) job = result.job; else if (result.state) job = result;
    else job = await getJSON('/api/job');
    newMolecules = null; showQueryInput(null); closeInspector();
    if ($('color-mode').value === 'neighbor_evidence') $('color-mode').value = 'activity';
    if (job) await showResult();
    status(t('u8.removed'));
  } catch (error) { fail(error); }
  finally { setBusy(false); }
}
function bindNewMolecules() {
$('open-new-molecules').addEventListener('click', () => { renderReferenceMap(); $('new-molecules-feedback').textContent = ''; $('new-molecules-dialog').showModal(); });
$('new-molecules-form').addEventListener('submit', insertNewMolecules);
$('new-activity-col').addEventListener('change',updateNewNavigation);
$('new-unit-col').addEventListener('change',updateNewNavigation);
$('open-new-list').addEventListener('click', () => { renderNewList(true); $('new-list-dialog').showModal(); });
for (const id of ['new-class-filter','new-label-filter']) $(id).addEventListener('change', () => renderNewList());
$('download-new-list').addEventListener('click',downloadNewList);
}
function newListColumns() {
  const columns = [['query_id','ID'],['label',t('u8.labelFilter')],['status',t('u8.status')],['neighbor_evidence',t('u8.neighborEvidence')],['nearest_reference_id',t('u8.referenceId')],['tanimoto','Tanimoto'],['structure_group',t('d.Structure group')]];
  if ((newMolecules?.rows || []).some(row => row.measurement_status === 'measured')) columns.push(['activity_pchembl',t('u8.ownMedian')],['activity_class',t('d.Activity class')],['neighbor_comparison',t('u8.comparisonColumn')]);
  return columns;
}
function visibleNewMolecules() {
  const category = $('new-class-filter').value, label = $('new-label-filter').value;
  const rows = (newMolecules?.rows || []).filter(row => (!category || row.neighbor_evidence === category) && (!label || String(row.label ?? '') === label));
  return rows.slice().sort((a,b) => {
    const x = a[newListSort.key] ?? '', y = b[newListSort.key] ?? '';
    const order = typeof x === 'number' && typeof y === 'number' ? x-y : String(x).localeCompare(String(y), undefined, {numeric:true});
    return (order || String(a.query_id).localeCompare(String(b.query_id))) * newListSort.direction;
  });
}
function newListValue(row, key) {
  const value = row[key];
  if (key === 'status' || key === 'neighbor_evidence' || key === 'neighbor_comparison') return value ? t('u8.' + value) : '—';
  if (key === 'activity_class') return value ? t('d.' + value) : '—';
  return format(value);
}
function renderNewList(resetFilters = false) {
  renderReferenceMap();
  if (resetFilters) {
    const previous = [$('new-class-filter').value, $('new-label-filter').value];
    for (const [id,values] of [['new-class-filter',['active','medium','inactive','mixed','insufficient','none']],['new-label-filter',[...new Set((newMolecules?.rows || []).map(row => row.label).filter(value => value != null && value !== ''))].sort()]]) {
      $(id).replaceChildren();
      for (const value of ['',...values]) { const option = document.createElement('option'); option.value = value; option.textContent = !value ? t('u8.all') : id === 'new-class-filter' ? t('u8.' + value) : String(value); $(id).append(option); }
    }
    $('new-class-filter').value = previous[0]; $('new-label-filter').value = previous[1];
  }
  $('new-list-counts').replaceChildren();
  for (const category of ['active','medium','inactive','mixed','insufficient','none']) {
    const count = (newMolecules?.rows || []).filter(row => row.status === 'mapped' && row.neighbor_evidence === category).length;
    const item = document.createElement('span'); item.textContent = t('d.{label}: {value}', {label:t('u8.'+category), value:count}); item.dataset.evidenceCategory = category; $('new-list-counts').append(item);
  }
  const rows = visibleNewMolecules(), columns = newListColumns(), container = $('new-list-table');
  container.replaceChildren(); const element = document.createElement('table'), head = element.createTHead().insertRow(), body = element.createTBody();
  for (const [key,label] of columns) {
    const th = document.createElement('th'), button = document.createElement('button');
    th.setAttribute('aria-sort', newListSort.key === key ? (newListSort.direction === 1 ? 'ascending' : 'descending') : 'none');
    button.type = 'button'; button.textContent = label; button.setAttribute('aria-label', t('u8.sort',{column:label}));
    button.onclick = () => { newListSort = {key, direction:newListSort.key === key ? -newListSort.direction : 1}; renderNewList(); };
    th.append(button); head.append(th);
  }
  for (const row of rows) {
    const tr = body.insertRow(); tr.dataset.queryId = row.query_id; tr.tabIndex = 0; tr.setAttribute('role','button');
    for (const [key] of columns) tr.insertCell().textContent = newListValue(row,key);
    const locate = () => {
      if (row.status === 'invalid_smiles') { $('new-list-note').textContent = t('u8.invalidNote'); return; }
      $('new-list-dialog').close(); if (view !== 'grid') setView('grid'); Promise.resolve(inspect('id='+encodeURIComponent(row.query_id))).then(centerOnCurrentCell);
    };
    tr.onclick = locate; tr.onkeydown = event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); locate(); } };
  }
  container.append(element); $('new-list-note').textContent = t('u8.displayed',{count:rows.length});
}
function centerOnCurrentCell() {
  if (!currentCell) return;
  const cell = mapCells.find(cell => Number(cell.dataset.gridRow) === currentCell.row && Number(cell.dataset.gridCol) === currentCell.col);
  if (!cell) return;
  const cellBox = cell.getBoundingClientRect(), viewportBox = $('map-viewport').getBoundingClientRect();
  transform.x += viewportBox.left + viewportBox.width / 2 - cellBox.left - cellBox.width / 2;
  transform.y += viewportBox.top + viewportBox.height / 2 - cellBox.top - cellBox.height / 2;
  applyTransform();
}
function downloadNewList() {
  const columns = newListColumns(), rows = visibleNewMolecules();
  const quote = value => '"' + String(value ?? '').replace(/"/g,'""') + '"';
  const csv = [columns.map(([key]) => key),...rows.map(row => columns.map(([key]) => row[key] ?? ''))].map(row => row.map(quote).join(',')).join('\r\n')+'\r\n';
  const url = URL.createObjectURL(new Blob(['\uFEFF'+csv],{type:'text/csv;charset=utf-8'})), link = document.createElement('a');
  link.href = url; link.download = 'new_molecules_filtered.csv'; link.click(); URL.revokeObjectURL(url);
}
function showSourceInput(input) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  sourceInputInfo = input;
  $('source-input-info').hidden = !input;
  $('source-csv').value = '';
  $('source-input-name').textContent = input ? input.filename : '';
  $('source-input-rows').textContent = input ? t('d.{count} rows', {count: input.rows.toLocaleString()}) : '';
}
async function loadSourceInput() {
  if (busy || !$('source-csv').files[0]) return;
  try {
    const file = $('source-csv').files[0]; setBusy(true);
    const input = await getJSON('/api/source-input', {method:'POST', body:file, headers:{'X-Filename':encodeURIComponent(file.name)}});
    resetResult(); showSourceInput(input); showPage('configure');
    status('Source lookup file loaded. Generate the map to apply this input.');
  } catch (error) { fail(error); }
  finally { setBusy(false); $('source-csv').value = ''; }
}
async function removeSourceInput() {
  if (busy) return;
  try {
    setBusy(true); await getJSON('/api/source-input/clear', {method:'POST'});
    resetResult(); showSourceInput(null); showPage('configure');
    status('Source lookup file removed. Generate the map to update the result.');
  } catch (error) { fail(error); }
  finally { setBusy(false); }
}
function updateGrouping() {
  const grouping = $('grouping').value;
  document.querySelectorAll('[data-butina-setting]').forEach(el => el.hidden = grouping !== 'butina');
  document.querySelectorAll('[data-external-setting]').forEach(el => el.hidden = grouping !== 'external');
  document.querySelectorAll('[data-core-setting]').forEach(el => el.hidden = grouping === 'none');
  document.querySelectorAll('[data-fmcs-setting]').forEach(el => el.hidden = grouping === 'none' || $('core-rule').value !== 'fmcs');
}
function options() {
  const config = {input_format: $('input-format').value === 'generic' ? 'molecule' : $('input-format').value, target_id: $('target').value.trim() || null,
    activity_type: $('activity-type').value, smiles_col: $('smiles-col').value,
    value_col: $('value-col').value || null, label_col: $('label-col').value || null,
    compound_column: $('input-format').value === 'molecule' ? null : $('compound-column').value || null, document_column: $('document-column').value || null,
    relation_column: $('relation-column').value || null,
    representation: $('representation').value, projection: $('projection').value,
    lower_threshold: Number($('lower').value), upper_threshold: Number($('upper').value),
    conflict_range_threshold: Number($('conflict-range').value), random_state: Number($('seed').value),
    grid_occupancy: Number($('occupancy').value) / 100, render_detail: $('detail').value,
    structure_policy: $('structure-policy').value, molecule_identity: $('identity').value,
    grouping: $('grouping').value, group_similarity: Number($('group-similarity').value),
    cluster_col: $('cluster-col').value.trim() || null, cluster_source: $('cluster-source').value.trim() || null,
    grid_refinement: $('grid-refinement').value, core_rule: $('core-rule').value,
    core_min_coverage: Number($('core-coverage-target').value), fmcs_timeout: Number($('fmcs-timeout').value),
    activity_confidence: $('activity-confidence').checked, evidence_model: $('evidence-model').value,
    merge_identical_documents: $('merge-identical-documents').checked,
    support_level: Number($('support-level').value), offset_display: $('offset-display').value,
    offset_min_shared: Number($('offset-min-shared').value),
    offset_min_difference: Number($('offset-min-difference').value),
    offset_sign_agreement: Number($('offset-sign-agreement').value),
    identity_conflicts: $('identity-conflicts').value,
    output_formats: $('raster').checked ? ['svg','png','pdf'] : ['svg']};
  if (config.input_format !== 'chembl') {
    if ($('unit-column').value) config.unit_column = $('unit-column').value;
    else if ($('value-unit').value !== 'p') config.value_unit = $('value-unit').value;
  }
  if ($('input-format').value === 'generic' && $('reference-column').value) config.reference_column = $('reference-column').value;
  if (!$('use-assumed-sigma').disabled && $('use-assumed-sigma').checked) config.assumed_sigma = Number($('assumed-sigma').value);
  return config;
}
async function build(event) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  event.preventDefault();
  if (busy) return;
  try {
    const config = options();
    await getJSON('/api/build', {method:'POST', body:JSON.stringify(config), headers:{'Content-Type':'application/json'}});
    resetResult(); setBusy(true); poll();
    showPage('explore'); status("Preparing input records...", 'running');
  } catch (error) { fail(error); }
}
async function poll() {
  try {
    const state = await getJSON('/api/job');
    const elapsed = state.started ? ` (${Math.round(Date.now()/1000 - state.started)} s)` : '';
    status(state.message + (state.state === 'running' ? elapsed : ''), state.state === 'error' ? 'error' : state.state === 'running' ? 'running' : '');
    if (state.state === 'running') { setTimeout(poll, 1000); return; }
    setBusy(false);
    if (state.state === 'done') { job = state; await showResult(); status(''); }
    else if (state.state === 'error') showPage('configure');
  } catch (error) { setBusy(false); fail(error); }
}
function stat(value, label) {
  const el = document.createElement('div'); el.className = 'stat';
  const number = document.createElement('b'); number.textContent = value == null ? '—' : value.toLocaleString();
  const text = document.createElement('span'); text.textContent = label; el.append(number, text); return el;
}
function showDataAudit(report) {
  var translateServer = globalThis.ChemGridMapI18n?.translateServer || (text => text);
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $('audit-body').replaceChildren();
  const note = document.createElement('p');
  note.textContent = job.options.input_format === 'chembl' ? t("d.All retained records are summarized by median pChEMBL. A wide record range is retained, not converted to Medium or silently removed. Aggregation does not harmonize assay conditions. An entry without a wide record range is not proof of experimental agreement.") : job.options.compound_column ? t("d.Retained activity records are converted to pActivity and summarized by compound. Evidence units describe the supplied document, batch or run identifiers.") : t("d.This map uses the supplied molecule-level annotations. It does not create an assay-record provenance trail.");
  $('audit-body').append(note);
  if (report.document_coverage) {
    const coverage = document.createElement('p'); coverage.dataset.documentCoverage = 'true';
    coverage.textContent = t(report.document_coverage.kind === 'chembl' ? 'u5.coverageChembl' : 'u5.coverageReference', report.document_coverage);
    $('audit-body').append(coverage);
  }
  if (report.source_lookup) {
    const lookup = report.source_lookup;
    const coverage = document.createElement('p'); coverage.dataset.sourceLookup = 'true';
    coverage.textContent = t('u7.coverage', lookup); $('audit-body').append(coverage);
    if (lookup.target_count > 1) {
      const warning = document.createElement('p'); warning.dataset.sourceLookupTargets = 'true';
      warning.textContent = t('u7.multipleTargets', {n: lookup.target_count}); $('audit-body').append(warning);
    }
  }
  if (job.options.assumed_sigma != null) {
    const sigma = document.createElement('p'); sigma.dataset.assumedSigma = 'true';
    const source = report.activity_confidence?.evidence_unit_source || job.options.document_column || t("d.one per compound");
    sigma.textContent = t("d.σ (assumed) {sigma} log units. Evidence units: {source}.", {sigma: format(job.options.assumed_sigma), source: source === 'one unit per compound' ? t('d.one per compound') : source});
    if (job.options.input_format !== 'chembl' && !job.options.document_column) sigma.textContent += t("d. Offset detection needs an evidence-unit column, so it was not run.");
    $('audit-body').append(sigma);
  }
  const generic = report.generic_activity;
  if (generic) {
    const labels = {
      non_exact_relation: t("d.censored or inexact values (relation other than =)"),
      invalid_activity_value: t("d.empty or non-numeric activity values"),
      nonpositive_activity_value: t("d.zero or negative values"),
      unrecognized_activity_unit: t("d.unrecognized units"),
      invalid_smiles: t("d.invalid SMILES"),
      missing_compound_id: t("d.missing compound IDs")
    };
    const counts = generic.excluded_reasons || {non_exact_relation: generic.excluded_non_exact_relation_records};
    const details = Object.entries(labels).filter(([key]) => counts[key] > 0).map(([key, label]) => t('d.{label}: {value}', {label, value: counts[key]})).join(t('d.; '));
    const excluded = document.createElement('p'); excluded.dataset.excludedRecords = 'true';
    excluded.textContent = t("d.Excluded records: {excluded} of {total}.", {excluded: generic.excluded_activity_records, total: generic.input_rows});
    if (details) excluded.textContent = t('d.{first} {second}', {first: excluded.textContent, second: t('d.{details}.', {details: details[0].toUpperCase() + details.slice(1)})});
    $('audit-body').append(excluded);
  }
  const exclusions = report.exclusion_counts || {};
  if (Object.keys(exclusions).length) {
    const reasons = document.createElement('p'); reasons.textContent = t("d.Exclusion reasons (may overlap): ") + Object.entries(exclusions).map(([k,v]) => t('d.{label}: {value}', {label: t('d.' + k), value: v})).join(t('d.; ')); $('audit-body').append(reasons);
  }
  if (report.query_molecules) {
    const section = document.createElement('section'), heading = document.createElement('h3');
    heading.textContent = t("u8.newMolecules"); section.append(heading);
    const reference = document.createElement('p'); reference.dataset.newMoleculesReference = 'true';
    reference.textContent = t('u8d.auditReference', {...referenceMapDetails(), filename: (typeof newMolecules !== 'undefined' && newMolecules?.input?.filename) || queryInputInfo?.filename || '—'});
    section.append(reference);
    for (const [key, label] of [['mapped',t('u8.mapped')],['invalid_smiles',t('u8.invalid_smiles')],['duplicate_query',t('u8.duplicate_query')],['in_reference',t('u8.in_reference')]]) {
      const count = document.createElement('p'); count.dataset.queryStatus = key;
      count.dataset.queryCount = report.query_molecules[key]; count.textContent = t('d.{label}: {value}', {label, value: report.query_molecules[key]}); section.append(count);
    }
    const matches = report.query_reference_matches || [];
    if (matches.length) {
      const heading = document.createElement('h4'), container = document.createElement('div');
      heading.textContent = t('u6.referenceMatches');
      container.className = 'table-scroll query-reference-matches';
      const visible = matches.slice(0, 50);
      table(visible.map(match => {
        const title = match.molecule_id, identity = match.parent_molecule_id || match.reference_id;
        return {
          [t('u6.queryId')]: match.query_id,
          [t('u6.measuredId')]: identity && title !== identity ? t('d.{title} (parent {identity})', {title, identity}) : title,
          [t('u6.sourceCount')]: match.source_count,
          [t('u6.linkedSourceCount')]: match.linked_source_count,
          ...((match.user_median ?? match.own_median) != null ? {[t('u8.ownMedian')]:match.user_median ?? match.own_median,[t('u8.referenceMedian')]:match.reference_median,[t('u8.difference')]:format(match.median_difference,true)} : {})
        };
      }), container);
      container.querySelectorAll('tbody tr').forEach((row, index) => {
        row.tabIndex = 0; row.setAttribute('role', 'button');
        row.dataset.queryId = visible[index].query_id;
        const locate = () => { $('audit-dialog').close(); inspect('id=' + encodeURIComponent(visible[index].query_id)); };
        row.onclick = locate;
        row.onkeydown = event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); locate(); } };
      });
      section.append(heading, container);
      if (matches.length > visible.length) {
        const remaining = document.createElement('p'); remaining.className = 'micro';
        remaining.textContent = t('u6.remainingMatches', {count: matches.length - visible.length});
        section.append(remaining);
      }
    }
    const added = report.query_molecules;
    if (added.input_rows != null) { const p=document.createElement('p');p.textContent=t('u8.inputRows',{count:added.input_rows});section.append(p); }
    for (const category of ['active','medium','inactive','mixed','insufficient','none']) {
      if (!added.neighbor_evidence_counts) break;
      const p=document.createElement('p');p.dataset.evidenceCategory=category;p.textContent=t('d.{label}: {value}',{label:t('u8.'+category),value:added.neighbor_evidence_counts[category] || 0});section.append(p);
    }
    if (added.reference_support_assessed === false) {const p=document.createElement('p');p.textContent=t('u8.supportNotAssessed');section.append(p);}
    const differences = matches.map(match => match.median_difference).filter(value => value != null).sort((a,b)=>a-b);
    if (differences.length) {
      const middle=Math.floor(differences.length/2),median=differences.length%2 ? differences[middle] : (differences[middle-1]+differences[middle])/2;
      const p=document.createElement('p');p.textContent=t('u8.differenceSummary',{median:format(median,true),minimum:format(differences[0],true),maximum:format(differences[differences.length-1],true)});section.append(p);
    }
    const excludedReasons = Object.entries(added.excluded_reasons || {}).filter(([,count]) => count > 0);
    if (excludedReasons.length) {
      const h=document.createElement('h4'),container=document.createElement('div');h.textContent=t('u8.excluded');
      container.className='new-excluded-reasons';
      table(excludedReasons.map(([reason,count]) => ({[t('u8.reason')]:t('d.'+reason),[t('d.Source records')]:count})),container);section.append(h,container);
    }
    if (added.excluded_records?.length) {
      const container=document.createElement('div');container.className='table-scroll new-excluded-records';
      const columns=added.excluded_record_columns || {}, defaults=added.excluded_record_defaults || {};
      table(added.excluded_records.slice(0,50).map(record => ({
        [t('d.Compound ID')]:record[columns.compound_id],
        [t('d.Evidence units')]:record[columns.evidence_unit],
        [t('u8.value')]:record[columns.value],
        [t('u8.unit')]:columns.unit ? record[columns.unit] : defaults.unit,
        [t('u8.relation')]:columns.relation ? record[columns.relation] : defaults.relation,
        [t('u8.reason')]:t('d.'+record.exclusion_reason)
      })),container);section.append(container);
      if (added.excluded_records.length > 50) {
        const more=document.createElement('p');more.className='micro';
        more.textContent=t('u8.excludedMore',{count:added.excluded_records.length-50});section.append(more);
      }
    }
    if (added.no_retained_activity_compounds?.length) { const p=document.createElement('p');p.textContent=t('u8.noRetained',{ids:added.no_retained_activity_compounds.join(globalThis.ChemGridMapI18n?.getLanguage() === 'zh' ? '、' : ', ')});section.append(p); }
    if (added.document_pair_offsets?.length) {
      const h=document.createElement('h4'),container=document.createElement('div');h.textContent=t('d.Evidence units and offsets');container.className='table-scroll new-offset-pairs';
      table(added.document_pair_offsets.map(pair=>({Pair:t('d.{a} vs {b}',{a:pair.document_a,b:pair.document_b}),'Shared entries':pair.shared_entries,'Median difference (B − A)':format(pair.median_offset,true,2),'Sign agreement':`${Math.round(100*pair.sign_agreement)}%`,Flagged:pair.systematic_offset?t('d.Yes'):t('d.No')})),container);section.append(h,container);
    }
    for (const message of report.query_molecules.warnings || []) {
      const warning = document.createElement('p'); warning.className = 'new-molecule-warning';
      warning.textContent = newMoleculeWarning(message); section.append(warning);
    }
    $('audit-body').append(section);
  }
  const confidence = report.activity_confidence;
  if (confidence?.document_pair_count != null && confidence.document_offset_detection !== 'not applicable' && (!generic || job.options.document_column)) {
    const section = document.createElement('section'), heading = document.createElement('h3');
    heading.textContent = t("d.Evidence units and offsets"); section.append(heading);
    const description = document.createElement('p');
    const source = job.options.input_format === 'chembl' ? t('d.document') : (globalThis.ChemGridMapI18n?.getLanguage() === 'zh' ? ' ' : '') + t('d.{column} column', {column: generic?.document_column || job.options.document_column});
    description.textContent = t("d.Evidence units come from {source}; {pairs} pairs with at least {minimum} shared entries were evaluated, {flagged} pairs were flagged, using |median difference| ≥ {threshold} log units and sign agreement ≥ {agreement}%.", {source, pairs: confidence.document_pair_count, minimum: confidence.offset_min_shared, flagged: confidence.discordant_document_pair_count, threshold: format(confidence.offset_threshold), agreement: Math.round(100 * confidence.offset_sign_agreement)});
    section.append(description);
    const pairs = confidence.document_offsets_top10 || [], container = document.createElement('div');
    container.className = 'table-scroll offset-pairs';
    table(pairs.map(pair => ({
      Pair: t('d.{a} vs {b}', {a: pair.document_a, b: pair.document_b}),
      'Shared entries': pair.shared_entries,
      'Median difference (B − A)': `${pair.median_offset >= 0 ? '+' : '−'}${Math.abs(pair.median_offset).toFixed(2)}`,
      'Sign agreement': `${Math.round(100 * pair.sign_agreement)}%`,
      Flagged: pair.systematic_offset ? t("d.Yes") : t("d.No"),
      'Entries touched': pair.touched_entry_count
    })), container);
    container.querySelectorAll('tbody tr').forEach((row, index) => { row.dataset.offsetFlagged = String(pairs[index].systematic_offset); });
    section.append(container);
    const display = document.createElement('p');
    display.textContent = (confidence.offset_display || job.options.offset_display || 'uncertain') === 'uncertain'
      ? t("d.Entries with records from a flagged pair are shown as uncertain.")
      : t("d.Display unchanged: entries with records from a flagged pair are marked ◣ on Activity class support.");
    const full = document.createElement('p'); full.textContent = t("d.Full table: map_document_offsets.csv (Export).");
    section.append(display, full); $('audit-body').append(section);
  }
  const warnings = document.createElement('ul');
  for (const text of report.warnings || []) { const item = document.createElement('li'); item.textContent = translateServer(text); warnings.append(item); }
  $('audit-body').append(warnings); $('audit').open = true;
  $('audit-count').textContent = (report.warnings || []).length;
  $('audit-count').hidden = !warnings.children.length;
}
async function showResult(languageOnly = false) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  if (!languageOnly) restoreOptions(job.options);
  const report = job.report;
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const newCount = report.query_molecules?.mapped || 0;
  const stats = [stat(report.input_rows, t("d.Input records")), stat(job.molecules-newCount, t(newCount ? 'u8.referenceMolecules' : 'd.Mapped molecules'))];
  if (report.query_molecules?.mapped) stats.push(stat(report.query_molecules.mapped, t('u8.newMolecules')));
  stats.push(stat(report.retained_activity_records, t("d.Retained records")), stat(report.repeated_conflicted_molecules, t("d.Wide record ranges")));
  $('summary').replaceChildren(...stats);
  $('summary').hidden = false; $('audit').hidden = false; if (!languageOnly) $('preview').open = false;
  const auditOpen = $('audit').open;
  if (!languageOnly) {
  newMolecules = await getJSON('/api/new-molecules');
  const svgFile = (newMolecules?.files || []).find(name => name.endsWith('map_with_new_molecules.svg')) || job.files.find(name => name === 'map.svg');
  gridSVG = await (await api('/api/file?name=' + encodeURIComponent(svgFile))).text();
  points = await getJSON('/api/points');
  $('measurement-filter-control').hidden = !points.some(point => point.is_new || point.measurement_status === 'not_measured');
  $('color-mode').querySelector('[value=neighbor_evidence]').disabled = newMolecules?.report?.input_rows == null && !points.some(point => point.is_new || point.measurement_status === 'not_measured');
  showResult.grouping = await getJSON('/api/groups');
  }
  showDataAudit(report);
  if (languageOnly) $('audit').open = auditOpen;
  const grouping = showResult.grouping; structureGroups = grouping.groups;
  $('group-controls').hidden = false;
  document.querySelectorAll('[data-group-only]').forEach(el => el.hidden = !structureGroups.length);
  const structureOption = $('color-mode').querySelector('[value=structure]');
  structureOption.disabled = structureOption.hidden = !structureGroups.length;
  if (structureOption.disabled && $('color-mode').value === 'structure') $('color-mode').value = 'activity';
  $('group-select').replaceChildren();
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const none = document.createElement('option'); none.value = ''; none.textContent = t("d.None"); $('group-select').append(none);
  for (const group of structureGroups) {
    const option = document.createElement('option'); option.value = group.structure_group;
    option.textContent = t("d.{id} | {members} | core {covered}/{size}", {id: group.structure_group, members: countLabel(group.size, 'member'), covered: group.core_coverage_count, size: group.size});
    $('group-select').append(option);
  }
  $('group-select').value = highlightedGroup;
  $('group-method').textContent = t("d.{method} groups are defined before projection. Colors may repeat; group IDs are authoritative.", {method: t('d.' + grouping.protocol.method)});
  $('neighbor-key').querySelector('.representation-only').textContent =
    job.options.representation === 'morgan' ? t("d.Fingerprint only") : t("d.Representation only");
  if (!languageOnly) { updateNavigation(); showPage('explore'); }
  $('run-caption').textContent = [referenceMapDetails().description, job.options.representation === 'morgan' ? t("d.Morgan fingerprint") : t("d.Molecular descriptors"), job.options.projection.toUpperCase(), t("d.{value}% occupancy", {value: Math.round(job.options.grid_occupancy * 100)})].join(' · ');
  renderReferenceMap();
  if (!languageOnly) { view = 'grid'; setView('grid'); }
  $('download-all').hidden = false;
  $('downloads').replaceChildren();
  $('figure-downloads').replaceChildren();
  for (const name of [...new Set([...job.files, ...(newMolecules?.files || [])])].filter(name => !['input.csv','chemgridmap_results.zip'].includes(name))) {
    const button = document.createElement('button');
    const figure = /^map(?:_(?:class_support|neighbor_recall|structure_group))?\.(svg|png|pdf)$/.test(name);
    if (figure) {
      const extension = name.split('.').pop(), mode = name.slice(0, -(extension.length + 1));
      const titles = {map: t("d.Observed activity"), map_class_support: t("d.Activity class support"),
        map_neighbor_recall: t("d.Neighbor reliability"), map_structure_group: t("d.Structure group")};
      button.textContent = `${titles[mode]} · ${extension.toUpperCase()}`;
      const caption = document.createElement('small'); caption.textContent = extension === 'png' ? t("d.Raster image") : t("d.Vector figure"); button.append(caption);
    } else button.textContent = name;
    button.addEventListener('click', () => download('/api/file?name=' + encodeURIComponent(name), name));
    $(figure ? 'figure-downloads' : 'downloads').append(button);
  }
  $('saved-location').textContent = t("d.Saved locally: ") + job.output_dir;
  table(Object.entries(job.metrics).map(([metric,value]) => ({Metric:metric,Value:value})), $('metrics-table')); $('metrics').hidden = false;
  const hasSupport = points.some(p => p.class_status != null), hasReliability = points.some(p => p.neighbor_recall != null);
  $('color-mode').querySelector('[value=support]').disabled = !hasSupport;
  $('color-mode').querySelector('[value=reliability]').disabled = !hasReliability;
  if ($('color-mode').selectedOptions[0]?.disabled) $('color-mode').value = 'activity';
  renderLegend($('color-mode').value);
}
function projectionSVG() {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg'); svg.setAttribute('viewBox','0 0 900 900'); svg.setAttribute('width','900'); svg.setAttribute('height','900');
  for (const point of points.filter(point => !point.is_new && point.measurement_status !== 'not_measured')) {
    const circle = document.createElementNS(ns,'circle');
    circle.setAttribute('cx',40 + point.projection_x * 820); circle.setAttribute('cy',860 - point.projection_y * 820);
    circle.setAttribute('r','4'); circle.setAttribute('fill',colors[point.activity_class] || '#b8c2cc');
    circle.dataset.gridRow = point.grid_row; circle.dataset.gridCol = point.grid_col;
    const title = document.createElementNS(ns,'title'); title.dataset.originalText = `${point.molecule_id || point.molecule_identity_key || 'Molecule'} · ${point.activity_class}`; title.textContent = title.dataset.originalText; circle.append(title); svg.append(circle);
  }
  return svg;
}
function setView(name) {
  view = name;
  for (const key of ['grid','projection']) { $(key+'-tab').classList.toggle('active',key === name); $(key+'-tab').setAttribute('aria-selected',key === name); }
  if (!job) return;
  $('empty-map').hidden = true;
  if (name === 'grid') $('map-stage').innerHTML = gridSVG;
  else $('map-stage').replaceChildren(projectionSVG());
  const svg = $('map-stage').querySelector('svg'), box = svg.viewBox.baseVal;
  // The browser legend follows the selected layer; exports keep their own legend.
  if (points.some(point => point.is_new || point.measurement_status === 'not_measured')) svg.querySelector('[data-color-legend]')?.remove();
  tileWidth = Number(svg.querySelector('g[data-grid-row] > rect')?.getAttribute('width') || 100);
  hasStructures = !!svg.querySelector('g[data-grid-row] > g');
  mapCells = Array.from(svg.querySelectorAll('[data-grid-row]'));
  drawingSize = {w:box.width,h:box.height}; svg.setAttribute('width',box.width); svg.setAttribute('height',box.height);
  colorMap(); fit(); highlight(); translateSVGText($('map-stage'));
}
function colorMap() {
  const mode = $('color-mode').value;
  const lookup = new Map(points.map(p => [`${p.grid_row},${p.grid_col}`, p]));
  $('map-stage').querySelectorAll('[data-grid-row]').forEach(el => {
    const point = lookup.get(`${el.dataset.gridRow},${el.dataset.gridCol}`);
    if (!point) return;
    const mark = el.tagName === 'circle' ? el : el.querySelector(':scope > rect');
    if (!mark) return;
    if (!mark.dataset.activityFill) mark.dataset.activityFill = mark.getAttribute('fill');
    if (!('activityStroke' in mark.dataset)) {
      mark.dataset.activityStroke = mark.getAttribute('stroke') || '';
      mark.dataset.activityStrokeWidth = mark.getAttribute('stroke-width') || '';
    }
    const unmeasured = point.measurement_status === 'not_measured', isNew = point.is_new || unmeasured;
    if (point.measurement_status) el.dataset.measurementStatus = point.measurement_status;
    el.classList.toggle('unmeasured-outline', unmeasured && mode === 'activity');
    el.classList.toggle('new-molecule-outline', !!isNew && !(unmeasured && mode === 'activity') && mode !== 'neighbor_evidence');
    el.classList.toggle('reference-evidence-faded', !isNew && mode === 'neighbor_evidence');
    el.classList.toggle('new-evidence-outline', !!isNew && mode === 'neighbor_evidence');
    let fill = mark.dataset.activityFill;
    if (mode === 'structure') fill = point.structure_group_color || '#B8C2CC';
    else if (mode === 'reliability') fill = reliabilityColor(point.neighbor_recall);
    else if (unmeasured) fill = '#FFFFFF';
    if (isNew && ['support','reliability'].includes(mode)) fill = '#dfe4e9';
    else if (isNew && ['activity','neighbor_evidence'].includes(mode)) fill = unmeasured ? '#FFFFFF' : colors[point.activity_class] || '#b8c2cc';
    else if (mode === 'support' && point.class_status === 'uncertain') fill = fadedColors[point.activity_class] || fill;
    mark.setAttribute('fill', fill);
    if (isNew && mode === 'neighbor_evidence') {
      const style = ({active:['#2a9d8f','',3],medium:['#e6a23c','',3],inactive:['#c44e52','',3],mixed:['#7a5195','6 4',3],insufficient:['#8a99a6','1 4',3],none:['#b8c2cc','',1.5]})[point.neighbor_evidence || 'none'];
      mark.style.setProperty('--evidence-color',style[0]); mark.style.setProperty('--evidence-width',String(style[2]));
      mark.style.setProperty('--evidence-dash',style[1] || 'none');
      mark.setAttribute('stroke',style[0]);
      mark.setAttribute('stroke-dasharray',style[1]);
    } else if (unmeasured && mode === 'activity') {
      mark.setAttribute('stroke', '#8A99A6'); mark.setAttribute('stroke-width', '1');
    } else {
      for (const [attribute, value] of [['stroke',mark.dataset.activityStroke],['stroke-width',mark.dataset.activityStrokeWidth]]) {
        if (value) mark.setAttribute(attribute, value); else mark.removeAttribute(attribute);
      }
    }
    if (mode !== 'neighbor_evidence') mark.removeAttribute('stroke-dasharray');
    let offsetMarker = el.querySelector(':scope > [data-offset-pair-marker]');
    if (!isNew && mode === 'support' && point.offset_pair_touched) {
      if (!offsetMarker) {
        const x = Number(mark.getAttribute('x') || 0), y = Number(mark.getAttribute('y') || 0);
        const width = Number(mark.getAttribute('width') || 100), height = Number(mark.getAttribute('height') || 100);
        const corner = width * 0.3;
        offsetMarker = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        offsetMarker.dataset.offsetPairMarker = 'true';
        offsetMarker.setAttribute('d', `M${x},${y + height - corner} L${x},${y + height} L${x + corner},${y + height} Z`);
        offsetMarker.setAttribute('fill', '#172536');
        el.append(offsetMarker);
      }
    } else if (offsetMarker) offsetMarker.remove();
    let conflictMarker = el.querySelector(':scope > [data-conflict-marker]');
    if ((point.wide_record_range ?? point.is_conflicted) && (mode === 'activity' || !isNew && mode === 'support')) {
      if (!conflictMarker) {
        const x = Number(mark.getAttribute('x') || 0), y = Number(mark.getAttribute('y') || 0);
        const width = Number(mark.getAttribute('width') || 100);
        const corner = width * 0.3;
        conflictMarker = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        conflictMarker.dataset.conflictMarker = 'true';
        conflictMarker.setAttribute('d', `M${x + width - corner},${y} L${x + width},${y} L${x + width},${y + corner} Z`);
        conflictMarker.setAttribute('fill', '#172536');
        el.append(conflictMarker);
      }
    } else if (conflictMarker) conflictMarker.remove();
  });
  renderLegend(mode);
}
function legendItem(color, text) {
  const item = document.createElement('span'), swatch = document.createElement('i');
  swatch.style.background = color; item.append(swatch, document.createTextNode(text)); return item;
}
function renderLegend(mode) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  renderColorLegend(mode);
  const hasQueries = points.some(point => point.is_new || point.measurement_status === 'not_measured');
  if (hasQueries && ['activity','support','reliability'].includes(mode)) {
    const items = mode === 'activity' ? [[t('u8.newUnmeasured'),'#FFFFFF','#8A99A6'],[t('u8.newMeasured'),colors.active,'#173e34']] : [[t('u8.notInLayer'),'#dfe4e9','#173e34']];
    for (const [label,fill,border] of items) { const item = legendItem(fill,label + ' ' + t('u8d.seeNeighborEvidence')); item.firstElementChild.style.border = '1.5px solid '+border; $('legend').append(item); }
  }
  const group = structureGroups.find(g => g.structure_group === highlightedGroup);
  if (group) {
    const note = document.createElement('span'); note.className = 'group-legend';
    const members = points.filter(point => point.structure_group === group.structure_group);
    const queryMembers = members.filter(point => point.is_new || point.measurement_status === 'not_measured').length;
    note.textContent = t("d.Outlined: {id}, {members}{unmeasured}. Others faded.", {id: group.structure_group, members: countLabel(members.length, 'molecule'), unmeasured: queryMembers ? t('u8.newGroupCount', {count: queryMembers}) : ''});
    $('legend').append(note); $('legend').hidden = false;
  }
  if (hasQueries) {
    const note = document.createElement('span'); note.className = 'measurement-legend';
    note.textContent = t("d.Show: {filter}{suffix}", {filter: $('measurement-filter').selectedOptions[0].textContent, suffix: $('measurement-filter').value === 'all' ? t('d..') : t('d.. Other molecules faded.')});
    $('legend').append(note); $('legend').hidden = false;
  }
}
function renderColorLegend(mode) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $('legend').replaceChildren(); $('legend').hidden = mode === 'structure';
  if (!job || mode === 'structure') return;
  if (mode === 'neighbor_evidence') {
    for (const [key,color,border,width] of [['active','#2a9d8f','solid',3],['medium','#e6a23c','solid',3],['inactive','#c44e52','solid',3],['mixed','#7a5195','dashed',3],['insufficient','#8a99a6','dotted',3],['none','#b8c2cc','solid',1.5]]) {
      const item = legendItem('#fff',t('u8.'+key)); item.dataset.evidenceCategory = key; item.firstElementChild.style.border = `${width}px ${border} ${color}`; $('legend').append(item);
    }
    const note = document.createElement('span'); note.textContent = t('u8.referenceFaded'); $('legend').append(note); return;
  }
  const lower = job.options.lower_threshold, upper = job.options.upper_threshold;
  if (job.options.input_format === 'molecule' && !job.options.compound_column && mode !== 'reliability') {
    $('legend').textContent = t("d.Colors follow supplied class labels, or the selected pActivity thresholds; gray = unlabelled.");
    return;
  }
  if (mode === 'reliability') {
    reliabilityBins.forEach(([low, color], i) => $('legend').append(legendItem(color,
      `${Math.round(low * 100)}–${Math.round((reliabilityBins[i + 1]?.[0] ?? 1) * 100)}%`)));
    const term = job.options.representation === 'morgan' ? 'fingerprint' : 'representation';
    $('legend').prepend(document.createTextNode(t("d.Share of {term} neighbors that are also grid neighbors:", {term: t('d.' + term)})));
    return;
  }
  for (const [label, text] of [['inactive', t("d.Inactive ≤{lower}", {lower})], ['medium', t("d.Medium ({lower}, {upper})", {lower, upper})], ['active', t("d.Active ≥{upper}", {upper})]]) {
    $('legend').append(legendItem(colors[label], text));
    if (mode === 'support') $('legend').append(legendItem(fadedColors[label], t('d.uncertain')));
  }
  if (mode === 'support') {
    // Jobs created before 99b57ad always showed touched entries as uncertain.
    const offsetDisplay = job.options.offset_display || 'uncertain';
    const suffix = offsetDisplay === 'uncertain' ? t("d., or an offset document pair") : '';
    $('legend').append(document.createTextNode(t("d.Pale: uncertain (class support below {level}{suffix}).", {level: job.options.support_level ?? .8, suffix})));
    if (points.some(point => point.offset_pair_touched)) {
      const offset = document.createElement('span');
      offset.textContent = t("d.◣ offset document pair{suffix}", {suffix: offsetDisplay === 'uncertain' ? '' : t('d. (display unchanged)')});
      $('legend').append(offset);
    }
  }
  const flag = document.createElement('span'); flag.textContent = t("d.◥ Wide record range"); $('legend').append(flag);
}
function applyTransform() {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $('map-stage').style.transform = `translate(${transform.x}px,${transform.y}px) scale(${transform.s})`;
  $('map-stage').style.setProperty('--zoom-inverse', String(1 / transform.s));
  // Suppress unreadable bonds only in the browser; exported SVGs remain untouched.
  const overview = view === 'grid' && tileWidth * transform.s < 24;
  $('map-viewport').classList.toggle('color-overview', overview);
  $('view-hint').textContent = view === 'projection' ? t("d.2D projection · Click a point to inspect") + (points.some(point => point.is_new || point.measurement_status === 'not_measured') ? ' · '+t('u8.projectionNote') : '') :
    overview && hasStructures ? t("d.Color overview · Zoom in for structures · Click a cell to inspect") :
    t("d.Click a cell to inspect · Scroll to zoom · Drag to pan");
}
function fit() {
  if (!job || page !== 'explore') return;
  const width = $('map-viewport').clientWidth, height = $('map-viewport').clientHeight;
  const s = Math.min((width - 35) / drawingSize.w,(height - 35) / drawingSize.h);
  transform = {s,x:(width-drawingSize.w*s)/2,y:(height-drawingSize.h*s)/2}; applyTransform();
}
function zoom(factor,x=$('map-viewport').clientWidth/2,y=$('map-viewport').clientHeight/2) {
  if (!job) return;
  const s = Math.min(15,Math.max(.005,transform.s*factor)), ratio = s/transform.s;
  transform = {s,x:x-(x-transform.x)*ratio,y:y-(y-transform.y)*ratio}; applyTransform();
}
function highlight() {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $('map-stage').querySelector('[data-map-outlines]')?.remove();
  const groups = new Map(points.map(p => [`${p.grid_row},${p.grid_col}`, p.structure_group]));
  const measurements = new Map(points.map(p => [`${p.grid_row},${p.grid_col}`, p.is_new || p.measurement_status === 'not_measured' ? 'new' : 'reference']));
  const newSelection = currentCell && points.find(p => (p.is_new || p.measurement_status === 'not_measured') && p.grid_row === currentCell.row && p.grid_col === currentCell.col);
  const relationFor = el => newSelection || $('show-neighbor-links').checked ? neighborRelations.get(`${el.dataset.gridRow},${el.dataset.gridCol}`) : null;
  const measurementFilter = $('measurement-filter').value;
  const layers = [[], [], [], []];
  const outlines = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  outlines.dataset.mapOutlines = 'true'; outlines.setAttribute('pointer-events', 'none');
  mapCells.forEach(el => {
    const member = !!highlightedGroup && groups.get(`${el.dataset.gridRow},${el.dataset.gridCol}`) === highlightedGroup;
    const unmeasured = measurements.get(`${el.dataset.gridRow},${el.dataset.gridCol}`) === 'new';
    const measurementFaded = measurementFilter !== 'all' && measurements.get(`${el.dataset.gridRow},${el.dataset.gridCol}`) !== measurementFilter;
    el.classList.toggle('group-highlight', member);
    el.classList.toggle('measurement-faded', measurementFaded);
    el.classList.toggle('group-faded', !!highlightedGroup && !member || measurementFaded);
    const selected = currentCell && Number(el.dataset.gridRow) === currentCell.row && Number(el.dataset.gridCol) === currentCell.col;
    el.classList.toggle('selected-cell',!!selected);
    const relation = relationFor(el);
    for (const category of ['shared','grid_only','representation_only']) el.classList.toggle('neighbor-' + category, relation === category);
    if (el.tagName === 'circle') {
      const queryOutline = unmeasured && ['activity','support'].includes($('color-mode').value);
      el.setAttribute('stroke', selected ? '#183e92' : queryOutline ? '#8A99A6' : 'none');
      el.setAttribute('stroke-width',selected ? '3' : queryOutline ? '1' : '3');
    }
    if (!el.classList.length) el.removeAttribute('class');
    layers[selected ? 3 : relation ? 2 : member ? 1 : 0].push(el);
  });
  // Paint group outlines above adjacent cells, then neighbors and the selection.
  layers.flat().forEach(el => {
    el.parentNode.append(el);
    if (!el.classList.contains('group-faded') && !el.classList.contains('reference-evidence-faded')) return;
    const relation = relationFor(el);
    if (!el.classList.contains('selected-cell') && !relation) return;
    // Keep inspection outlines opaque while the entire underlying cell is faded.
    const outline = (el.tagName === 'circle' ? el : el.querySelector(':scope > rect')).cloneNode(false);
    for (const attribute of ['class', 'data-grid-row', 'data-grid-col']) outline.removeAttribute(attribute);
    outline.setAttribute('fill', 'none'); outline.dataset.outline = el.classList.contains('selected-cell') ? 'selected' : relation;
    // 保留参照 SVG 包装层的平移；选框统一绘制在根层。
    const transforms = [];
    for (let parent = el; parent && parent.tagName.toLowerCase() !== 'svg'; parent = parent.parentNode) {
      const value = parent.getAttribute('transform'); if (value) transforms.unshift(value);
    }
    if (outline.getAttribute('transform')) transforms.push(outline.getAttribute('transform'));
    if (transforms.length) outline.setAttribute('transform',transforms.join(' '));
    outlines.append(outline);
  });
  const svg = $('map-stage').querySelector('svg');
  if (newSelection && view === 'grid') {
    const nearestId = newSelection.nearest_reference_id || newSelection.nearest_measured_id;
    const nearest = points.find(p => !p.is_new && p.measurement_status !== 'not_measured' && (p.molecule_identity_key || p.molecule_id) === nearestId);
    const cellFor = point => point && mapCells.find(el => Number(el.dataset.gridRow) === point.grid_row && Number(el.dataset.gridCol) === point.grid_col);
    const from = cellFor(newSelection), to = cellFor(nearest);
    if (from && to) {
      // 把格子中心变换到根 SVG 坐标，保留参照包装层的平移。
      const center = cell => {
        const rect = cell.querySelector(':scope > rect'), box = rect.getBBox(), point = svg.createSVGPoint();
        point.x = box.x+box.width/2; point.y = box.y+box.height/2;
        return point.matrixTransform(svg.getCTM().inverse().multiply(rect.getCTM()));
      };
      const start=center(from), end=center(to), line=document.createElementNS('http://www.w3.org/2000/svg','line');
      line.dataset.newReferenceLink='true';
      for (const [key,value] of Object.entries({x1:start.x,y1:start.y,x2:end.x,y2:end.y})) line.setAttribute(key,value);
      outlines.prepend(line);
    }
  }
  if (outlines.children.length) svg.append(outlines);
  const selectedGroup = currentCell && groups.get(`${currentCell.row},${currentCell.col}`);
  $('selected-group-highlight').textContent = selectedGroup && selectedGroup === highlightedGroup ? t("d.Clear highlight") : t("d.Highlight on map");
}
function coreRule(group) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  return group.core_rule === 'coverage'
    ? t("d.largest scaffold or ring fragment matching at least {coverage}% of members ({status}).", {coverage: Math.round(100 * group.core_coverage_target), status: group.core_target_met ? t('d.target met') : t('d.target not met; highest-coverage core shown')})
    : group.core_rule === 'fmcs'
      ? t("d.FMCS matching at least {coverage}% of members{timeout}.", {coverage: Math.round(100 * group.core_coverage_target), timeout: group.fmcs_timed_out ? t('d.; the search timed out, so this partial core may vary with machine speed and load') : ''})
      : t("d.most frequent ring scaffold.");
}
function countLabel(count, noun) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  if (noun === 'reference member') return t(count === 1 ? 'u8.referenceMember' : 'u8.referenceMembers',{count});
  return t(`d.{count} ${noun}${count === 1 ? '' : 's'}`, {count});
}
async function selectGroup(groupId, cached = null) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const request = cached ? groupRequestNumber : ++groupRequestNumber;
  const group = structureGroups.find(g => g.structure_group === groupId);
  highlightedGroup = group ? groupId : ''; $('group-select').value = highlightedGroup;
  $('group-core-card').hidden = !group; $('group-core-svg').replaceChildren();
  highlight(); renderLegend($('color-mode').value); updateNavigation();
  if (!group) return;
  const queryMembers = points.filter(point => (point.is_new || point.measurement_status === 'not_measured') && point.structure_group === groupId).length;
  const members = countLabel(group.size, queryMembers ? 'reference member' : 'member');
  $('group-core-title').textContent = t("d.{id} · {members}{unmeasured}", {id: group.structure_group, members, unmeasured: queryMembers ? t('u8.newGroupAdded', {count: queryMembers}) : ''});
  $('group-core-coverage').textContent = group.core_smiles
    ? t("d.Core in {covered} of {members} · exact scaffold in {exact}", {covered: group.core_coverage_count, members, exact: group.core_exact_scaffold_count})
    : t("d.{members} · no ring scaffold", {members});
  $('group-core-rule').textContent = coreRule(group);
  $('group-core-svg').textContent = group.core_smiles ? t("d.Loading core...") : t("d.No ring scaffold");
  try {
    const result = cached || await getJSON('/api/group?group=' + encodeURIComponent(groupId));
    if (request !== groupRequestNumber) return;
    selectGroup.lastResult = result;
    if (result.core_svg) $('group-core-svg').innerHTML = result.core_svg;
  } catch (error) { if (request === groupRequestNumber) { $('group-core-svg').textContent = t("d.Core unavailable"); fail(error); } }
}
function format(value, signed = false, precision = null) {
  if (value == null) return '—';
  if (typeof value !== 'number') return String(value);
  if (signed) return (value < 0 ? '−' : '+') + (precision == null ? Number(Math.abs(value).toFixed(3)).toString() : Math.abs(value).toFixed(precision));
  return Number(value.toFixed(3)).toString();
}
function ordinal(n) {
  const lastTwo = n % 100;
  return `${n}${lastTwo >= 11 && lastTwo <= 13 ? 'th' : ({1: 'st', 2: 'nd', 3: 'rd'})[n % 10] || 'th'}`;
}
function percentage(value) {
  if (value == null || Number.isNaN(value)) return '—';
  if (value >= 0.995) return '> 99%';
  if (value < 0.005) return '< 1%';
  return `${Math.round(100 * value)}%`;
}
function measuredLabel(identity) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  // Neighbor identities remain the click target; the visible ID matches the inspector.
  const point = points.find(point => (point.molecule_identity_key || point.molecule_id) === identity);
  const title = point?.molecule_id || identity;
  return title === identity ? title : t("d.{title} (parent {identity})", {title, identity});
}
function showQuerySelection(detail) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const row = detail.molecule;
  const stats = document.createElement('dl');
  const fields = [
    [t('u8.referenceId'), `${measuredLabel(row.nearest_measured_id)} · Tanimoto ${format(row.nearest_measured_similarity)}`],
    [t("d.Structure group"), row.structure_group || t('u8.noStructureGroup', {similarity: format(detail.group_similarity ?? job.options.group_similarity)})],
    [t("d.Contains the group core"), row.contains_group_core ? t("d.Yes") : t("d.No")]
  ];
  for (const [label, value] of fields) {
    const dt = document.createElement('dt'), dd = document.createElement('dd');
    dt.textContent = label; dd.textContent = value; stats.append(dt, dd);
  }
  $('new-placement').textContent = row.placement_distance == null ? '' : t('u8.placement',{id:measuredLabel(row.nearest_reference_id || row.nearest_measured_id),distance:Number(row.placement_distance).toFixed(1)});
  $('new-map-markers').hidden = false;
  $('new-support-note').hidden = row.reference_support_assessed !== false;
  $('new-support-note').textContent = t('u8.supportNotAssessed');
  const extra = [[t('u8.labelFilter'),row.label || '—'],[t('u8.neighborEvidence'),t('u8.'+(row.neighbor_evidence || 'none'))]];
  if (row.similar_reference_count != null) {
    const rawClasses = row.reference_support_assessed === false ? row.evidence_reference_classes : row.supported_reference_classes;
    const classes = Array.isArray(rawClasses) ? rawClasses : String(rawClasses || '').split(/[;,]/).filter(Boolean);
    let countsKey = 'u8.neighborCounts';
    if (row.reference_support_assessed === false) countsKey += 'Unassessed';
    else {
      if (row.supported_reference_count === 0) countsKey += 'NoSupported';
      if (row.similar_reference_count === 1) countsKey += 'One';
    }
    extra.push([t('u8.similarReferenceMolecules'),t(countsKey,{similar:row.similar_reference_count,supported:row.supported_reference_count,classes:classes.map(value => t('d.'+value)).join(' / ') || '—'})]);
  }
  if (row.measurement_status === 'measured') {
    extra.push([t('u8.ownMedian'),format(row.activity_pchembl)],[t('d.{term} range',{term:'pActivity'}),format(row.pchembl_min)+'–'+format(row.pchembl_max)],
      [t('d.Source records'),format(row.n_records)],[t('d.Evidence units'),format(row.n_evidence_units)],[t('d.Activity class'),t('d.'+row.activity_class)]);
    if (row.wide_record_range || row.is_conflicted) extra.push([t('d.◥ Wide record range'),t('d.Yes')]);
    if (row.neighbor_comparison) extra.push([t('u8.comparisonColumn'),t('u8.comparison',{own:t('d.'+row.activity_class),neighbor:t('d.'+row.neighbor_evidence),comparison:t('u8.'+row.neighbor_comparison)})]);
  }
  for (const [label,value] of extra) { const dt=document.createElement('dt'),dd=document.createElement('dd'); dt.textContent=label;dd.textContent=value;stats.append(dt,dd); }
  $('query-summary').replaceChildren(stats);
  $('new-records').hidden = row.measurement_status !== 'measured';
  if (!$('new-records').hidden) table((detail.records || []).map(record => ({
    [t('d.Evidence units')]:record.document_chembl_id ?? record.evidence_unit,
    [t('u8.value')]:record.raw_activity_value ?? record.pchembl_value,
    [t('u8.unit')]:record.raw_activity_unit ?? 'pActivity',
    [t('u8.relation')]:record.record_relation ?? record.relation ?? record.standard_relation ?? '='
  })), $('new-records'));

  const rows = (detail.query_neighbors || []).map(neighbor => ({
    Rank:neighbor.rank, ID:measuredLabel(neighbor.measured_id), Tanimoto:format(neighbor.tanimoto),
    'Display class':neighbor.activity_class, pActivity:format(neighbor.activity_pchembl),
    [t("d.Class support")]:percentage(neighbor.class_confidence), 'Class status':neighbor.class_status,
    [t("d.Evidence units")]:neighbor.n_evidence_units
  }));
  table(rows, $('query-neighbor-table'));
  $('query-neighbor-table').querySelectorAll('tbody tr').forEach((tr, i) => {
    tr.tabIndex = 0; tr.setAttribute('role','button');
    const locate = () => {
      const id = detail.query_neighbors[i].measured_id;
      const reference = points.find(point => !point.is_new && point.measurement_status !== 'not_measured' && (point.molecule_identity_key || point.molecule_id) === id);
      inspect(reference && reference.grid_row != null ? `row=${reference.grid_row}&col=${reference.grid_col}` : 'id='+encodeURIComponent(id));
    };
    tr.onclick = locate; tr.onkeydown = event => { if (event.key === 'Enter') locate(); };
  });
}
function showSources(sources = [], lookup = {}) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  for (const possible of [false, true]) {
    const entries = possible ? lookup.possible_sources || [] : sources;
    const section = $(possible ? 'possible-sources' : 'sources');
    const list = $(possible ? 'possible-source-list' : 'source-list');
    list.replaceChildren();
    section.hidden = possible ? !Object.hasOwn(lookup, 'source_lookup_file') : !entries.length;
    if (possible) {
      $('possible-sources-note').textContent = section.hidden ? '' : t('u7.sourceNote', {file: lookup.source_lookup_file});
      $('possible-sources-empty').hidden = section.hidden || !!entries.length;
      $('possible-sources-empty').textContent = section.hidden ? '' : t('u7.noMatch', {file: lookup.source_lookup_file});
    }
    for (const source of entries) {
      const row = document.createElement('li'), heading = document.createElement('strong');
      row.dataset.sourceId = source.id; heading.textContent = source.id; row.append(heading);
      const metadata = document.createElement('p'); metadata.className = 'micro';
      metadata.textContent = [source.doc_type ? t('u5.type.' + source.doc_type) : '', source.journal, source.year,
        possible ? t('u7.chemblRecords', {n: source.n_records}) : t('u5.retainedRecords', {count: source.n_records})].filter(value => value != null && value !== '').join(' · ');
      row.append(metadata);
      if (possible && source.same_values?.length) {
        const same = document.createElement('p'); same.className = 'micro source-same-value';
        same.textContent = t('u7.sameValue', {values: source.same_values.join(t('u.idSeparator'))}); row.append(same);
      }
      const links = document.createElement('div'); links.className = 'source-links';
      for (const link of source.links || []) {
        // 仅产生用户主动点击的 http(s) 导航，不发起文献查询或预加载。
        let destination;
        try { destination = new URL(link.url); } catch (_) { destination = null; }
        if (!destination || !['http:', 'https:'].includes(destination.protocol)) continue;
        const anchor = document.createElement('a');
        const label = link.kind === 'patent' ? t('u5.patentLink', {id: link.label.replace(/^Patent /, '')}) : link.label;
        anchor.textContent = label; anchor.setAttribute('href', link.url);
        anchor.setAttribute('target', '_blank'); anchor.setAttribute('rel', 'noopener noreferrer');
        anchor.setAttribute('aria-label', label + ' — ' + link.url);
        links.append(anchor);
      }
      for (const reference of source.references || []) {
        if (reference.url) continue;
        const text = document.createElement('span'); text.className = 'source-reference';
        text.textContent = reference.value; links.append(text);
      }
      if (!links.querySelector('a')) {
        const unavailable = document.createElement('span'); unavailable.className = 'micro';
        unavailable.textContent = t('u5.noLinks'); links.append(unavailable);
      }
      row.append(links); list.append(row);
    }
  }
}
async function inspect(query) {
  // 语言切换可复用第二个参数中的检查结果，普通调用仍只传查询字符串。
  const cached = arguments[1] || null;
  var translateServer = globalThis.ChemGridMapI18n?.translateServer || (text => text);
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  const request = cached ? requestNumber : ++requestNumber;
  try {
    if (!job) return;
    const detail = cached || await getJSON('/api/inspect?' + query);
    if (request !== requestNumber) return;
    inspect.lastDetail = detail;
    if (!cached) status('');
    const row = detail.molecule, summary = detail.summary;
    const unmeasured = row.is_new || row.measurement_status === 'not_measured';
    currentCell = {row:row.grid_row,col:row.grid_col};
    if (unmeasured) {
      const references = new Map(points.filter(p => !p.is_new && p.measurement_status !== 'not_measured').map(p => [p.molecule_identity_key || p.molecule_id,p]));
      // 直接使用后端近邻表，包含第五名的全部并列分子。
      neighborRelations = new Map((detail.query_neighbors || []).map(p => references.get(p.measured_id)).filter(Boolean).map(p => [`${p.grid_row},${p.grid_col}`,'representation_only']));
    } else neighborRelations = new Map((detail.neighbors || []).map(p => [`${p.grid_row},${p.grid_col}`, p.neighbor_relation]));
    highlight();
    $('neighbor-key').hidden = unmeasured;
    const group = structureGroups.find(g => g.structure_group === row.structure_group);
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
    $('selected-group').hidden = !group;
    if (group) {
      $('selected-group-summary').textContent = t("d.{id} | {members} | core covers {covered}/{size}", {id: group.structure_group, members: countLabel(group.size, 'member'), covered: group.core_coverage_count, size: group.size});
      $('selected-group-button').onclick = () => compareGroup(group.structure_group);
      $('selected-group-highlight').onclick = () => selectGroup(highlightedGroup === group.structure_group ? '' : group.structure_group);
    }
    const wasClosed = $('inspector').hidden;
    $('inspector').hidden = false; $('map-layout').classList.add('has-selection');
    if (wasClosed) requestAnimationFrame(fit);
    $('selection').hidden = false; $('selection-hint').hidden = true;
    $('molecule-svg').innerHTML = detail.svg;
    globalThis.ChemGridMapI18n?.translateSVG?.($('molecule-svg'));
    $('selected-id').textContent = row.query_id || row.molecule_id || row.molecule_identity_key || row.canonical_smiles;
    const match = detail.query_match, matchNote = $('query-match-note');
    matchNote.replaceChildren();
    matchNote.hidden = !match || !['in_reference', 'duplicate_query'].includes(match.status);
    if (!matchNote.hidden) {
      const parts = [];
      if (match.status === 'duplicate_query') parts.push(t('u6.duplicateQuery', {query: match.query_id, first: match.first_query_id}));
      if (match.status === 'in_reference' || match.reference_id) {
        const title = $('selected-id').textContent, identity = row.parent_molecule_id || match.reference_id;
        const measured = identity && title !== identity ? t('d.{title} (parent {identity})', {title, identity}) : title;
        parts.push(t('u6.inReference', {query: match.status === 'duplicate_query' ? match.first_query_id : match.query_id, measured}));
      }
      const measurement = detail.new_molecule_measurement || match;
      if ((measurement.user_median ?? measurement.own_median) != null) parts.push(t('u8.matchMeasurement',{own:format(measurement.user_median ?? measurement.own_median),count:measurement.n_records ?? measurement.own_n_records,reference:format(measurement.reference_median),difference:format(measurement.median_difference,true)}));
      matchNote.textContent = parts.reduce((first, second) => t('d.{first} {second}', {first, second}));
    }
    $('cell-location').textContent = t("d.Grid row {row} · column {column} (zero-based)", {row: row.grid_row, column: row.grid_col});
    $('measurement-badge').hidden = !unmeasured; $('query-selection').hidden = !unmeasured;
    $('measured-selection').hidden = unmeasured;
    if (unmeasured) {
      showSources();
      $('evidence-summary').replaceChildren(); $('record-table').replaceChildren(); $('record-section').hidden = true;
      $('view-records').hidden = $('download-evidence').hidden = $('download-audit').hidden = true;
      showQuerySelection(detail); return;
    }
    $('neighborhood-svg').innerHTML = detail.neighborhood_svg;
    globalThis.ChemGridMapI18n?.translateSVG?.($('neighborhood-svg'));
    $('evidence-summary').replaceChildren();
    const stats = document.createElement('dl');
    const neighbors = summary.neighborhood;
    const fingerprintNeighbors = job.options.representation === 'morgan';
    const activityTerm = job.options.input_format === 'chembl' ? 'pChEMBL' : 'pActivity';
    const assumedSigma = summary.method === 'assumed_sigma' || summary.assumed_sigma != null;
    const genericRecords = summary.input_format === 'generic_activity';
    const fields = [[t("d.Median {term}", {term: activityTerm}),summary.median_pchembl],[t("d.{term} range", {term: activityTerm}),summary.minimum_pchembl == null ? null : `${format(summary.minimum_pchembl)}–${format(summary.maximum_pchembl)}`],[t("d.Source records"),summary.n_records],[t("d.Distinct assays"),summary.n_assays],[t("d.Distinct documents"),summary.n_documents]];
    if (genericRecords || !summary.median_verified) {
      const auditKey = summary.n_records === 1 ? "d.Median recomputed from {count} listed record" : "d.Median recomputed from {count} listed records";
      fields.splice(0,fields.length,[t("d.Activity class"),t('d.' + row.activity_class)],[t("d.Source-record audit"),genericRecords && summary.median_verified ? t(auditKey, {count: summary.n_records}) : t("d.Not available")]);
      if (job.options.input_format !== 'chembl') {
        const activity = row.activity_pchembl ?? row[job.options.value_col];
        if (job.options.compound_column) {
          fields.unshift([t("d.Median pActivity"), activity], [t("d.pActivity range"), row.pchembl_min == null ? null : `${format(row.pchembl_min)}–${format(row.pchembl_max)}`], [t("d.Source records"), row.n_records]);
        } else if (activity != null) fields.unshift([t("d.Activity (pActivity)"), activity]);
      }
    }
    if (row.class_confidence != null) {
      fields.push([t("d.Evidence units"), assumedSigma ? (summary.evidence_unit_source === 'one unit per compound' ? t("d.one per compound") : summary.evidence_unit_source) : t("d.{units} ({documents} documents)", {units: row.n_evidence_units, documents: row.n_documents})],
      [t("d.Document-balanced {term}", {term: activityTerm}), row.document_balanced_pchembl],
      [t("d.Class support"), `${percentage(row.class_confidence)} ${t('d.' + row.activity_class)} · ${t('d.' + row.class_status)}`]);
      if (summary.status_reason != null) {
        const threshold = summary.confidence_threshold == null ? (job.options.support_level ?? 0.8) : summary.confidence_threshold;
        const supportText = summary.status_reason === 'supported' ? t('d.supported') :
          summary.status_reason === 'low_support' ? t("d.low support (below {threshold})", {threshold}) :
          summary.status_reason === 'offset_pair' ? t("d.offset document pair") : summary.status_reason;
        fields.push([t("d.Support status"), supportText]);
      }
    }
    if (assumedSigma) fields.push([t("d.σ (assumed)"), format(summary.assumed_sigma ?? summary.sigma_doc)]);
    else if (summary.sigma_doc != null) fields.push(
      ['σ_doc', format(summary.sigma_doc)],
      [t("d.σ_doc degrees of freedom"), summary.sigma_doc_degrees_of_freedom == null ? t("d.not estimated (fallback)") : summary.sigma_doc_degrees_of_freedom],
      [t("d.Entries with multiple units"), summary.multi_unit_entry_count]);
    if (row.neighbor_recall != null) fields.push(
      [t("d.Neighbor recall (tie-inclusive)"), `${row.neighbor_shared_count}/${row.representation_neighbor_count} (${percentage(row.neighbor_recall)})`]);
    if (row.neighbor_recall_fixed_k != null) fields.push(
      [t("d.Neighbor recall (fixed k)"), percentage(row.neighbor_recall_fixed_k)]);
    if (summary.neighbor_kth_similarity != null) fields.push(
      [`${ordinal(neighbors.effective_k)} fingerprint neighbor`, `Tanimoto ${summary.neighbor_kth_similarity.toFixed(2)}`]);
    if (summary.neighbor_kth_distance != null) fields.push(
      [`${ordinal(neighbors.effective_k)} representation neighbor`, `distance ${summary.neighbor_kth_distance.toFixed(2)}`]);
    for (const [label,value] of fields) {
      const dt = document.createElement('dt'), dd = document.createElement('dd');
      dt.textContent = translateServer(label);
      // 证据单元来源可以是用户列名；只翻译代码生成的距离说明。
      dd.textContent = label.endsWith(' representation neighbor') ? translateServer(format(value)) : format(value);
      stats.append(dt,dd);
    }
    $('evidence-summary').append(stats);
    const marker = document.createElement('p'); marker.className = genericRecords ? 'micro' : !summary.quality_annotations_verified && summary.median_verified ? 'flag' : summary.wide_record_range ? 'note' : 'micro';
    marker.textContent = genericRecords ? summary.interpretation : !summary.quality_annotations_verified && summary.median_verified ? t("d.Quality flags not verified: original annotation rules are unavailable.") : summary.wide_record_range ? t("d.Records span more than one log unit or lie on both sides of the class thresholds. Class support already accounts for agreement among documents.") : summary.median_verified ? t("d.No wide record range rule triggered. This does not establish assay equivalence.") : summary.interpretation;
    marker.textContent = translateServer(marker.textContent);
    $('evidence-summary').append(marker);
    for (const text of summary.warnings || []) {
      if (!text.startsWith('Offset document pair:') && globalThis.ChemGridMapI18n?.getLanguage() !== 'zh') continue;
      const warning = document.createElement('p'); warning.className = 'flag';
      warning.dataset.documentOffsetWarning = String(text.startsWith('Offset document pair:')); warning.textContent = translateServer(text);
      $('evidence-summary').append(warning);
    }
    if (summary.median_verified && !genericRecords) { const verified = document.createElement('p'); verified.className = 'verified'; verified.textContent = summary.quality_annotations_verified ? t("d.Identity, counts, median, activity class and quality flags checked against retained records and saved rules.") : t("d.Record count and median verified; quality annotations remain unverified."); $('evidence-summary').append(verified); }
    showSources(detail.sources, detail);
    $('neighbor-title').textContent = `${fingerprintNeighbors ? 'Fingerprint' : 'Representation'} vs grid neighbors`;
    $('neighbor-title').textContent = translateServer($('neighbor-title').textContent);
    $('neighbor-summary').textContent = `${neighbors.representation} (${neighbors.representation_metric}), k=${neighbors.effective_k}: ${neighbors.shared_count} shared neighbors out of ${neighbors.representation_count} ${fingerprintNeighbors ? 'fingerprint' : 'representation'} neighbors and ${neighbors.grid_count} grid neighbors. Recall ${(100 * neighbors.recall).toFixed(1)}%.`;
    $('neighbor-summary').textContent = translateServer($('neighbor-summary').textContent);
    table(detail.neighbors, $('neighbor-table'));
    $('download-audit').hidden = !detail.audit_file;
    $('download-audit').onclick = () => download('/api/file?name=' + encodeURIComponent(detail.audit_file), detail.audit_file);
    $('download-evidence').hidden = !detail.total_records;
    $('view-records').hidden = !detail.total_records;
    $('record-section').hidden = !detail.total_records;
    const priority = ['source_row','molecule_identity_key','activity_pchembl_raw','assay_chembl_id','document_chembl_id','assay_description'];
    const ordered = genericRecords ? detail.records : detail.records.map(record => Object.fromEntries([...priority.filter(key => key in record),...Object.keys(record).filter(key => !priority.includes(key))].map(key => [key,record[key]])));
    table(ordered,$('record-table'),999,genericRecords ? detail.record_columns : null);
    $('record-count').textContent = t("d.{count} records{suffix}", {count: detail.total_records, suffix: !genericRecords && detail.total_records > 100 ? t('d. · first 100 shown; CSV contains all records') : ''});
  } catch (error) { if (request === requestNumber) fail(error); }
}
async function download(path,name) {
  try { const blob = await (await api(path)).blob(), url = URL.createObjectURL(blob); const a = document.createElement('a'); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url),10000); }
  catch (error) { fail(error); }
}
async function compareGroup(groupId, cached = null) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  try {
    const result = cached || await getJSON('/api/compare?group=' + encodeURIComponent(groupId));
    compareGroup.lastResult = result;
    const group = result.summary;
    $('comparison-title').textContent = `${group.structure_group}: ${countLabel(group.size, 'molecule')}`;
    $('core-svg').innerHTML = result.drawings.core_smiles || '';
    if (!result.drawings.core_smiles) $('core-svg').textContent = t('d.No ring scaffold');
    $('representative-svg').innerHTML = result.drawings.representative_smiles;
    const rule = t("d.Core rule: ") + coreRule(group) + (group.core_rule === 'coverage'
      ? t("d. The most frequent scaffold matches {covered}/{size}.", {covered: group.most_frequent_scaffold_coverage_count, size: group.size}) : '');
    $('core-coverage').textContent = t("d.{status}. This core matches {covered}/{members}; {exact} have exactly this scaffold. {distinct} distinct ring scaffolds; {acyclic}. {rule}", {status: t('d.' + group.core_status), covered: group.core_coverage_count, members: countLabel(group.size, 'member'), exact: group.core_exact_scaffold_count, distinct: group.distinct_scaffolds, acyclic: countLabel(group.acyclic_members, 'acyclic member'), rule});
    $('comparison-rule').textContent = t("d.{selected}/{members} shown. {rule} A scaffold summary is not a complete explanation of fingerprint similarity.", {selected: result.selected_count, members: countLabel(group.size, 'member'), rule: t('d.' + result.selection_rule)});
    $('comparison-svg').innerHTML = result.svg;
    const activityTerm = job.options.input_format === 'chembl' ? 'pChEMBL' : 'pActivity';
    const rows = result.members.map(m => ({ID:m.comparison_id, [activityTerm]:m.activity_pchembl, Class:m.activity_class, Core:group.core_member_ids.includes(m.comparison_id), Scaffold:m.scaffold_smiles, Row:m.grid_row, Column:m.grid_col}));
    table(rows, $('group-members'));
    $('group-members').querySelectorAll('tbody tr').forEach((tr, i) => {
      tr.tabIndex = 0; tr.setAttribute('role','button');
      const locate = () => { $('comparison-dialog').close(); inspect(`row=${rows[i].Row}&col=${rows[i].Column}`); };
      tr.onclick = locate; tr.onkeydown = e => { if (e.key === 'Enter') locate(); };
    });
    $('download-comparison').onclick = () => download('/api/file?name=' + encodeURIComponent(result.archive), result.archive);
    if (!$('comparison-dialog').open) $('comparison-dialog').showModal();
    translateSVGText($('comparison-svg'));
  } catch (error) { fail(error); }
}
$('compare-group').addEventListener('click', () => compareGroup($('group-select').value));
$('group-select').addEventListener('change', () => selectGroup($('group-select').value));
$('measurement-filter').addEventListener('change', () => { highlight(); renderLegend($('color-mode').value); });
$('color-mode').addEventListener('change', () => { colorMap(); highlight(); });
$('show-neighbor-links').addEventListener('change', highlight);
$('demo').addEventListener('click', () => loadInput(true));
$('csv').addEventListener('change', () => loadInput());
$('query-csv').addEventListener('change', loadQueryInput);
$('remove-query-input').addEventListener('click', removeQueryInput);
bindNewMolecules();
$('source-csv').addEventListener('change', loadSourceInput);
$('remove-source-input').addEventListener('click', removeSourceInput);
$('input-format').addEventListener('change',updateInputType);
$('unit-column').addEventListener('change', updateEvidenceOptions);
$('use-assumed-sigma').addEventListener('change', updateEvidenceOptions);
$('assumed-sigma').addEventListener('input', updateEvidenceOptions);
$('grouping').addEventListener('change', updateGrouping);
$('occupancy').addEventListener('input', densityPreview);
$('options').addEventListener('submit',build);
$('grid-tab').addEventListener('click', () => setView('grid'));
$('projection-tab').addEventListener('click', () => setView('projection'));
$('fit').addEventListener('click',fit);
$('zoom-in').addEventListener('click', () => zoom(1.4));
$('zoom-out').addEventListener('click', () => zoom(1/1.4));
$('download-all').addEventListener('click', () => download('/api/file?name=chemgridmap_results.zip','chemgridmap_results.zip'));
$('download-evidence').addEventListener('click', () => currentCell && download(`/api/evidence.csv?row=${currentCell.row}&col=${currentCell.col}`,'source_records.csv'));
$('find-molecule').addEventListener('submit',event => { event.preventDefault(); inspect('id=' + encodeURIComponent($('molecule-id').value.trim())); });
$('quit').addEventListener('click',async () => { try { const result = await getJSON('/api/shutdown',{method:'POST'}); status(result.message); document.querySelectorAll('button,input,select').forEach(el => el.disabled = true); } catch (error) { fail(error); } });
const viewport = $('map-viewport');
let drag = null;
viewport.addEventListener('wheel',event => { if (!job) return; event.preventDefault(); const rect = viewport.getBoundingClientRect(); zoom(Math.exp(-event.deltaY*.002),event.clientX-rect.left,event.clientY-rect.top); },{passive:false});
viewport.addEventListener('pointerdown',event => { if (!job || event.button !== 0) return; drag = {x:event.clientX,y:event.clientY,tx:transform.x,ty:transform.y,moved:false,target:event.target}; viewport.setPointerCapture(event.pointerId); });
viewport.addEventListener('pointermove',event => { if (!drag) return; const dx=event.clientX-drag.x,dy=event.clientY-drag.y; drag.moved ||= Math.abs(dx)+Math.abs(dy)>5; if (drag.moved) {transform.x=drag.tx+dx;transform.y=drag.ty+dy;applyTransform();} });
viewport.addEventListener('pointerup',() => { if (!drag) return; if (!drag.moved) { const cell=drag.target.closest('[data-grid-row]'); if (cell) inspect(`row=${cell.dataset.gridRow}&col=${cell.dataset.gridCol}`); } drag=null; });
viewport.addEventListener('pointercancel',() => drag=null);
viewport.addEventListener('keydown',event => { if (['+','=','-','0'].includes(event.key)) {event.preventDefault(); event.key==='0'?fit():zoom(event.key==='-'?1/1.4:1.4);} });
window.addEventListener('resize',fit);
for (const name of ['import', 'configure', 'explore']) $('step-' + name).addEventListener('click', () => { status(''); showPage(name); });
$('home').addEventListener('click', event => { event.preventDefault(); if (!busy) { status(''); showPage('import'); } });
$('replace-input').addEventListener('click', () => { status(''); showPage('import'); });
$('adjust-settings').addEventListener('click', () => { status(''); showPage('configure'); });
$('close-inspector').addEventListener('click', closeInspector);
for (const [button, dialog] of [['open-audit','audit-dialog'],['open-export','export-dialog'],['view-records','records-dialog']]) {
  $(button).addEventListener('click', () => $(dialog).showModal());
}
document.querySelectorAll('[data-close-dialog]').forEach(button => button.addEventListener('click', () => button.closest('dialog').close()));
$('options').addEventListener('invalid', event => {
  const details = event.target.closest('details'); if (details) details.open = true;
}, true);
const dropZone = $('drop-zone');
dropZone.addEventListener('dragover', event => { event.preventDefault(); if (!busy) dropZone.classList.add('drag-over'); });
dropZone.addEventListener('dragleave', () => dropZone.classList.remove('drag-over'));
dropZone.addEventListener('drop', event => {
  event.preventDefault(); dropZone.classList.remove('drag-over');
  if (busy) return;
  const files = event.dataTransfer.files;
  if (files.length !== 1) { fail(new Error("Choose one CSV file at a time.")); return; }
  loadInput(false, files[0]);
});
document.querySelector('.tabs').addEventListener('keydown', event => {
  if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
    event.preventDefault(); const next = view === 'grid' ? 'projection' : 'grid'; setView(next); $(next + '-tab').focus();
  }
});
densityPreview(); updateNavigation(); updateGrouping();
function restoreOptions(config) {
  var t = globalThis.ChemGridMapI18n?.t || ((key, values = {}) => key.slice(2).replace(/\{(\w+)\}/g, (_, name) => values[name] ?? `{${name}}`));
  $('grouping').value = config.grouping || 'butina';
  $('group-similarity').value = config.group_similarity ?? .6;
  $('cluster-col').value = config.cluster_col || ''; $('cluster-source').value = config.cluster_source || '';
  $('grid-refinement').value = config.grid_refinement || 'neighbors'; $('core-rule').value = config.core_rule || 'coverage';
  $('core-coverage-target').value = config.core_min_coverage ?? .8; $('fmcs-timeout').value = config.fmcs_timeout ?? 1;
  $('activity-confidence').checked = config.activity_confidence ?? true; $('evidence-model').value = config.evidence_model || 'gaussian';
  $('merge-identical-documents').checked = config.merge_identical_documents ?? true;
  $('support-level').value = config.support_level ?? .8;
  $('offset-display').value = config.offset_display || 'none';
  $('offset-min-shared').value = config.offset_min_shared ?? 3;
  $('offset-min-difference').value = config.offset_min_difference ?? 1;
  $('offset-sign-agreement').value = config.offset_sign_agreement ?? .8;
  $('identity-conflicts').value = config.identity_conflicts || 'error';
  const fields = {'input-format':'input_format',target:'target_id','activity-type':'activity_type',representation:'representation',projection:'projection',lower:'lower_threshold',upper:'upper_threshold',seed:'random_state',detail:'render_detail','conflict-range':'conflict_range_threshold',identity:'molecule_identity','structure-policy':'structure_policy','compound-column':'compound_column','smiles-col':'smiles_col','value-col':'value_col','document-column':'document_column','relation-column':'relation_column','label-col':'label_col'};
  fields['reference-column'] = 'reference_column';
  for (const [id,key] of Object.entries(fields)) $(id).value = config[key] ?? '';
  if (config.input_format === 'molecule' && config.compound_column) $('input-format').value = 'generic';
  $('value-unit').value = config.value_unit || 'p'; $('unit-column').value = config.unit_column || '';
  $('use-assumed-sigma').checked = config.assumed_sigma != null; $('assumed-sigma').value = config.assumed_sigma ?? .64;
  $('occupancy').value = config.grid_occupancy*100; densityPreview();
  $('raster').checked = config.output_formats.includes('png'); updateInputType();
  updateGrouping();
}
$('core-rule').addEventListener('change', updateGrouping);

function translateSVGText(container) {
  for (const node of container.querySelectorAll('text,title')) {
    if (!node.dataset.originalText) node.dataset.originalText = node.textContent;
    // Compound IDs and chemical drawing atom labels are not interface text.
    const original = node.dataset.originalText;
    const card = node.closest('[data-comparison-id]');
    if (card && original === card.dataset.comparisonId) continue;
    node.textContent = translateServer(original);
  }
}
function refreshValidation(input) {
  if (!input.setCustomValidity) return;
  input.setCustomValidity('');
  if (input.id === 'assumed-sigma' && !(Number(input.value) > 0)) {
    input.setCustomValidity(t('d.Assumed σ must be greater than zero.')); return;
  }
  if (window.ChemGridMapI18n.getLanguage() !== 'zh' || input.disabled) return;
  const validity = input.validity;
  if (validity.valueMissing) input.setCustomValidity(t(input.tagName === 'SELECT' ? 'd.Please select an item in the list.' : 'd.Please fill out this field.'));
  else if (validity.badInput) input.setCustomValidity(t('d.Please enter a number.'));
  else if (validity.rangeUnderflow) input.setCustomValidity(t('d.Value must be greater than or equal to {value}.', {value: input.min}));
  else if (validity.rangeOverflow) input.setCustomValidity(t('d.Value must be less than or equal to {value}.', {value: input.max}));
  else if (!validity.valid) input.setCustomValidity(t('d.Please enter a valid value.'));
}
async function refreshLanguage() {
  const focus = document.activeElement;
  const scroll = [window.scrollX, window.scrollY];
  updateInputType();
  $('build').textContent = t(busy ? 'd.Building map...' : 'd.Generate map →');
  $('selection-hint').textContent = t('d.Select a molecule on the map.');
  if (!inspect.lastDetail || $('inspector').hidden) $('neighbor-title').textContent = t('d.Representation vs grid neighbors');
  if (inputInfo) {
    $('input-summary').textContent = t('d.CSV loaded · stored only on this computer');
    $('input-stats').replaceChildren(stat(inputInfo.rows, t('d.Input rows')), stat(inputInfo.columns.length, t('d.Columns')));
    if (!inputInfo.preview.length) $('preview-table').textContent = t('d.No records available.');
  } else $('input-name').textContent = t('d.Your dataset');
  if (queryInputInfo) $('query-input-rows').textContent = t('d.{count} rows', {count: queryInputInfo.rows.toLocaleString()});
  if (sourceInputInfo) $('source-input-rows').textContent = t('d.{count} rows', {count: sourceInputInfo.rows.toLocaleString()});
  for (const select of document.querySelectorAll('#molecule-options select,#query-input-info select')) {
    for (const option of select.options) if (!option.value) option.textContent = t('d.None');
  }
  if (umapAvailable === false) $('projection').querySelector('[value=umap]').textContent = t('d.UMAP (not installed)');
  if (job && showResult.grouping) await showResult(true);
  if (currentCell && inspect.lastDetail) await inspect(null, inspect.lastDetail);
  if (highlightedGroup) {
    // 新组请求尚在运行时，更新标签但不借用上一个结构组的核心图。
    const cachedGroup = selectGroup.lastResult?.summary?.structure_group === highlightedGroup ? selectGroup.lastResult : {};
    await selectGroup(highlightedGroup, cachedGroup);
  }
  if (typeof newMolecules !== 'undefined' && newMolecules) renderNewList(true);
  if (typeof updateNewNavigation === 'function') updateNewNavigation();
  if ($('comparison-dialog').open && compareGroup.lastResult) await compareGroup(null, compareGroup.lastResult);
  applyTransform();
  for (const id of ['map-stage','molecule-svg','neighborhood-svg','comparison-svg']) translateSVGText($(id));
  if (status.last) status(status.last.message, status.last.kind);
  else $('build-progress').textContent = t('d.Preparing input records...');
  for (const input of document.querySelectorAll('#options input,#options select')) refreshValidation(input);
  window.ChemGridMapI18n.updateHelpLabels?.();
  if (focus?.isConnected && document.activeElement !== focus) focus.focus({preventScroll:true});
  window.scrollTo({left:scroll[0], top:scroll[1]});
}
window.ChemGridMapI18n.translateSVG = translateSVGText;
document.addEventListener('chemgridmap:languagechange', refreshLanguage);
$('options').addEventListener('input', event => refreshValidation(event.target));
$('options').addEventListener('change', event => refreshValidation(event.target));
$('options').addEventListener('invalid', event => refreshValidation(event.target), true);
refreshLanguage();

getJSON('/api/session').then(async info => {
  $('version').textContent = 'v' + info.version;
  umapAvailable = info.umap; demoAvailable = info.demo_available;
  if (!info.umap) { $('projection').querySelector('[value=umap]').disabled=true; $('projection').querySelector('[value=umap]').textContent=t("d.UMAP (not installed)"); }
  $('demo').disabled = !info.demo_available;
  if (info.input) { showInput(info.input); showPage('configure'); }
  newMolecules = info.new_molecules || null;
  if (info.query_input) showQueryInput(info.query_input, newMolecules?.options || {});
  if (info.source_input) showSourceInput(info.source_input);
  updateNavigation();
  if (info.job.state === 'running') { setBusy(true); showPage('explore'); poll(); }
  else if (info.job.state === 'done') {
    job = info.job;
    await showResult(); status('');
  }
  else if (info.job.state === 'error') status(info.job.message, 'error');
}).catch(fail);

}
initializeWorkspace().catch(error => {
  const status = document.getElementById('status');
  status.hidden = false; status.className = 'error'; status.textContent = error.message;
});
