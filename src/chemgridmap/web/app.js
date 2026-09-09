'use strict';
const $ = id => document.getElementById(id);
const token = location.hash.slice(1) || sessionStorage.getItem('chemgridmap-token');
if (token) sessionStorage.setItem('chemgridmap-token', token);
history.replaceState(null, '', location.pathname);
let busy = false, job = null, points = [], gridSVG = '', view = 'grid', currentCell = null, requestNumber = 0;
let transform = {s: 1, x: 0, y: 0}, drawingSize = {w: 1, h: 1};
let inputInfo = null, page = 'import', umapAvailable = false, demoAvailable = false;
const colors = {active: '#2a9d8f', medium: '#e6a23c', inactive: '#c44e52', unlabelled: '#b8c2cc'};

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
  $('build').disabled = busy || !inputInfo;
  $('demo').disabled = busy || !demoAvailable;
  $('projection').querySelector('[value=umap]').disabled = !umapAvailable;
  document.querySelectorAll('#find-molecule input,#find-molecule button,.tabs button,.zoom-tools button').forEach(el => el.disabled = !job || busy);
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
  currentCell = null; requestNumber++; highlight(); requestAnimationFrame(fit);
}

async function api(path, options = {}) {
  const response = await fetch(path, {...options, headers: {...options.headers, 'X-ChemGridMap-Token': token || ''}});
  if (!response.ok) {
    const data = await response.json();
    throw new Error(data.error || 'The request failed.');
  }
  return response;
}
const getJSON = async (path, options) => (await api(path, options)).json();
function status(message, kind = '') {
  $('status').textContent = message; $('status').className = kind;
  $('status').hidden = !message || kind === 'running';
  if (kind === 'running') $('build-progress').textContent = message;
}
function fail(error) { status(error.message, 'error'); }
function setBusy(value) {
  busy = value;
  document.querySelectorAll('#options input,#options select,#options button,#csv,#demo').forEach(el => el.disabled = value);
  $('quit').disabled = value;
  $('build').textContent = value ? 'Building map...' : 'Generate map →';
  updateNavigation();
}
function table(rows, container, columnLimit = 999) {
  container.replaceChildren();
  if (!rows.length) { container.textContent = 'No records available.'; return; }
  const columns = Object.keys(rows[0]).slice(0, columnLimit), el = document.createElement('table');
  const head = el.createTHead().insertRow();
  for (const col of columns) { const th = document.createElement('th'); th.textContent = col; head.append(th); }
  const body = el.createTBody();
  for (const row of rows) {
    const tr = body.insertRow();
    for (const col of columns) { const td = tr.insertCell(); td.textContent = row[col] ?? '—'; td.title = td.textContent; }
  }
  container.append(el);
}
function selectColumns(id, columns, preferred, optional = false) {
  $(id).replaceChildren();
  for (const name of (optional ? [''].concat(columns) : columns)) {
    const option = document.createElement('option'); option.value = name; option.textContent = name || 'None'; $(id).append(option);
  }
  $(id).value = columns.includes(preferred) ? preferred : (optional ? '' : columns[0]);
}
function resetResult() {
  job = null; currentCell = null; points = []; gridSVG = ''; requestNumber++;
  ['summary','audit','metrics','selection','record-section','download-all'].forEach(id => $(id).hidden = true);
  $('selection-hint').hidden = false; $('downloads').replaceChildren(); $('map-stage').replaceChildren();
  $('empty-map').hidden = false; $('saved-location').textContent = '';
  $('figure-downloads').replaceChildren(); $('legend').replaceChildren();
  $('run-caption').textContent = ''; closeInspector(); updateNavigation();
}
async function loadInput(demo = false, droppedFile = null) {
  try {
    if (busy) return;
    const file = droppedFile || $('csv').files[0];
    if (!demo && !file) return;
    setBusy(true); status('Reading the CSV...');
    const input = await getJSON(demo ? '/api/demo' : '/api/input', {
      method: 'POST', body: demo ? '' : file, headers: demo ? {} : {'X-Filename': encodeURIComponent(file.name)}
    });
    resetResult(); showInput(input, demo);
    status(input.targets.length > 1 ? 'Multiple targets detected. Select a target before building; they will not be pooled.' : '');
    setBusy(false);
    showPage('configure');
  } catch (error) { setBusy(false); $('build').disabled = !job && $('input-summary').hidden; fail(error); }
}
function showInput(input, demo = false) {
    inputInfo = input;
    $('input-summary').hidden = false;
    $('input-name').textContent = input.filename;
    $('input-summary').textContent = 'CSV loaded · stored only on this computer';
    $('input-stats').replaceChildren(stat(input.rows, 'Input rows'),stat(input.columns.length, 'Columns'));
    $('target').value = input.targets.length === 1 ? input.targets[0] : '';
    $('target-list').replaceChildren();
    input.targets.forEach(target => { const option = document.createElement('option'); option.value = target; $('target-list').append(option); });
    selectColumns('smiles-col', input.columns, input.detected.smiles || 'smiles');
    selectColumns('value-col', input.columns, input.detected.pchembl || 'activity_pchembl', true);
    selectColumns('label-col', input.columns, 'activity_class', true);
    $('input-format').value = demo || (input.detected.pchembl && (input.detected.standard_type || input.detected.target_id)) ? 'chembl' : 'molecule';
    updateInputType();
    table(input.preview, $('preview-table'), 12); $('preview').hidden = false;
    updateNavigation();
}
function updateInputType() {
  const chembl = $('input-format').value === 'chembl';
  $('chembl-options').hidden = !chembl; $('molecule-options').hidden = chembl;
  $('find-molecule').hidden = !chembl;
}
function options() {
  return {input_format: $('input-format').value, target_id: $('target').value.trim() || null,
    activity_type: $('activity-type').value, smiles_col: $('smiles-col').value,
    value_col: $('value-col').value || null, label_col: $('label-col').value || null,
    representation: $('representation').value, projection: $('projection').value,
    lower_threshold: Number($('lower').value), upper_threshold: Number($('upper').value),
    conflict_range_threshold: Number($('conflict-range').value), random_state: Number($('seed').value),
    grid_occupancy: Number($('occupancy').value) / 100, render_detail: $('detail').value,
    structure_policy: $('structure-policy').value, molecule_identity: $('identity').value,
    output_formats: $('raster').checked ? ['svg','png','pdf'] : ['svg']};
}
async function build(event) {
  event.preventDefault();
  if (busy) return;
  try {
    const config = options();
    await getJSON('/api/build', {method:'POST', body:JSON.stringify(config), headers:{'Content-Type':'application/json'}});
    resetResult(); setBusy(true); poll();
    showPage('explore'); status('Preparing input records...', 'running');
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
async function showResult() {
  restoreOptions(job.options);
  const report = job.report;
  $('summary').replaceChildren(stat(report.input_rows, 'Input records'),stat(job.molecules, 'Mapped molecules'),
    stat(report.retained_activity_records, 'Retained records'),stat(report.repeated_conflicted_molecules, 'Flagged molecules'));
  $('summary').hidden = false; $('audit').hidden = false; $('preview').open = false;
  $('audit-body').replaceChildren();
  const note = document.createElement('p');
  note.textContent = job.options.input_format === 'chembl' ? 'All retained records are summarized by median pChEMBL. A conflict flag is retained, not converted to Medium or silently removed. Aggregation does not harmonize assay conditions. An unflagged entry is not proof of experimental agreement.' : 'This map uses the supplied molecule-level annotations. It does not create an assay-record provenance trail.';
  $('audit-body').append(note);
  const exclusions = report.exclusion_counts || {};
  if (Object.keys(exclusions).length) {
    const reasons = document.createElement('p'); reasons.textContent = 'Exclusion reasons (may overlap): ' + Object.entries(exclusions).map(([k,v]) => `${k}: ${v}`).join('; '); $('audit-body').append(reasons);
  }
  const warnings = document.createElement('ul');
  for (const text of report.warnings || []) { const item = document.createElement('li'); item.textContent = text; warnings.append(item); }
  $('audit-body').append(warnings); $('audit').open = true;
  $('audit-count').textContent = (report.warnings || []).length;
  $('audit-count').hidden = !warnings.children.length;
  const svgFile = job.files.find(name => name === 'map.svg');
  gridSVG = await (await api('/api/file?name=' + encodeURIComponent(svgFile))).text();
  points = await getJSON('/api/points');
  updateNavigation(); showPage('explore');
  $('run-caption').textContent = [job.options.target_id || inputInfo.filename, job.options.representation === 'morgan' ? 'Morgan fingerprint' : 'Molecular descriptors', job.options.projection.toUpperCase(), `${Math.round(job.options.grid_occupancy * 100)}% occupancy`].join(' · ');
  view = 'grid'; setView('grid');
  $('download-all').hidden = false;
  $('downloads').replaceChildren();
  $('figure-downloads').replaceChildren();
  for (const name of job.files.filter(name => !['input.csv','chemgridmap_results.zip'].includes(name))) {
    const button = document.createElement('button');
    const figure = /^map\.(svg|png|pdf)$/.test(name);
    if (figure) {
      const extension = name.split('.').pop(); button.textContent = extension.toUpperCase();
      const caption = document.createElement('small'); caption.textContent = extension === 'png' ? 'Raster image' : 'Vector figure'; button.append(caption);
    } else button.textContent = name;
    button.addEventListener('click', () => download('/api/file?name=' + encodeURIComponent(name), name));
    $(figure ? 'figure-downloads' : 'downloads').append(button);
  }
  $('saved-location').textContent = 'Saved locally: ' + job.output_dir;
  table(Object.entries(job.metrics).map(([metric,value]) => ({Metric:metric,Value:value})), $('metrics-table')); $('metrics').hidden = false;
  const lower = job.options.lower_threshold, upper = job.options.upper_threshold;
  $('legend').replaceChildren();
  for (const [label,text] of [['inactive',`Inactive ≤${lower}`],['medium',`Medium (${lower}, ${upper})`],['active',`Active ≥${upper}`]]) {
    const item = document.createElement('span'), swatch = document.createElement('i'); swatch.className = label; item.append(swatch, document.createTextNode(text)); $('legend').append(item);
  }
  if (job.options.input_format === 'molecule') {
    $('legend').textContent = 'Colors follow supplied class labels, or the selected pChEMBL thresholds; gray = unlabelled.';
  } else { const flag = document.createElement('span'); flag.textContent = '◥ Conflict flag'; $('legend').append(flag); }
}
function projectionSVG() {
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg'); svg.setAttribute('viewBox','0 0 900 900'); svg.setAttribute('width','900'); svg.setAttribute('height','900');
  for (const point of points) {
    const circle = document.createElementNS(ns,'circle');
    circle.setAttribute('cx',40 + point.projection_x * 820); circle.setAttribute('cy',860 - point.projection_y * 820);
    circle.setAttribute('r','4'); circle.setAttribute('fill',colors[point.activity_class] || '#b8c2cc');
    circle.dataset.gridRow = point.grid_row; circle.dataset.gridCol = point.grid_col;
    const title = document.createElementNS(ns,'title'); title.textContent = `Grid cell (${point.grid_row}, ${point.grid_col}) · ${point.activity_class}`; circle.append(title); svg.append(circle);
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
  drawingSize = {w:box.width,h:box.height}; svg.setAttribute('width',box.width); svg.setAttribute('height',box.height);
  fit(); highlight();
}
function applyTransform() { $('map-stage').style.transform = `translate(${transform.x}px,${transform.y}px) scale(${transform.s})`; }
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
  $('map-stage').querySelectorAll('[data-grid-row]').forEach(el => {
    const selected = currentCell && Number(el.dataset.gridRow) === currentCell.row && Number(el.dataset.gridCol) === currentCell.col;
    el.classList.toggle('selected-cell',!!selected);
    if (el.tagName === 'circle') { el.setAttribute('stroke',selected ? '#183e92' : 'none'); el.setAttribute('stroke-width','3'); }
  });
}
function format(value) { return value == null ? '—' : typeof value === 'number' ? Number(value.toFixed(3)).toString() : String(value); }
async function inspect(query) {
  const request = ++requestNumber;
  try {
    if (!job) return;
    const detail = await getJSON('/api/inspect?' + query);
    if (request !== requestNumber) return;
    status('');
    const row = detail.molecule, summary = detail.summary;
    currentCell = {row:row.grid_row,col:row.grid_col}; highlight();
    const wasClosed = $('inspector').hidden;
    $('inspector').hidden = false; $('map-layout').classList.add('has-selection');
    if (wasClosed) requestAnimationFrame(fit);
    $('selection').hidden = false; $('selection-hint').hidden = true;
    $('molecule-svg').innerHTML = detail.svg; $('neighborhood-svg').innerHTML = detail.neighborhood_svg;
    $('selected-id').textContent = row.molecule_id || row.molecule_identity_key || row.canonical_smiles;
    $('cell-location').textContent = `Grid row ${row.grid_row} · column ${row.grid_col} (zero-based)`;
    $('evidence-summary').replaceChildren();
    const stats = document.createElement('dl');
    const fields = [['Median pChEMBL',summary.median_pchembl],['pChEMBL range',summary.minimum_pchembl == null ? null : `${format(summary.minimum_pchembl)}–${format(summary.maximum_pchembl)}`],['Source records',summary.n_records],['Distinct assays',summary.n_assays],['Distinct documents',summary.n_documents]];
    if (!summary.median_verified) fields.splice(0,fields.length,['Activity class',row.activity_class],['Source-record audit','Not available']);
    for (const [label,value] of fields) { const dt = document.createElement('dt'), dd = document.createElement('dd'); dt.textContent = label; dd.textContent = format(value); stats.append(dt,dd); }
    $('evidence-summary').append(stats);
    const marker = document.createElement('p'); marker.className = summary.is_conflicted ? 'flag' : 'micro';
    marker.textContent = summary.is_conflicted ? 'Conflict flagged. Review the underlying assays before interpretation.' : summary.median_verified ? 'No conflict rule triggered. This does not establish assay equivalence.' : summary.interpretation;
    $('evidence-summary').append(marker);
    if (summary.median_verified) { const verified = document.createElement('p'); verified.className = 'verified'; verified.textContent = 'Record count and median verified against source rows.'; $('evidence-summary').append(verified); }
    $('download-evidence').hidden = !detail.total_records;
    $('view-records').hidden = !detail.total_records;
    $('record-section').hidden = !detail.total_records;
    const priority = ['source_row','molecule_identity_key','activity_pchembl_raw','assay_chembl_id','document_chembl_id','assay_description'];
    const ordered = detail.records.map(record => Object.fromEntries([...priority.filter(key => key in record),...Object.keys(record).filter(key => !priority.includes(key))].map(key => [key,record[key]])));
    table(ordered,$('record-table'));
    $('record-count').textContent = `${detail.total_records} records${detail.total_records > 100 ? ' · first 100 shown; CSV contains all records' : ''}`;
  } catch (error) { if (request === requestNumber) fail(error); }
}
async function download(path,name) {
  try { const blob = await (await api(path)).blob(), url = URL.createObjectURL(blob); const a = document.createElement('a'); a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url),10000); }
  catch (error) { fail(error); }
}
$('demo').addEventListener('click', () => loadInput(true));
$('csv').addEventListener('change', () => loadInput());
$('input-format').addEventListener('change',updateInputType);
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
  if (files.length !== 1) { fail(new Error('Choose one CSV file at a time.')); return; }
  loadInput(false, files[0]);
});
document.querySelector('.tabs').addEventListener('keydown', event => {
  if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
    event.preventDefault(); const next = view === 'grid' ? 'projection' : 'grid'; setView(next); $(next + '-tab').focus();
  }
});
densityPreview(); updateNavigation();
function restoreOptions(config) {
  const fields = {'input-format':'input_format',target:'target_id','activity-type':'activity_type',representation:'representation',projection:'projection',lower:'lower_threshold',upper:'upper_threshold',seed:'random_state',detail:'render_detail','conflict-range':'conflict_range_threshold',identity:'molecule_identity','structure-policy':'structure_policy','smiles-col':'smiles_col','value-col':'value_col','label-col':'label_col'};
  for (const [id,key] of Object.entries(fields)) $(id).value = config[key] ?? '';
  $('occupancy').value = config.grid_occupancy*100; densityPreview();
  $('raster').checked = config.output_formats.includes('png'); updateInputType();
}
getJSON('/api/session').then(async info => {
  $('version').textContent = 'v' + info.version;
  umapAvailable = info.umap; demoAvailable = info.demo_available;
  if (!info.umap) { $('projection').querySelector('[value=umap]').disabled=true; $('projection').querySelector('[value=umap]').textContent='UMAP (not installed)'; }
  $('demo').disabled = !info.demo_available;
  if (info.input) { showInput(info.input); showPage('configure'); }
  updateNavigation();
  if (info.job.state === 'running') { setBusy(true); showPage('explore'); poll(); }
  else if (info.job.state === 'done') {
    job = info.job;
    await showResult(); status('');
  }
  else if (info.job.state === 'error') status(info.job.message, 'error');
}).catch(fail);
