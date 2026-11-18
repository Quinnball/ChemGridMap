// Fig. 7 capture: builds the thrombin IC50 reference map (records up to 2015) in the browser, inserts the later compounds and takes the screenshots.
// Adapted from the U10 real-browser check (u10_browser_check.mjs), which was adapted from the U9a script. Node 22 + Chromium, built-in modules only.
// Usage: CHROME=<Chromium binary> FIG7_INPUTS=<folder with the two thrombin CSV files> node paper/figures/source/fig7_capture.mjs
// Optional: WT=<ChemGridMap checkout> (default: this checkout), PYTHON (default: python3), CHECK_ROOT (default: paper/figures/fig7_work); SOFT=1 continues past failures.
import {spawn} from 'node:child_process';
import {mkdir, mkdtemp, writeFile, readFile, readdir, stat} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {join, relative, resolve, dirname} from 'node:path';
import {fileURLToPath} from 'node:url';
import {setTimeout as sleep} from 'node:timers/promises';
import assert from 'node:assert/strict';

const HERE = dirname(fileURLToPath(import.meta.url));
const ROOT = process.env.CHECK_ROOT || resolve(HERE, '../fig7_work');
const WT = process.env.WT || resolve(HERE, '../../..');
const SCREENS = `${ROOT}/fig7_screens`;
const PYTHON = process.env.PYTHON || 'python3';
const CHROME = process.env.CHROME;
if (!CHROME) throw new Error('Set CHROME to a Chromium or Chrome binary.');
const HEIGHT = 2000; const SCALE = 2;
const results = {node: process.version, started: new Date().toISOString(), status: 'running', checks: [], screenshots: [], observations: {}};
const processes = [];
let cdp, apiJSON;
const redact = value => String(value).replace(/(https?:\/\/127\.0\.0\.1:\d+\/?)#[^\s]+/g, '$1#[REDACTED_SECRET]');
const save = () => writeFile(`${ROOT}/fig7_results.json`, JSON.stringify(results, null, 2) + '\n');
const check = (name, condition, data = {}) => {
  if (process.env.SOFT && !condition) { (results.failures ||= []).push({name, data}); console.log('FAIL: ' + name + ' :: ' + JSON.stringify(data).slice(0, 600)); return; }
  assert.ok(condition, `${name}: ${JSON.stringify(data).slice(0, 1500)}`); results.checks.push({check: name, pass: true, ...data}); console.log('ok  ' + name);
};
const observe = (name, data) => { results.observations[name] = data; };

class CDP {
  constructor(ws) {
    this.ws = ws; this.next = 1; this.pending = new Map(); this.exceptions = [];
    ws.addEventListener('message', event => {
      const message = JSON.parse(event.data);
      if (message.id) { const r = this.pending.get(message.id); if (!r) return; clearTimeout(r.timer); this.pending.delete(message.id);
        message.error ? r.reject(new Error(JSON.stringify(message.error))) : r.resolve(message.result); }
      else if (message.method === 'Runtime.exceptionThrown') this.exceptions.push(message.params.exceptionDetails.exception?.description || message.params.exceptionDetails.text);
    });
  }
  static async open(url) { const ws = new WebSocket(url); await new Promise((res, rej) => { ws.addEventListener('open', res, {once: true}); ws.addEventListener('error', () => rej(new Error('ws failed')), {once: true}); }); return new CDP(ws); }
  call(method, params = {}, timeout = 45000) { const id = this.next++; return new Promise((resolve, reject) => { const timer = setTimeout(() => { this.pending.delete(id); reject(new Error(`CDP timeout: ${method}`)); }, timeout); this.pending.set(id, {resolve, reject, timer}); this.ws.send(JSON.stringify({id, method, params})); }); }
  async evaluate(expression, timeout = 20000) { const v = await this.call('Runtime.evaluate', {expression, awaitPromise: true, returnByValue: true}, timeout); if (v.exceptionDetails) throw new Error(v.exceptionDetails.exception?.description || v.exceptionDetails.text); return v.result.value; }
}
async function waitFor(expression, timeout = 30000) {
  const started = Date.now();
  while (Date.now() - started < timeout) {
    if (await cdp.evaluate(expression)) return Date.now() - started;
    const error = await cdp.evaluate(`document.getElementById('status')?.classList.contains('error') ? document.getElementById('status').textContent.slice(0,800) : ''`);
    if (error) throw new Error(`Browser application error: ${error}`);
    await sleep(200);
  }
  throw new Error(`Timed out waiting for: ${expression}`);
}
async function launch(name, command, args, env, pattern) {
  let output = '', exited = false;
  const child = spawn(command, args, {cwd: WT, env, stdio: ['ignore', 'pipe', 'pipe']}); processes.push(child);
  child.stdout.on('data', c => { output += c; }); child.stderr.on('data', c => { output += c; }); child.on('exit', () => { exited = true; });
  const started = Date.now();
  while (Date.now() - started < 30000 && !exited && !output.match(pattern)) await sleep(100);
  const match = output.match(pattern); if (!match) throw new Error(`${name} startup failed: ${redact(output).slice(-3000)}`); return match[0];
}
async function fetchJSON(url, options = {}) { const r = await fetch(url, {...options, signal: AbortSignal.timeout(60000)}); assert.equal(r.ok, true, `HTTP ${r.status} ${new URL(url).pathname}`); return r.json(); }
async function click(selector) {
  const box = await cdp.evaluate(`(() => { const el=document.querySelector(${JSON.stringify(selector)}); if (!el) throw new Error('Missing control: '+${JSON.stringify(selector)});
    el.scrollIntoView({block:'center'}); const r=el.getBoundingClientRect(); if (!r.width||!r.height) throw new Error('Hidden control: '+${JSON.stringify(selector)}); return {x:r.x+r.width/2,y:r.y+r.height/2}; })()`);
  for (const type of ['mouseMoved', 'mousePressed', 'mouseReleased']) await cdp.call('Input.dispatchMouseEvent', {type, button: 'left', clickCount: 1, ...box});
}
const select = (id, value) => cdp.evaluate(`(() => { const el=document.getElementById(${JSON.stringify(id)}); el.value=${JSON.stringify(value)}; if (el.value!==${JSON.stringify(value)}) throw new Error('Missing option ${value}'); el.dispatchEvent(new Event('change',{bubbles:true})); })()`);
const input = (id, value) => cdp.evaluate(`(() => { const el=document.getElementById(${JSON.stringify(id)}); el.value=${JSON.stringify(value)}; el.dispatchEvent(new Event('input',{bubbles:true})); el.dispatchEvent(new Event('change',{bubbles:true})); })()`);
const settle = () => cdp.evaluate('new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)))');
const parkPointer = async () => { await cdp.call('Input.dispatchMouseEvent', {type: 'mouseMoved', x: 8, y: 8}); await settle(); };
async function shot(name, focus = null) {
  if (focus) await cdp.evaluate(`document.querySelector(${JSON.stringify(focus)})?.scrollIntoView({block:'center'})`); else await cdp.evaluate('window.scrollTo({top:0})');
  await settle(); const {data} = await cdp.call('Page.captureScreenshot', {format: 'png'});
  const file = `${name}_1600x${HEIGHT}.png`; await writeFile(`${SCREENS}/${file}`, Buffer.from(data, 'base64')); results.screenshots.push(file);
}
const language = async lang => { await click(`[data-language="${lang}"]`); await waitFor(`document.documentElement.lang === ${JSON.stringify(lang === 'zh' ? 'zh-CN' : 'en')}`); await settle(); };
const findId = async id => {
  await input('molecule-id', id); await click('#find-molecule button');
  await waitFor(`!document.getElementById('inspector').hidden && document.getElementById('selected-id').textContent === ${JSON.stringify(id)} && !document.getElementById('sources').hidden`, 60000);
  await settle();
};
const highlightState = () => cdp.evaluate(`(() => { const cells=[...document.querySelectorAll('#map-stage [data-grid-row]')];
  return {cells: cells.length, source: cells.filter(e=>e.classList.contains('source-highlight')).length, group: cells.filter(e=>e.classList.contains('group-highlight')).length,
    faded: cells.filter(e=>e.classList.contains('group-faded')).length, status: $('source-highlight-status').hidden ? '' : $('source-highlight-label').textContent,
    dashed: cells.filter(e=>e.classList.contains('source-highlight')).filter(e => { const m = e.tagName==='circle'?e:e.querySelector(':scope > rect'); return m && getComputedStyle(m).strokeDasharray !== 'none'; }).length}; })()`);
const sourceRows = () => cdp.evaluate(`[...document.querySelectorAll('#source-list > li')].map(li => ({id: li.dataset.sourceId, button: li.querySelector('button[id^="highlight-source-"]')?.id || '', text: li.querySelector('button[id^="highlight-source-"]')?.textContent || li.textContent}))`);
async function hashTree(dir) {
  const out = {};
  async function walk(d) { for (const e of await readdir(d, {withFileTypes: true})) { const p = join(d, e.name); if (e.isDirectory()) await walk(p); else out[relative(dir, p)] = createHash('sha256').update(await readFile(p)).digest('hex'); } }
  await walk(dir); return out;
}
const cssEsc = s => s.replace(/"/g, '\\"');



async function upload(id, path) {
  const {root} = await cdp.call('DOM.getDocument', {depth:0});
  const {nodeId} = await cdp.call('DOM.querySelector', {nodeId:root.nodeId, selector:`#${id}`});
  await cdp.call('DOM.setFileInputFiles', {nodeId, files:[path]});
}
async function shotClip(name, selector, pad = 8) {
  await cdp.evaluate(`document.querySelector(${JSON.stringify(selector)}).scrollIntoView({block:'start'})`); await settle();
  const r = await cdp.evaluate(`(() => { const r = document.querySelector(${JSON.stringify(selector)}).getBoundingClientRect(); return {x: r.x + scrollX, y: r.y + scrollY, width: r.width, height: r.height}; })()`);
  const {data} = await cdp.call('Page.captureScreenshot', {format: 'png', clip: {x: Math.max(0, r.x - pad), y: Math.max(0, r.y - pad), width: r.width + 2 * pad, height: r.height + 2 * pad, scale: 1}});
  await writeFile(`${SCREENS}/${name}.png`, Buffer.from(data, 'base64')); results.screenshots.push(name);
}
const IN = process.env.FIG7_INPUTS;
if (!IN) throw new Error('Set FIG7_INPUTS to the folder with the two thrombin CSV files.');
try {
  await mkdir(SCREENS, {recursive: true});
  const runRoot = await mkdtemp(`${ROOT}/fig7-run-`), profile = await mkdtemp(`${ROOT}/chrome-profile-`);
  const appUrl = await launch('app', PYTHON, ['-m', 'chemgridmap.app', '--no-browser', '--port', '0', '--output-dir', `${runRoot}/runs`],
    {...process.env, PYTHONPATH: 'src', PYTHONUNBUFFERED: '1', PYTHONDONTWRITEBYTECODE: '1', MPLCONFIGDIR: `${ROOT}/matplotlib`, NUMBA_CACHE_DIR: `${ROOT}/numba_cache`}, /http:\/\/127\.0\.0\.1:\d+\/#[^\s]+/);
  const app = new URL(appUrl), secret = app.hash.slice(1);
  apiJSON = (path, options = {}) => fetchJSON(`${app.origin}${path}`, {...options, headers: {'X-ChemGridMap-Token': secret, ...(options.headers || {})}});
  const websocket = await launch('chrome', CHROME, ['--headless=new', '--no-sandbox', '--remote-debugging-port=0', `--user-data-dir=${profile}`, '--no-first-run', '--no-default-browser-check', '--disable-background-networking', '--disable-extensions', '--lang=en-US', `--window-size=1600,${HEIGHT}`, 'about:blank'], process.env, /ws:\/\/127\.0\.0\.1:\d+\/devtools\/browser\/[^\s]+/);
  const target = await fetchJSON(`http://${new URL(websocket).host}/json/new?about:blank`, {method: 'PUT'});
  cdp = await CDP.open(target.webSocketDebuggerUrl); await cdp.call('Page.enable'); await cdp.call('Runtime.enable'); await cdp.call('DOM.enable');
  await cdp.call('Page.addScriptToEvaluateOnNewDocument', {source: "window.$ = id => document.getElementById(id);"});
  await cdp.call('Emulation.setDeviceMetricsOverride', {width: 1600, height: HEIGHT, deviceScaleFactor: SCALE, mobile: false});
  await cdp.call('Page.navigate', {url: appUrl}); await waitFor(`document.readyState==='complete' && !document.getElementById('demo').disabled`);
  await language('en');
  await upload('csv', `${IN}/thrombin_chembl37_ic50_up_to_2015.csv`);
  await waitFor(`!document.getElementById('page-configure').hidden && !document.getElementById('build').disabled`, 60000);
  await click('#build'); await waitFor(`!!document.querySelector('#map-stage [data-grid-row]') && !document.getElementById('summary').hidden && !document.getElementById('build').disabled`, 1500000);
  await settle(); await parkPointer();
  const refCells = await cdp.evaluate(`document.querySelectorAll('#map-stage [data-grid-row]').length`);
  check('thrombin reference map (records up to 2015) rendered', refCells === 1429, {refCells});
  await click('#open-new-molecules'); await waitFor(`$('new-molecules-dialog').open`);
  await upload('query-csv', `${IN}/thrombin_compounds_first_reported_after_2015.csv`);
  await waitFor(`!$('query-input-info').hidden`);
  for (const [id, value] of Object.entries({'query-smiles-col': 'canonical_smiles', 'query-id-col': 'query_id', 'new-activity-col': 'activity_pchembl_raw'})) await select(id, value);
  await select('new-unit', 'p'); await select('new-document-col', 'document_chembl_id');
  await waitFor(`!$('insert-new-molecules').disabled`);
  await shotClip('fig7_dialog', '#new-molecules-dialog');
  const t0 = Date.now(); await click('#insert-new-molecules');
  await waitFor(`!$('new-molecules-dialog').open && !$('open-new-list').hidden && !!document.querySelector('#map-stage [data-new-molecule]')`, 300000);
  observe('insert_ms', Date.now() - t0); await settle(); await parkPointer();
  const state = await cdp.evaluate(`({mode: $('color-mode').value, newCells: document.querySelectorAll('#map-stage [data-new-molecule]').length, status: $('status')?.textContent})`);
  observe('after_insert', state);
  check('464 later compounds inserted; layer switched to neighbor evidence', state.newCells === 464 && state.mode === 'neighbor_evidence', state);
  await click('#fit'); await settle(); await sleep(500); await settle(); await parkPointer();
  observe('rects', await cdp.evaluate(`Object.fromEntries(['map-viewport','map-stage','inspector','selection'].map(id => { const e = $(id); if (!e) return [id, null]; const r = e.getBoundingClientRect(); return [id, {x: r.x, y: r.y, w: r.width, h: r.height}]; }))`));
  await shot('fig7_map_neighbor_evidence');
  await shotClip('fig7_legend', '#legend', 4);
  await shotClip('fig7_map_viewport', '#map-viewport', 2);
  const points = await apiJSON('/api/points');
  const newPoints = points.filter(p => p.is_new);
  observe('new_points', newPoints.length);
  observe('placement_distance', (() => { const d = newPoints.map(p => Number(p.placement_distance)).filter(Number.isFinite).sort((a, b) => a - b);
    return d.length ? {n: d.length, median: d[Math.floor((d.length - 1) / 2)] / 2 + d[Math.ceil((d.length - 1) / 2)] / 2, max: d[d.length - 1]} : null; })());
  const pick = newPoints.filter(p => p.neighbor_evidence === 'active' && p.activity_class === 'active' && p.tanimoto < 0.95 && p.tanimoto >= 0.7).sort((a, b) => (b.tanimoto ?? 0) - (a.tanimoto ?? 0) || String(a.molecule_id).localeCompare(String(b.molecule_id)))[0];
  const pickId = String(pick?.molecule_id ?? pick?.query_id);
  observe('picked', pick ? {id: pickId, row: pick.grid_row, col: pick.grid_col, tanimoto: pick.tanimoto, measured: pick.activity_class} : null);
  await click('#open-new-list'); await waitFor(`$('new-list-dialog').open`); await settle(); await sleep(300);
  await shotClip('fig7_new_list', '#new-list-dialog', 4);
  const found = await cdp.evaluate(`(() => { const rows = Array.from(document.querySelectorAll('#new-list-table tbody tr')); const row = rows.find(tr => Array.from(tr.cells).some(td => td.textContent.trim() === ${JSON.stringify(pickId)})) || rows[0]; row.click(); return row.cells[0].textContent; })()`);
  observe('clicked_row', found);
  await waitFor(`!$('new-list-dialog').open && !$('inspector').hidden`, 60000); await settle(); await sleep(500); await parkPointer();
  observe('selected', await cdp.evaluate(`$('selected-id').textContent`));
  observe('rects_selected', await cdp.evaluate(`Object.fromEntries(['map-viewport','inspector','selection'].map(id => { const r = $(id).getBoundingClientRect(); return [id, {x: r.x, y: r.y, w: r.width, h: r.height, sh: $(id).scrollHeight}]; }))`));
  await shot('fig7_selected_full');
  await shotClip('fig7_inspector', '#inspector', 2);
  await shotClip('fig7_query_summary', '#query-summary', 4);
  await shotClip('fig7_query_neighbors', '#query-neighbor-table', 4);
  observe('inspector_text', await cdp.evaluate(`$('inspector').innerText.slice(0, 4000)`));
  observe('legend_text', await cdp.evaluate(`$('legend').innerText`));
  observe('chromium', await cdp.evaluate('navigator.userAgent'));
  await shotClip('fig7_map_selected', '#map-viewport', 2);
  check('no JavaScript exceptions on the page', cdp.exceptions.length === 0, {exceptions: cdp.exceptions});
  results.status = 'passed';
} catch (error) {
  results.status = 'failed'; results.error = redact(error.stack || error);
  try { const {data} = await cdp.call('Page.captureScreenshot', {format: 'png'}); await writeFile(`${ROOT}/fig7_failure_v07.png`, Buffer.from(data, 'base64')); } catch {}
  console.error(results.error); process.exitCode = 1;
} finally {
  results.exceptions = cdp?.exceptions || []; await save();
  console.log(JSON.stringify({status: results.status, checks: results.checks.length, screenshots: results.screenshots, obs: results.observations}).slice(0, 2000));
  cdp?.ws.close(); for (const child of processes) if (child.exitCode === null) child.kill('SIGTERM');
}
