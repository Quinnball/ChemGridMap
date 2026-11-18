'use strict';
// 执行生产函数；DOM 只覆盖本契约访问的表面，不替代浏览器排版和原生按键验收。
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const {runtime, realFunction} = require('./u4_web_contract.cjs');
const assets = process.argv[2];
const source = fs.readFileSync(path.join(assets, 'app.js'), 'utf8');
const sample = {run_name:'saved-1', created_at:'2026-11-12T12:32:36.551872+00:00', filename:'CAII.csv',
  reference_map:{label:'CHEMBL205 · Carbonic anhydrase 2 · IC50'}, molecules:1145,
  core_rule:'coverage_brics', new_molecules:294};
const point = {query_id:'M001', molecule_id:'M001', is_new:true, grid_row:0, grid_col:1,
  measurement_status:'measured', activity_class:'inactive', measured_class:'inactive',
  neighbor_evidence:'active', neighbor_comparison:'disagrees'};

async function setup() {
  const env = await runtime({saved:'zh'}), {context, node, document} = env;
  const $ = id => { const value = document.getElementById(id); assert.ok(value, id); return value; };
  const requests = [], grouping = {groups:[], protocol:{method:'none'}};
  const job = {state:'done', files:['map.svg'], molecules:2, metrics:{}, output_dir:'/local/saved',
    options:{representation:'morgan', projection:'pca', grid_occupancy:.7, core_rule:'coverage_brics'},
    report:{input_rows:2, retained_activity_records:1, repeated_conflicted_molecules:0, warnings:[]}};
  Object.assign(context, {$, t:env.api.t, translateServer:env.api.translateServer, job, busy:false,
    points:[{...point}], mapCells:[], newMolecules:{rows:[point], files:['map_with_new_molecules.svg'], options:{}},
    view:'grid', gridSVG:'<svg/>', transform:{s:1, x:0, y:0}, tileWidth:100, hasStructures:true,
    inputInfo:null, queryInputInfo:null, sourceInputInfo:null, currentCell:null, highlightedGroup:'', structureGroups:[],
    umapAvailable:true, demoAvailable:true, colors:{active:'#2a9d8f', inactive:'#c44e52'}, fadedColors:{},
    updateNewNavigation(){}, updateEvidenceOptions(){}, updateInputType(){}, refreshValidation(){},
    renderSourceHighlight(){}, renderCoreRuleStatus(){}, renderNewList(){}, renderLegend(){},
    renderCoreRuleOptions(){}, renderGroupOptions(){}, renderReferenceMap(){}, showDataAudit(){},
    showInput(){}, showQueryInput(){}, showSourceInput(){}, restoreOptions(){}, showPage(){},
    selectGroup(){}, compareGroup(){}, clearSourceHighlight(){}, highlight(){},
    referenceMapDetails:() => ({description:sample.reference_map.label}),
    stat(value, label){const item=node('span');item.textContent=`${value} ${label}`;return item;}, table(){},
    inspect:async query => requests.push({inspect:query}),
    resetResult(){ $('map-stage').replaceChildren(); context.mapCells=[]; },
    fit(){ context.applyTransform(); },
    api:async () => ({text:async () => '<svg/>'}),
    getJSON:async (url, options) => {
      requests.push({url, options});
      if(url === '/api/new-molecules') return context.newMolecules;
      if(url === '/api/points') return context.points;
      if(url === '/api/groups') return grouping;
      if(url === '/api/previous-runs') return {runs:[sample]};
      if(url === '/api/open-run') return {job, new_molecules:context.newMolecules, restored:{new_set:'manual', total:1, recorded:1}};
      if(url === '/api/new-molecules/result') return {job, point:context.points[0]};
      throw Error(`Unexpected request ${url}`);
    }});
  vm.runInContext(['refreshLanguage','applyTransform','zoom','colorMap','setView','showResult',
    'translateSVGText','updateNavigation','setBusy','status','fail'].map(realFunction).join('\n'), context);
  context.showResult.grouping = grouping;
  // 原生 DOM 在浏览器中解析 SVG；这里用等价节点替身覆盖解析这一外部边界。
  function mapFixture() {
    const svg=node('svg'), reference=node('g'), cell=node('g'), rect=node('rect');
    svg.viewBox={baseVal:{x:0,y:0,width:200,height:100}};
    reference.setAttribute('data-reference-map','true');
    cell.setAttribute('data-grid-row','0');cell.setAttribute('data-grid-col','1');cell.setAttribute('data-new-molecule','true');
    rect.setAttribute('x','100');rect.setAttribute('y','0');rect.setAttribute('width','100');rect.setAttribute('fill','#c44e52');
    const originalQuery=cell.querySelector.bind(cell);
    cell.querySelector=selector => selector.startsWith(':scope > ')
      ? cell.children.find(child => child.matches(selector.slice(9))) || null : originalQuery(selector);
    cell.replaceWith=function(fresh){const parent=this.parentNode,index=parent.childNodes.indexOf(this);parent.childNodes[index]=fresh;fresh.parentNode=parent;this.parentNode=null;};
    cell.append(rect);svg.append(reference,cell);
    return {svg,cell};
  }
  function mount() {const fixture=mapFixture();$('map-stage').replaceChildren(fixture.svg);context.mapCells=[fixture.cell];return fixture;}
  Object.defineProperty($('map-stage'),'innerHTML',{get:() => '<svg/>',set:mount});
  context.DOMParser=class {parseFromString(){return mapFixture().svg;}};
  document.importNode=value=>value;
  mount();$('color-mode').value='activity';context.colorMap();
  // 使用 app.js 的实际事件绑定；包装仅收集异步任务，便于等待断言。
  const refreshLanguage=context.refreshLanguage;let languageWork;
  context.refreshLanguage=()=>languageWork=refreshLanguage();
  const binding="document.addEventListener('chemgridmap:languagechange', refreshLanguage);";
  assert.ok(source.includes(binding));vm.runInContext(binding,context);
  return {...env, $, requests, job, grouping, mount, languageWork:()=>languageWork};
}

function assertMarkers(env, measured='inactive', neighbor='active') {
  const expected=env.api.t('u9c.mismatchTooltip',{measured:env.api.t('d.'+measured),neighbor:env.api.t('d.'+neighbor)});
  const original=env.context.mapCells[0].querySelector(':scope > [data-neighbor-mismatch]');
  const overlays=env.$('map-stage').querySelectorAll('[data-mismatch-overlay]');
  assert.ok(original);assert.equal(overlays.length,1);assert.equal(overlays[0].children.length,1);
  for(const marker of [original,overlays[0].children[0]]) {
    assert.equal(marker.querySelector('title').textContent,expected);
    assert.equal(marker.getAttribute('aria-label'),expected);
  }
}

const scenarios = {
  async language() {
    const env=await setup(),{context,api}=env;assertMarkers(env);
    const originalShow=context.showResult;
    for(const language of ['en','zh']) {
      let complete;context.showResult=()=>new Promise(resolve=>{complete=resolve;});context.showResult.grouping=env.grouping;
      api.setLanguage(language);
      assert.equal(typeof complete,'function','语言事件已经进入异步摘要刷新');
      assertMarkers(env); // 摘要 await 尚未结束，菱形的 title 和 aria-label 必须已经更新。
      complete();await env.languageWork();assertMarkers(env);
    }
    context.showResult=originalShow;
  },
  async redraw() {
    const env=await setup(),{context,api,$,requests}=env;
    for(const language of ['en','zh']) {
      api.setLanguage(language);await env.languageWork();
      for(const scale of [.02,.3,1,5]) {context.transform.s=scale;context.applyTransform();assertMarkers(env);}
      for(const mode of ['activity','support','neighbor_evidence','structure']) {$('color-mode').value=mode;context.colorMap();assertMarkers(env);}
      // 实际保存入口 -> 刷新接口 -> 单格替换 -> colorMap/applyTransform，覆盖首次记录和修改。
      for(const [measured,neighbor] of [['inactive','active'],['active','inactive']]) {
        context.points=[{...point,measured_class:measured,activity_class:measured,neighbor_evidence:neighbor}];
        context.renderMeasuredResult.detail={molecule:{...point}};
        $('result-mode').value='class';$('result-class').value=measured;
        await context.saveMeasuredResult({preventDefault(){}});
        assert.equal($('result-feedback').classList.contains('error'),false,$('result-feedback').textContent);
        const saved=requests.findLast(request=>request.url==='/api/new-molecules/result');assert.ok(saved);
        assert.equal(JSON.parse(saved.options.body).measurement.measured_class,measured);
        assertMarkers(env,measured,neighbor);
      }
      // 实际打开入口 -> 快照恢复 -> showResult -> setView -> colorMap。
      $('previous-run-select').value=sample.run_name;
      await context.openSelectedRun();assert.equal($('status').classList.contains('error'),false,$('status').textContent);
      assert.equal(JSON.parse(requests.findLast(request=>request.url==='/api/open-run').options.body).run_name,sample.run_name);
      assertMarkers(env,'active','inactive');
      context.points=[{...point}];context.colorMap();
    }
  },
  async time() {
    const {context,api}=await setup();
    assert.equal(new Date(sample.created_at).getTimezoneOffset(),-480,'测试子进程必须在 UTC+8');
    for(const language of ['en','zh']) {
      api.setLanguage(language);
      assert.equal(context.previousRunTime(sample.created_at),'2026-11-12 20:32');
      assert.equal(context.previousRunTime('2026-01-02T16:05:59Z'),'2026-01-03 00:05');
      assert.match(context.previousRunLabel(sample),/2026-11-12 20:32/);
      assert.doesNotMatch(context.previousRunLabel(sample),/T12:32|\+00:00/);
    }
  },
  async list() {
    const env=await setup(),{context,api,$,document}=env;
    const runs=Array.from({length:25},(_,index)=>({...sample,run_name:'saved-'+index}));
    for(const language of ['en','zh']) {
      api.setLanguage(language);await env.languageWork();context.renderPreviousRuns({runs});
      const list=$('previous-run-select'),rows=list.querySelectorAll('label'),radios=list.querySelectorAll('input');
      assert.equal(list.tagName,'FIELDSET');assert.equal(list.getAttribute('aria-label'),api.t('u9c.previousList'));
      assert.equal(rows.length,20);assert.equal(radios.length,20);
      for(const row of rows) {
        const lines=row.querySelector('.previous-run-content').children;
        assert.deepEqual(lines.map(line=>line.textContent),['2026-11-12 20:32',sample.reference_map.label,
          api.t('u9c.previousCounts',{molecules:1145,newCount:294}),context.coreRuleName('coverage_brics')]);
        assert.equal(row.title,context.previousRunLabel(sample));assert.ok(row.title.includes('CAII.csv'));
      }
      for(const radio of radios) {
        assert.equal(radio.type,'radio');assert.equal(radio.name,'previous-run');
        assert.ok(radio.closest('label').textContent.trim(),'隐式 label 提供无障碍名称');
        assert.notEqual(radio.getAttribute('tabindex'),'-1');assert.equal(radio.events.has('keydown'),false,'保留原生方向键语义');
      }
      assert.equal(radios.filter(radio=>radio.checked).length,1);
      radios.forEach((radio,index)=>{radio.checked=index===3;});radios[3].dispatchEvent({type:'change',bubbles:true});
      assert.equal(list.value,'saved-3');assert.equal($('open-selected-run').disabled,false);
      context.renderPreviousRuns({runs});assert.equal(list.value,'saved-3');
      assert.equal(list.querySelectorAll('input').find(radio=>radio.checked).value,'saved-3');
      context.setBusy(true);assert.equal(list.disabled,true);assert.equal($('open-selected-run').disabled,true);
      context.setBusy(false);assert.equal(list.disabled,false);
      context.renderPreviousRuns({runs:[]});assert.equal(list.value,'');assert.equal($('open-selected-run').disabled,true);
      assert.equal($('previous-results-status').textContent,api.t('u9c.previousEmpty'));
      assert.equal(context.previousRunReference({reference_map:{description:'description'}}),'description');
      assert.equal(context.previousRunReference({}),'—');
    }
    assert.equal(document.querySelectorAll('[data-help]').length,81);
  },
  async layout() {
    const css=fs.readFileSync(path.join(assets,'app.css'),'utf8');
    const rule=selector=>{const block=css.split('}').find(item=>item.slice(item.lastIndexOf('{')+1)&&item.slice(0,item.lastIndexOf('{')).trim().endsWith(selector));assert.ok(block,selector);return block.slice(block.lastIndexOf('{')+1);};
    const list=rule('#previous-run-select'),row=rule('#previous-run-select .previous-run-option');
    const height=value=>Number(value.match(/(?:^|;)\s*height:\s*(\d+)px/)?.[1]);
    assert.match(css,/\*\s*\{[^}]*box-sizing:\s*border-box/);
    assert.match(list,/padding:\s*0/);assert.match(list,/overflow-y:\s*auto/);
    assert.match(row,/margin:\s*0/);assert.match(row,/border-bottom:\s*1px/);
    assert.match(rule('.previous-run-content > *'),/white-space:\s*nowrap/);
    assert.equal((css.match(/#previous-run-select\s*\{/g)||[]).length,1,'宽度媒体查询没有覆盖列表固定高度');
    for(const width of [1280,1600]) {
      assert.ok((height(list)-2)/height(row)>=5,`${width}px 窗口下固定 CSS 尺寸至少容纳五条`);
    }
    assert.match(css,/#previous-run-select \.previous-run-option:focus-within\s*\{[^}]*outline:/);
  }
};
(async()=>{const name=process.argv[3];assert.ok(name in scenarios);await scenarios[name]();console.log(`U9c-R1 ${name}: passed`);})()
  .catch(error=>{console.error(error);process.exitCode=1;});
