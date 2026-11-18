'use strict';
const assert=require('node:assert/strict'),vm=require('node:vm'),fs=require('node:fs'),path=require('node:path');
const {runtime,realFunction}=require('./u4_web_contract.cjs');
const sample={run_name:'20261112_190300',created_at:'2026-11-12T19:03:00+08:00',filename:'CAII.csv',reference_map:{label:'CHEMBL205 · Carbonic anhydrase 2 · IC50'},molecules:1145,core_rule:'coverage_brics',new_molecules:294};
async function setup(){
 const env=await runtime(),$=id=>env.document.getElementById(id),calls=[];
 Object.assign(env.context,{$,t:env.api.t,translateServer:env.api.translateServer,busy:false,mapCells:[],points:[],transform:{s:1,x:0,y:0},job:null,newMolecules:null,
  updateNavigation(){ $('open-selected-run').disabled=env.context.busy || !$('previous-run-select').value; },setBusy(value){env.context.busy=value;env.context.updateNavigation();},
  resetResult(){calls.push('reset');},showInput(value){calls.push(['input',value]);},showQueryInput(...value){calls.push(['query',...value]);},showSourceInput(value){calls.push(['source',value]);},showResult:async()=>calls.push('show'),
  getJSON:async(url,options)=>{calls.push([url,options]);return {runs:[sample],unavailable:{missing_snapshot:2,incompatible:1,modified:3}};},
  inspect:async(value)=>calls.push(['inspect',value])});
 vm.runInContext(['previousRunLabel','renderPreviousRuns','loadPreviousRuns','restoredSessionMessage','applySessionSnapshot','openSelectedRun','status','fail'].map(realFunction).join('\n'),env.context);
 return {...env,$,calls};
}
const scenarios={
 async progress(){
  const {context,$,api}=await setup();let callback,delay;
  context.setTimeout=(fn,ms)=>{callback=fn;delay=ms;return 42;};
  for(const language of ['en','zh']){
   api.setLanguage(language);context.switchCoreRule.busy=true;context.switchCoreRule.request=3;context.switchCoreRule.message={rule:'coverage_brics',n:100,done:0,total:100};
   context.getJSON=async(url)=>{assert.equal(url,'/api/core-progress');return {running:true,rule:'coverage_brics',done:27,total:100};};
   await context.pollCoreProgress(3);assert.equal(delay,1000);assert.equal(typeof callback,'function');
   assert.equal($('core-rule-status').textContent,api.t('u9c.coreProgress',{rule:context.coreRuleName('coverage_brics'),done:27,total:100}));
   let finish;context.getJSON=()=>new Promise(resolve=>finish=resolve);const pending=context.pollCoreProgress(3);context.switchCoreRule.busy=false;finish({rule:'coverage_brics',done:99,total:100});await pending;
   assert.equal(context.switchCoreRule.message.done,27,'过期进度响应不能覆盖完成状态');
  }
 },
 async list(){
  const {context,$,api,calls}=await setup();
  for(const language of ['en','zh']){
   api.setLanguage(language);await context.loadPreviousRuns();assert.equal($('previous-results-panel').hidden,false);assert.equal(calls.at(-1)[0],'/api/previous-runs');
   assert.equal($('previous-run-select').querySelectorAll('input').length,1);assert.equal($('previous-run-select').value,sample.run_name);
   const row=$('previous-run-select').querySelector('label'),label=row.textContent;for(const value of ['2026-11-12','CHEMBL205','1145','294',context.coreRuleName('coverage_brics')])assert.ok(label.includes(value),value);assert.ok(row.title.includes('CAII.csv'));
   assert.equal($('previous-results-unavailable').textContent,api.t('u9c.unavailable',{count:6,missing:2,incompatible:1,modified:3}));assert.equal($('open-selected-run').disabled,false);
   context.renderPreviousRuns({runs:Array.from({length:25},(_,i)=>({...sample,run_name:'r'+i})),unavailable:{}});assert.equal($('previous-run-select').querySelectorAll('input').length,20);
   context.renderPreviousRuns({runs:[],unavailable:{}});assert.equal($('previous-results-status').textContent,api.t('u9c.previousEmpty'));assert.equal($('open-selected-run').disabled,true);
  }
 },
 async restore(){
  const {context,$,api,calls}=await setup();
  const info={input:{filename:'CAII.csv'},query_input:{filename:'candidates.csv'},source_input:{filename:'sources.csv'},new_molecules:{options:{smiles_col:'smiles'}},job:{state:'done'},restored:{new_set:'sim_generated_molecules',total:294,manual:1,recorded:2}};
  for(const language of ['en','zh']){
   api.setLanguage(language);$('previous-run-select').value=sample.run_name;
   context.getJSON=async(url,options)=>{assert.equal(url,'/api/open-run');assert.deepEqual(JSON.parse(options.body),{run_name:sample.run_name});return info;};
   await context.openSelectedRun();assert.equal(context.busy,false);assert.equal(context.job,info.job);assert.equal(context.newMolecules,info.new_molecules);
   assert.ok(calls.includes('reset'));assert.ok(calls.includes('show'));assert.ok(calls.some(call=>Array.isArray(call)&&call[0]==='source'&&call[1]===info.source_input));
   assert.equal($('status').textContent,api.t('u9c.restoredNew',{name:'sim_generated_molecules',total:294,manual:1,recorded:2}));
   assert.equal(context.restoredSessionMessage({new_set:null}),api.t('u9c.restoredRun'));
   const oldJob=context.job;context.getJSON=async()=>{throw Error(api.t('u9c.error.session_modified',{files:'map.csv'}));};await context.openSelectedRun();assert.equal(context.job,oldJob);assert.equal(context.busy,false);assert.ok($('status').textContent.includes('map.csv'));
   context.busy=true;context.getJSON=async()=>{throw Error('must not run');};await context.openSelectedRun();context.busy=false;
  }
 },
 async markers(){
  const {context,$,api,node,calls}=await setup(),svg=node('svg'),wrapper=node('g'),cell=node('g'),rect=node('rect'),marker=node('path'),outlines=node('g');
  wrapper.setAttribute('transform','translate(200,100)');cell.setAttribute('transform','translate(-50,20)');cell.dataset.gridRow='2';cell.dataset.gridCol='3';rect.setAttribute('x','100');rect.setAttribute('y','200');rect.setAttribute('width','100');marker.setAttribute('data-neighbor-mismatch','true');marker.dataset.measuredClass='inactive';marker.dataset.neighborClass='active';
  cell.append(rect,marker);wrapper.append(cell);svg.append(wrapper,outlines);$('map-stage').append(svg);context.mapCells=[cell];
  cell.querySelector=selector=>selector===':scope > rect'?rect:selector.includes('data-neighbor-mismatch')?marker:null;
  for(const language of ['en','zh'])for(const scale of [.02,.1,.2,.3,.5,1,5]){
   api.setLanguage(language);context.transform.s=scale;context.renderMismatchOverlay();
   const size=context.mismatchMarkerSize(100,scale);assert.ok(size<=45);if(100*scale>=20)assert.ok(size*scale>=9);else assert.equal(size,45);
   assert.equal(marker.querySelector('title').textContent,api.t('u9c.mismatchTooltip',{measured:api.t('d.inactive'),neighbor:api.t('d.active')}));
   const overlay=svg.children.at(-1);assert.equal(overlay.dataset.mismatchOverlay,'true');assert.equal(overlay.children.length,1);assert.equal(overlay.children[0].getAttribute('transform'),'translate(200,100) translate(-50,20)');
   assert.equal(overlay.children[0].querySelector('title').textContent,marker.querySelector('title').textContent);overlay.children[0].onclick();assert.deepEqual(calls.at(-1),['inspect','row=2&col=3']);
   assert.equal(svg.querySelectorAll('[data-mismatch-overlay]').length,1);
  }
  const css=fs.readFileSync(path.join(process.argv[2],'app.css'),'utf8');assert.match(css,/data-mismatch-overlay\] path[^}]+pointer-events: visiblePainted/);
 },
 async controls(){
  const {context,$,api,document}=await setup();assert.equal(document.querySelectorAll('[data-help]').length,81);
  for(const language of ['en','zh']){
   api.setLanguage(language);
   for(const id of ['open-previous-results','previous-run-select','open-selected-run']){
    assert.ok($(id));const help=$(id+'-help');help.click();assert.equal($('parameter-tooltip').textContent,api.t(help.dataset.help));help.click();
    assert.match(api.t(help.dataset.help),language==='zh'?/本启动器输出目录.*版本不同.*重新生成.*磁盘空间.*整个运行文件夹/:/launcher.*version changes.*regeneration.*disk space.*entire run folder/);
   }
   for(const id of ['map-core-rule','core-rule','group-core-comparison'])assert.match(api.t('h.'+id),language==='zh'?/FMCS 始终串行/:/FMCS always remains serial/);
   assert.doesNotMatch(api.t('h.measured-result'),/不会自动恢复|do not automatically reload/);
   for(const code of ['session_missing','session_version','session_modified','session_unsafe','session_invalid','run_not_listed','map_running'])assert.ok(!api.t('u9c.error.'+code,{files:'map.csv',details:'Python'}).startsWith('u9c.'));
  }
 }
};
(async()=>{const name=process.argv[3];assert.ok(name in scenarios);await scenarios[name]();console.log('U9c '+name+': passed (English and Chinese)');})().catch(error=>{console.error(error);process.exitCode=1;});
