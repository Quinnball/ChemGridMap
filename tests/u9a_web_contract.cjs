'use strict';
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const {runtime,realFunction,dictionaries} = require('./u4_web_contract.cjs');
const language = process.argv[4], assets = process.argv[2];
const src = fs.readFileSync(path.join(assets,'app.js'),'utf8');
const sourceA = {id:'DOC / A',doc_type:'PATENT',journal:'Patent journal',year:2024,map_molecule_count:2,n_records:4,links:[]};
const sourceB = {id:'DOC B',doc_type:'PUBLICATION',journal:'Medicinal Chemistry',year:2025,map_molecule_count:1,n_records:1,links:[]};
const group = {structure_group:'G001',size:3,core_smiles:'c1ccccc1',core_rule:'coverage',core_heavy_atoms:6,core_coverage_count:3,core_coverage_fraction:1,core_coverage_target:.8,core_target_met:true,core_exact_scaffold_count:2};
const reference = {molecule_id:'R1',molecule_identity_key:'R1',grid_row:0,grid_col:0,structure_group:'G001'};
const second = {molecule_id:'R2',grid_row:0,grid_col:1,structure_group:'G001'};
const added = {molecule_id:'N1',query_id:'N1',grid_row:1,grid_col:0,is_new:true,structure_group:'G001',measurement_status:'not_measured'};
const newOther = {...added,molecule_id:'N2',query_id:'N2',grid_row:1,grid_col:1};
function comparison() {
  return {group:'G001',rule:'coverage',file:'core_rules/group_123456789abc_core_comparison.csv',fmcs_timeout:10,
    cores:['most_frequent','coverage','coverage_brics','fmcs'].map((rule,i)=>({...group,core_rule:rule,core_svg:'<svg data-rule="'+rule+'"></svg>',core_heavy_atoms:6+i,seconds:i*.1,fmcs_timed_out:rule==='fmcs'}))};
}
async function setup(){
  const env = await runtime();
  const $ = id=>env.document.getElementById(id);
  const calls=[],downloads=[],inspections=[];
  const job={options:{core_rule:'coverage',projection:'pca',input_format:'chembl',representation:'morgan'},current_core_rule:'coverage',files:[],report:{}};
  Object.assign(env.context,{$,t:env.api.t,translateServer:env.api.translateServer,job,structureGroups:[group],points:[reference,second,added,newOther],highlightedGroup:'',groupRequestNumber:0,requestNumber:0,currentCell:null,neighborRelations:new Map(),mapCells:[],view:'projection',busy:false,inputInfo:{},demoAvailable:true,umapAvailable:true,newMolecules:null,
    setView(value){env.context.view=value;},showResult(){},inspect(...args){inspections.push(args);},compareGroup(){},highlight(){},fit(){},renderLegend(){},updateNewNavigation(){},download(...args){downloads.push(args);},fail(error){throw error;},
    getJSON:async(url,options)=>{
      calls.push([url,options]);
      if(url.startsWith('/api/source-members?')) {
        if(url.includes('scope=new'))return {scope:'new',sources:[{id:'Batch 1',map_molecule_count:2}],members:[added,newOther],map_molecule_count:2};
        return {scope:'reference',sources:url.includes('source=DOC%20B')?[sourceA,sourceB]:[sourceA],members:[reference,second],map_molecule_count:2};
      }
      if(url.startsWith('/api/group?'))return {summary:env.context.structureGroups[0],core_svg:'<svg></svg>'};
      if(url.startsWith('/api/group-cores?'))return comparison();
      throw Error('Unexpected URL: '+url);
    }
  });
  const names=['format','percentage','countLabel','showSources','selectGroup','coreRule','updateNavigation','showQuerySelection','measuredLabel','table'];
  vm.runInContext(names.map(realFunction).join('\n'),env.context);
  env.api.setLanguage(language);env.context.showResult.grouping={groups:[group],protocol:{method:'Butina'}};
  env.context.renderCoreRuleOptions();
  // 真实事件绑定；VM 只模拟本任务访问的 DOM，不启动浏览器。
  for(const id of ['map-core-rule','group-core-comparison','clear-source-highlight']) {
    const line=src.split('\n').find(line=>line.startsWith(`$('${id}').addEventListener`));assert.ok(line,id);vm.runInContext(line,env.context);
  }
  return {...env,$,calls,downloads,inspections};
}
function map(env,tag='g'){
  const svg=env.node('svg');env.$('map-stage').replaceChildren(svg);
  const cells=env.context.points.map(point=>{
    const cell=env.node(tag);cell.tagName=tag;cell.dataset.gridRow=String(point.grid_row);cell.dataset.gridCol=String(point.grid_col);
    cell.cloneNode=()=>{const copy=env.node(tag);copy.tagName=tag;return copy;};
    if(tag==='g'){const rect=env.node('rect');rect.cloneNode=()=>env.node('rect');cell.append(rect);const original=cell.querySelector.bind(cell);cell.querySelector=selector=>selector===':scope > rect'?rect:original(selector);}
    svg.append(cell);return cell;
  });
  env.context.mapCells=cells;env.$('measurement-filter').value='all';env.$('show-neighbor-links').checked=false;
  vm.runInContext(realFunction('highlight'),env.context);return cells;
}
const scenarios={
  async sources(){
    const {context,$}=await setup();
    context.showSources([sourceA,sourceB],{molecule:reference,source_union_molecule_count:2,source_lookup_file:'lookup.csv',possible_sources:[sourceA]});
    assert.equal($('source-list').children.length,2);
    assert.equal($('highlight-source-0').tagName,'BUTTON');
    assert.equal($('highlight-source-0').textContent,language==='zh'?'高亮该来源的 2 个分子':'Highlight the 2 molecules from this source');
    assert.equal($('highlight-source-1').tagName,'P');
    assert.equal($('highlight-source-1').textContent,language==='zh'?'地图上只有本分子来自该来源':'Only this molecule on the map comes from this source.');
    assert.equal($('highlight-source-union').hidden,false);
    assert.equal($('highlight-source-union').textContent,language==='zh'?'高亮与本分子共享来源的分子（2 个）':'Highlight molecules sharing any of these sources (2)');
    assert.equal($('possible-source-list').querySelectorAll('button').length,0);
    context.showSources([sourceA]);assert.equal($('highlight-source-union').hidden,true);
    context.showSources();assert.equal($('sources').hidden,true);
  },
  async source_highlight(){
    const {context,$,calls}=await setup();context.highlightedGroup='G001';
    await context.highlightSources([sourceA.id]);
    assert.equal(context.highlightedGroup,'');assert.equal($('group-select').value,'');
    assert.equal($('source-highlight-status').hidden,false);
    assert.equal($('source-highlight-label').textContent,language==='zh'?'来源：DOC / A · 专利 · 2024 · 2 个分子':'Source: DOC / A · Patent · 2024 · 2 molecules');
    assert.ok(calls[0][0].includes('source=DOC%20%2F%20A'));
    await context.highlightSources([sourceA.id,sourceB.id],{molecule:'R1'});
    assert.equal($('source-highlight-label').textContent,language==='zh'?'与 R1 共享来源：2 个来源 · 2 个分子':'Sources shared with R1: 2 sources · 2 molecules');
    await context.selectGroup('G001');assert.equal(context.highlight.source,null);assert.equal($('source-highlight-status').hidden,true);assert.equal(context.highlightedGroup,'G001');
    await context.highlightSources([sourceA.id]);$('clear-source-highlight').click();assert.equal(context.highlight.source,null);assert.equal($('source-highlight-status').hidden,true);
    let finish;context.getJSON=()=>new Promise(resolve=>finish=resolve);const pending=context.highlightSources([sourceA.id]);
    context.clearSourceHighlight();finish({scope:'reference',sources:[sourceA],members:[reference],map_molecule_count:1});await pending;assert.equal(context.highlight.source,null,'清除后不应用过期请求');
  },
  async layers(){
    const env=await setup(),{context,$}=env;
    for(const tag of ['g','circle'])for(const mode of ['activity','support','reliability','structure','neighbor_evidence']){
      const cells=map(env,tag);$('color-mode').value=mode;
      context.highlight.source={cells:new Set(['0,0','0,1'])};context.currentCell={row:1,col:0};context.highlight();
      assert.ok(cells[0].classList.contains('source-highlight'),tag+' '+mode);assert.ok(!cells[0].classList.contains('group-faded'));
      assert.ok(cells[2].classList.contains('group-faded'));assert.ok(cells[2].classList.contains('selected-cell'));
      assert.ok($('map-stage').querySelector('svg').children.some(node=>node.dataset.mapOutlines));
      context.clearSourceHighlight();assert.ok(!cells[0].classList.contains('source-highlight'));assert.ok(!cells[2].classList.contains('group-faded'));
    }
    const cells=map(env,'circle');
    cells[0].classList.add('reference-evidence-faded');cells[1].classList.add('reference-evidence-faded');
    context.highlight.source={cells:new Set(['0,0','0,1'])};context.currentCell={row:0,col:1};
    context.neighborRelations=new Map([['0,0','representation_only']]);$('show-neighbor-links').checked=true;context.highlight();
    const outline=$('map-stage').querySelector('svg').children.find(node=>node.dataset.mapOutlines);
    assert.equal(outline.children.length,1,'来源成员的近邻框不覆盖来源虚线');assert.equal(outline.children[0].dataset.outline,'selected','保留选中分子实线框');
    const css=fs.readFileSync(path.join(assets,'app.css'),'utf8');
    assert.match(css,/#map-stage \.source-highlight>rect, #map-stage circle\.source-highlight[^}]+stroke-dasharray[^}]+!important/);
    assert.match(css,/selected-cell[^}]+stroke-dasharray: none !important/);
  },
  async new_batches(){
    const {context,$,calls}=await setup();context.renderNewBatchSources([{id:'Batch 1',map_molecule_count:2},{id:'Solo',map_molecule_count:1}]);
    assert.equal($('new-batch-sources').hidden,false);
    assert.equal($('highlight-new-batch-0').textContent,language==='zh'?'高亮同一批次的 2 个新分子':'Highlight the 2 new molecules from this batch');
    assert.equal($('highlight-new-batch-1').textContent,language==='zh'?'高亮同一批次的 1 个新分子':'Highlight the 1 new molecule from this batch');
    await $('highlight-new-batch-0').onclick();assert.ok(calls.at(-1)[0].endsWith('&scope=new'));
    assert.deepEqual([...context.highlight.source.cells],['1,0','1,1']);
    assert.equal(context.view,'grid');
    assert.equal($('source-highlight-label').textContent,language==='zh'?'批次：Batch 1 · 2 个新分子 · 已切换到网格图；新分子仅有网格位置。':'Batch: Batch 1 · 2 new molecules · Switched to the grid; new molecules have grid positions only.');
    context.renderNewBatchSources([]);assert.equal($('new-batch-sources').hidden,true);
  },
  async rules(){
    const {context,$,document}=await setup();
    assert.deepEqual($('map-core-rule').options.map(o=>o.value),['most_frequent','coverage','coverage_brics']);
    assert.equal($('map-core-rule').value,'coverage');assert.ok($('map-core-rule').closest('[data-group-only]'));
    context.job.options.core_rule='fmcs';context.job.current_core_rule='fmcs';context.renderCoreRuleOptions();
    assert.deepEqual($('map-core-rule').options.map(o=>o.value),['most_frequent','coverage','coverage_brics','fmcs']);
    assert.equal($('map-core-rule').options.at(-1).textContent,language==='zh'?'FMCS（建图时计算）':'FMCS (computed when building the map)');
    assert.equal($('map-core-rule').value,'fmcs');
    assert.ok($('core-rule').options.some(o=>o.value==='coverage_brics'));
    assert.equal(document.querySelectorAll('[data-help]').length,60);
    assert.equal($('map-core-rule-help').dataset.help,'h.map-core-rule');
    assert.equal($('group-core-comparison-help').dataset.help,'h.group-core-comparison');
    assert.ok(dictionaries[language]['h.core-rule'].includes('BRICS'));
    assert.ok(dictionaries[language]['h.group-core-comparison'].includes('10'));
    assert.ok(context.coreRule({...group,core_rule:'coverage_brics'}).includes('BRICS'));
  },
  async busy(){
    const {context,$,inspections}=await setup();let finish;const oldGet=context.getJSON;
    context.getJSON=(url,options)=>url==='/api/core-rule'?new Promise(resolve=>finish=resolve):oldGet(url,options);
    context.highlightedGroup='G001';context.currentCell={row:0,col:0};context.inspect.lastDetail={molecule:reference};
    $('map-core-rule').value='coverage_brics';const pending=context.switchCoreRule();
    assert.equal(context.busy,false);assert.equal(context.switchCoreRule.busy,true);
    assert.equal($('map-core-rule').disabled,true);assert.equal($('compare-group').disabled,true);assert.equal($('selected-group-button').disabled,true);
    assert.equal($('fit').disabled,false);assert.equal($('zoom-in').disabled,false);assert.equal($('grid-tab').disabled,false);
    context.renderCoreRuleOptions();assert.equal($('map-core-rule').value,'coverage_brics','语言重绘期间保持正在请求的规则');
    assert.equal($('core-rule-status').textContent,language==='zh'?'正在按「覆盖目标比例的最大核心骨架（含 BRICS 片段）」计算 1 个组的核心…':'Computing cores for 1 groups using “Largest core covering the target share, including BRICS fragments”…');
    const changed={...group,core_rule:'coverage_brics',core_coverage_count:2,core_coverage_fraction:2/3};
    finish({rule:'coverage_brics',groups:[changed],files:['core_rules/run_cores_coverage_brics.csv']});await pending;
    assert.equal(context.job.current_core_rule,'coverage_brics');assert.equal(context.switchCoreRule.busy,false);assert.equal($('map-core-rule').disabled,false);assert.equal($('compare-group').disabled,false);
    assert.ok($('group-select').options.at(-1).textContent.includes('2/3'));assert.ok($('group-core-rule').textContent.includes('BRICS'));
    assert.equal(inspections.at(-1)[0],null,'参照分子只从缓存重绘，避免重写审计ZIP');assert.equal(inspections.at(-1)[1].group_core.core_rule,'coverage_brics');
    assert.ok(context.job.files.includes('core_rules/run_cores_coverage_brics.csv'));
    context.inspect.lastDetail={molecule:added};context.currentCell={row:1,col:0};$('map-core-rule').value='coverage';
    const next=context.switchCoreRule();finish({rule:'coverage',groups:[group],files:[]});await next;assert.equal(inspections.at(-1)[0],'row=1&col=0');
  },
  async comparison(){
    const {context,$,downloads}=await setup();context.highlightedGroup='G001';$('group-core-comparison').open=true;
    await context.loadCoreComparison();const rows=$('core-comparison-table').querySelectorAll('tbody tr');
    assert.equal(rows.length,4);assert.deepEqual(rows.map(row=>row.dataset.coreRule),['most_frequent','coverage','coverage_brics','fmcs']);
    assert.equal(rows[1].dataset.currentRule,'true');assert.equal(rows[0].children[4].textContent,'—');
    assert.equal(rows[2].children[3].textContent,'3/3 (100%)');assert.ok(rows.every(row=>row.children[1].innerHTML.includes('<svg')));
    assert.ok(rows[3].textContent.includes(language==='zh'?'超时，结果可能不完整':'Timed out; the result may be incomplete'));
    assert.equal($('core-comparison-status').textContent,language==='zh'?'FMCS 搜索预算：每组 10 秒。':'FMCS search budget: 10 s per group.');
    assert.equal($('download-core-comparison').hidden,false);$('download-core-comparison').click();assert.equal(downloads[0][1],comparison().file);
    context.job.current_core_rule='coverage_brics';context.renderCoreComparison(context.loadCoreComparison.lastResult);
    const next=$('core-comparison-table').querySelectorAll('tbody tr');assert.equal(next[2].dataset.currentRule,'true');assert.equal(next[1].dataset.currentRule,'false');
    context.showQuerySelection({molecule:{...added,contains_group_core:true,contains_group_core_rule:'coverage_brics'}});
    assert.ok($('query-summary').textContent.includes('BRICS'));
  },
  async stale_comparison(){
    const {context,$}=await setup();vm.runInContext(realFunction('compareGroup'),context);
    for(const change of ['rule','run']){
      let finish;context.getJSON=()=>new Promise(resolve=>finish=resolve);
      const pending=context.compareGroup('G001');
      if(change==='rule')context.job.current_core_rule='coverage_brics';else context.job={...context.job};
      finish({summary:{...group},drawings:{}});await pending;
      assert.equal(context.compareGroup.lastResult,undefined);assert.equal($('comparison-dialog').open,false);
    }
  },
  async language(){
    const {context,$,api}=await setup();await context.highlightSources([sourceA.id]);context.renderCoreComparison(comparison());
    const other=language==='zh'?'en':'zh';api.setLanguage(other);context.renderSourceHighlight();context.renderCoreRuleOptions();context.renderCoreComparison(comparison());
    assert.ok($('source-highlight-label').textContent.startsWith(other==='zh'?'来源：':'Source:'));
    assert.equal($('core-comparison-status').textContent,dictionaries[other]['u9.fmcsBudget']);assert.equal($('clear-source-highlight').textContent,dictionaries[other]['u9.clearHighlight']);
  },
  async errors(){
    const {context,$,api}=await setup();
    for(const key of ['server.unknownSource','server.chooseSources','server.buildMapFirst','server.unavailableCoreRule','server.noStructureGroups','server.coreMapChanged']){
      const english=dictionaries.en[key].replace('{id}','X');assert.equal(api.translateServer(english),dictionaries[language][key].replace('{id}','X'));
    }
    const errors=[];context.fail=error=>errors.push(error.message);context.getJSON=async()=>{throw Error('Unknown or unavailable core rule.');};
    $('map-core-rule').value='coverage_brics';await context.switchCoreRule();assert.equal($('map-core-rule').value,'coverage');assert.equal(context.switchCoreRule.busy,false);assert.equal(errors.length,1);assert.equal($('core-rule-status').hidden,true);
  },
  async reset(){
    const {context,$}=await setup();await context.highlightSources([sourceA.id]);assert.equal($('source-highlight-status').hidden,false);
    // 建图/插入/移除共享 showResult 或 resetResult 的清理入口。
    assert.match(realFunction('showResult'),/if \(!languageOnly\) \{ clearSourceHighlight\(\)/);
    assert.match(realFunction('resetResult'),/clearSourceHighlight\(\)/);
    assert.match(realFunction('insertNewMolecules'),/await showResult\(\)/);
    assert.match(realFunction('removeQueryInput'),/await showResult\(\)/);
    context.clearSourceHighlight();assert.equal(context.highlight.source,null);
  }
};
assert.ok(['en','zh'].includes(language));
scenarios[process.argv[3]]().then(()=>console.log(`U9a ${process.argv[3]} ${language}: passed`)).catch(error=>{console.error(error);process.exitCode=1;});
