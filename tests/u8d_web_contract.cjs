'use strict';
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const {runtime, realFunction} = require('./u4_web_contract.cjs');
const language = process.argv[4];
const description = 'CHEMBL205 · Carbonic anhydrase 2 · IC50';
const reference = {molecule_id:'R1', grid_row:0, grid_col:0, class_status:'supported', neighbor_recall:.8, activity_class:'active'};
const added = {molecule_id:'N1',query_id:'N1',grid_row:0,grid_col:1,is_new:true,status:'mapped',measurement_status:'not_measured',neighbor_evidence:'active'};
const input = {filename:'new.csv',rows:289,columns:['smiles','id'],detected:{smiles:'smiles'}};
const warningFixtures = [
  {code:'structural_mismatch',placed:10,none:9,median:.25},
  {code:'target_mismatch',targets:['CHEMBL204','CHEMBL240'],reference:'CHEMBL205'}
];
const expected = {
  zh:{counts:'近邻一致·活性 70、近邻证据不足 4、无相近的参照分子 215',
    inserted:counts=>`已插入 289 个新分子，参照分子的位置未改变。近邻证据：${counts}。已切换到“近邻证据”图层，外框颜色表示近邻证据。`,
    reference:description=>`参照图：${description} · 参照分子 900 个`,
    audit:(description,file)=>`参照图：${description}；新分子文件：${file}`,
    legend:'（近邻证据见“近邻证据”图层）',
    warnings:['大多数新分子与参照分子差别很大：10 个中有 9 个没有相近的参照分子，与最相似参照分子的 Tanimoto 中位数为 0.25。请确认新分子与参照图是同一个靶点。','新分子文件中的靶点（CHEMBL204、CHEMBL240）与参照图的靶点（CHEMBL205）不同。']},
  en:{counts:'Consistent neighbors: active 70, Too little neighbor evidence 4, No similar reference molecule 215',
    inserted:counts=>`Inserted 289 new molecules; reference molecules have not moved. Neighbor evidence: ${counts}. Switched to the Neighbor evidence layer; outline colors show the neighbor evidence.`,
    reference:description=>`Reference map: ${description} · 900 reference molecules`,
    audit:(description,file)=>`Reference map: ${description}; new-molecule file: ${file}`,
    legend:'(see the Neighbor evidence layer)',
    warnings:['Most new molecules differ strongly from the reference molecules: 9 of 10 have no similar reference molecule, and the median Tanimoto to the most similar reference molecule is 0.25. Check that the new molecules and the reference map concern the same target.',"The new-molecule file lists targets (CHEMBL204, CHEMBL240) that differ from the reference map's target (CHEMBL205)."]}
};
async function setup(){
  const env = await runtime();
  const $ = id=>{const element=env.document.getElementById(id);assert.ok(element,id);return element;};
  const report={input_rows:1000,retained_activity_records:900,repeated_conflicted_molecules:0,warnings:[],query_molecules:{input_rows:289,mapped:289,invalid_smiles:0,duplicate_query:0,in_reference:0,neighbor_evidence_counts:{none:215,mixed:0,active:70,insufficient:4,inactive:0,medium:0},warnings:[]}};
  const job={state:'done',options:{input_format:'chembl',target_id:'CHEMBL205',representation:'morgan',projection:'pca',grid_occupancy:.7,lower_threshold:6,upper_threshold:7},reference_map:{description,count:900,target_id:'CHEMBL205'},report,molecules:1189,files:['map.svg'],metrics:{},output_dir:'/tmp/example'};
  Object.assign(env.context,{$,t:env.api.t,translateServer:env.api.translateServer,job,points:[reference,added],structureGroups:[],highlightedGroup:'',newMolecules:{input,rows:[added],report:report.query_molecules,files:['map_with_new_molecules.svg']},newListSort:{key:'query_id',direction:1},queryInputInfo:input,inputInfo:{filename:'reference.csv',rows:1000,columns:[],preview:[]},sourceInputInfo:null,umapAvailable:true,busy:false,view:'grid',currentCell:null,requestNumber:0,mapCells:[],
    colors:{active:'#2a9d8f',medium:'#e6a23c',inactive:'#c44e52',unlabelled:'#b8c2cc'},fadedColors:{active:'#b8e0db'},reliabilityBins:[[0,'#cde1f5']],
    setBusy(value){env.context.busy=value;},closeInspector(){},highlight(){},fit(){},updateInputType(){},applyTransform(){},translateSVGText(){},refreshValidation(){},restoreOptions(){},updateNavigation(){},showPage(){},inspect(){},selectGroup(){},compareGroup(){},
    setView(value){env.context.view=value;env.context.colorMap();},
    api:async()=>({text:async()=>'<svg></svg>'}),
    getJSON:async(url,options)=>{
      if(url==='/api/new-molecules/clear'){
        delete env.context.job.report.query_molecules;
        env.context.job.molecules=env.context.job.reference_map.count;
        env.context.newMolecules={rows:[],report:{},options:{},files:[]};env.context.points=[reference];
        return {job:env.context.job};
      }
      if(url==='/api/new-molecules')return options?.method==='POST'?env.context.job:env.context.newMolecules;
      if(url==='/api/points')return env.context.points;
      if(url==='/api/groups')return {groups:[],protocol:{method:'none'}};
      throw Error(`Unexpected URL ${url}`);
    }
  });
  const names=['updateNewNavigation','referenceMapDetails','renderReferenceMap','newMoleculeWarning','newMoleculeInsertedMessage','selectColumns','showQueryInput','newMoleculeOptions','insertNewMolecules','removeQueryInput','bindNewMolecules','newListColumns','visibleNewMolecules','newListValue','renderNewList','downloadNewList','format','stat','countLabel','percentage','table','measuredLabel','showDataAudit','showResult','status','fail','legendItem','renderColorLegend','renderLegend','colorMap','reliabilityColor','refreshLanguage'];
  vm.runInContext(names.map(realFunction).join('\n'),env.context);
  env.api.setLanguage(language);env.context.showQueryInput(input);env.context.bindNewMolecules();
  return {...env,$};
}
function auditReference($){return $('audit-body').querySelectorAll('p').find(item=>item.dataset.newMoleculesReference);}
async function insert(context){await context.insertNewMolecules({preventDefault(){}});}
const scenarios={
  async workflow(){
    const {context,$,api}=await setup();
    $('color-mode').value='support';await insert(context);
    assert.equal($('color-mode').value,'neighbor_evidence');
    assert.equal($('new-molecules-dialog').open,false);
    $('color-mode').value='reliability';await context.showResult();
    assert.equal($('color-mode').value,'reliability','再次显示已有结果不重置图层');
    await context.refreshLanguage();assert.equal($('color-mode').value,'reliability');
    const other=language==='zh'?'en':'zh';api.setLanguage(other);await context.refreshLanguage();
    assert.equal($('color-mode').value,'reliability','语言切换不重置图层');
    assert.equal($('status').textContent,expected[other].inserted(expected[other].counts));
    api.setLanguage(language);
    const savedJob=context.job;context.job=null;context.busy=true;await context.refreshLanguage();
    assert.equal($('color-mode').value,'reliability','重新建图运行中切换语言不改图层');assert.equal($('new-list-reference').textContent,'');
    context.job=savedJob;context.busy=false;
    await insert(context);assert.equal($('color-mode').value,'neighbor_evidence','下一次成功插入再次切换');
    await context.removeQueryInput();assert.equal($('color-mode').value,'activity');
    assert.equal($('color-mode').querySelector('[value=neighbor_evidence]').disabled,true);
    $('color-mode').value='support';await context.showResult();assert.equal($('color-mode').value,'support','重新建图不自动切到近邻证据');
    await context.removeQueryInput();assert.equal($('color-mode').value,'support','移除时其他图层保留');
    context.queryInputInfo=input;context.getJSON=async()=>{throw Error('Insertion failed');};
    await insert(context);assert.equal($('color-mode').value,'support','失败插入不切换');assert.equal($('status').className,'error');
  },
  async zeroPlaced(){
    const {context,$}=await setup();
    delete context.job.report.query_molecules;context.newMolecules={input,rows:[],report:{},options:{},files:[]};context.points=[reference];
    await context.showResult();assert.equal($('color-mode').querySelector('[value=neighbor_evidence]').disabled,true,'仅上传待插入文件不启用近邻证据');
    context.job.report.query_molecules={input_rows:1,mapped:0,neighbor_evidence_counts:{active:0,medium:0,inactive:0,mixed:0,insufficient:0,none:0},warnings:[]};
    context.newMolecules={input,rows:[],files:['map_with_new_molecules.svg'],report:context.job.report.query_molecules};context.points=[reference];
    await insert(context);assert.equal($('color-mode').value,'neighbor_evidence');assert.equal($('color-mode').querySelector('[value=neighbor_evidence]').disabled,false);
    assert.ok($('status').textContent.includes(language==='zh'?'近邻证据：0。':'Neighbor evidence: 0.'));
    $('color-mode').value='activity';await context.showResult();assert.equal($('color-mode').value,'activity','恢复0个放置结果也不自动切换');
  },
  async counts(){
    const {context,$}=await setup();await insert(context);
    assert.equal($('status').textContent,expected[language].inserted(expected[language].counts));
    const all={none:6,insufficient:5,mixed:4,inactive:3,medium:2,active:1};context.job.report.query_molecules.neighbor_evidence_counts=all;await insert(context);
    const ordered=language==='zh'?'近邻一致·活性 1、近邻一致·中等活性 2、近邻一致·非活性 3、近邻不一致 4、近邻证据不足 5、无相近的参照分子 6':'Consistent neighbors: active 1, Consistent neighbors: medium 2, Consistent neighbors: inactive 3, Conflicting neighbors 4, Too little neighbor evidence 5, No similar reference molecule 6';
    assert.equal($('status').textContent,expected[language].inserted(ordered));assert.equal($('status').className,'');
  },
  async legends(){
    const {context,$,document}=await setup();
    for(const mode of ['activity','support','reliability']){
      context.renderLegend(mode);
      const labels=$('legend').children.filter(item=>item.textContent.includes(expected[language].legend));
      assert.equal(labels.length,mode==='activity'?2:1);
      for(const label of labels)assert.ok(label.textContent.endsWith(expected[language].legend));
    }
    context.renderLegend('neighbor_evidence');assert.equal($('legend').children.filter(item=>item.dataset.evidenceCategory).length,6);
    assert.ok(!$('legend').textContent.includes(expected[language].legend));
    assert.equal(document.querySelectorAll('[data-help]').length,80);
  },
  async reference(){
    const {context,$}=await setup();
    for(const text of [description,'CHEMBL205 · IC50','laboratory.csv · potency']){
      context.job.reference_map.description=text;await context.showResult();
      assert.ok($('run-caption').textContent.startsWith(text+' · '));
      $('open-new-molecules').click();assert.equal($('new-molecules-reference').textContent,expected[language].reference(text));
      $('open-new-list').click();assert.equal($('new-list-reference').textContent,expected[language].reference(text));
      const line=$('audit-body').children.find(item=>item.tagName==='SECTION' && item.children[0]?.textContent===(language==='zh'?'新分子':'New molecules')).children[1];
      assert.equal(line.textContent,expected[language].audit(text,'new.csv'));
      context.queryInputInfo={...input,filename:'not_inserted_yet.csv'};context.showDataAudit(context.job.report);
      assert.equal(auditReference($).textContent,expected[language].audit(text,'new.csv'),'审计显示已插入文件，不借用待插入文件名');
    }
  },
  async warnings(){
    const {context,$,api}=await setup();
    for(const warnings of [[warningFixtures[0]],[warningFixtures[1]],warningFixtures]){
      context.job.report.query_molecules.warnings=warnings;await insert(context);
      assert.equal($('color-mode').value,'neighbor_evidence');assert.equal($('status').className,'warning');assert.equal($('status').hidden,false);
      const texts=warnings.map(warning=>expected[language].warnings[warningFixtures.indexOf(warning)]);
      assert.equal($('status').textContent,[expected[language].inserted(expected[language].counts),...texts].join(' '));
      assert.deepEqual($('audit-body').querySelectorAll('.new-molecule-warning').map(item=>item.textContent),texts);
    }
    const other=language==='zh'?'en':'zh';api.setLanguage(other);await context.refreshLanguage();
    assert.equal($('status').textContent,[expected[other].inserted(expected[other].counts),...expected[other].warnings].join(' '));
    assert.deepEqual($('audit-body').querySelectorAll('.new-molecule-warning').map(item=>item.textContent),expected[other].warnings);
    assert.equal($('status').className,'warning');
    const css=fs.readFileSync(path.join(process.argv[2],'app.css'),'utf8');
    assert.match(css,/#status\.warning,\.new-molecule-warning\s*\{[^}]*#76521b[^}]*#fff7ea/);
    assert.match(css,/#status\.error\s*\{[^}]*#962f37/);
  },
  async payload(){
    const {context,$}=await setup();const payload=JSON.parse(fs.readFileSync(process.argv[5],'utf8'));
    Object.assign(context,{job:payload.job,newMolecules:payload.new,points:payload.points,queryInputInfo:payload.new.input});
    await insert(context);
    const details=payload.job.reference_map;
    assert.equal(details.description,description);assert.equal(details.count,4);
    assert.ok($('run-caption').textContent.startsWith(description+' · '));
    const ref=language==='zh'?`参照图：${description} · 参照分子 4 个`:`Reference map: ${description} · 4 reference molecules`;
    assert.equal($('new-molecules-reference').textContent,ref);$('open-new-list').click();assert.equal($('new-list-reference').textContent,ref);
    assert.equal(auditReference($).textContent,expected[language].audit(description,'silicon.csv'));
    assert.deepEqual(payload.new.report.warnings.map(w=>w.code),['structural_mismatch','target_mismatch']);
    assert.equal($('status').className,'warning');assert.equal($('color-mode').value,'neighbor_evidence');
    const auditWarnings=$('audit-body').querySelectorAll('.new-molecule-warning').map(item=>item.textContent);
    assert.equal(auditWarnings.length,2);for(const text of auditWarnings){assert.ok($('status').textContent.includes(text));assert.ok(text);assert.ok(!text.includes("[object Object]"));}
    assert.ok(auditWarnings[1].includes('CHEMBL204'));assert.ok(auditWarnings[1].includes('CHEMBL205'));
  }
};
assert.ok(['en','zh'].includes(language));
scenarios[process.argv[3]]().then(()=>console.log(`U8d ${process.argv[3]} ${language}: passed`)).catch(error=>{console.error(error);process.exitCode=1;});
