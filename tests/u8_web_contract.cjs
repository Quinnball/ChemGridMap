'use strict';
const assert = require('node:assert/strict');
const vm = require('node:vm');
const fs = require('node:fs');
const path = require('node:path');
const {runtime,realFunction} = require('./u4_web_contract.cjs');
const reference = {molecule_id:'R1',molecule_identity_key:'R1',grid_row:0,grid_col:0,activity_class:'active',class_status:'supported',projection_x:.5,projection_y:.5};
const added = {is_new:true,query_id:'N2',molecule_id:'N2',grid_row:0,grid_col:1,status:'mapped',label:'round 2',measurement_status:'measured',activity_class:'inactive',activity_pchembl:5,pchembl_min:4.8,pchembl_max:5.2,n_records:3,n_evidence_units:2,nearest_measured_id:'R1',nearest_measured_similarity:.9,nearest_reference_id:'R1',tanimoto:.9,placement_distance:1,neighbor_evidence:'active',similar_reference_count:3,supported_reference_count:2,supported_reference_classes:'active',reference_support_assessed:true,neighbor_comparison:'disagrees',structure_group:'G1'};
const input={filename:'new.csv',rows:4,columns:['smiles','gen_id','value','units','relation','batch','label'],detected:{smiles:'smiles'}};
const rows=[added,{...added,query_id:'N1',label:'round 1',neighbor_evidence:'mixed',tanimoto:.7},{...added,query_id:'N3',label:'round 1',status:'in_reference',neighbor_evidence:''},{...added,query_id:'N4',status:'invalid_smiles',neighbor_evidence:''}];
async function setup() {
  const env=await runtime(), requests=[], downloaded=[];
  const $=id=>{const item=env.document.getElementById(id);assert.ok(item,id);return item;};
  const job={state:'done',options:{representation:'morgan',input_format:'chembl',group_similarity:.6,lower_threshold:6,upper_threshold:7},files:['map.svg'],report:{},molecules:1};
  Object.assign(env.context,{$,t:env.api.t,translateServer:env.api.translateServer,job,points:[reference,added],structureGroups:[],highlightedGroup:'',newMolecules:{rows,files:['map_with_new_molecules.svg']},newListSort:{key:'query_id',direction:1},queryInputInfo:null,inputInfo:null,sourceInputInfo:null,umapAvailable:true,busy:false,view:'grid',currentCell:null,requestNumber:0,mapCells:[],
    colors:{active:'#2a9d8f',medium:'#e6a23c',inactive:'#c44e52',unlabelled:'#b8c2cc'},fadedColors:{active:'#b8e0db'},reliabilityBins:[[0,'#cde1f5']],
    setBusy(value){env.context.busy=value;env.context.updateNewNavigation();},showResult:async()=>{},closeInspector(){},highlight(){},fit(){},updateInputType(){},applyTransform(){},translateSVGText(){},refreshValidation(){},
    inspect(query){requests.push({query});},setView(value){env.context.view=value;},
    getJSON:async(url,options)=>{requests.push({url,options});return url.endsWith('/input')?input:url.endsWith('/clear')?{job}:job;},
    Blob,URL:{createObjectURL(blob){downloaded.push(blob);return 'blob:test';},revokeObjectURL(){}},
  });
  const names=['updateNewNavigation','showQueryInput','selectColumns','loadQueryInput','newMoleculeOptions','insertNewMolecules','removeQueryInput','bindNewMolecules','newListColumns','visibleNewMolecules','newListValue','renderNewList','downloadNewList','centerOnCurrentCell','format','countLabel','percentage','table','measuredLabel','showQuerySelection','referenceMapDetails', 'renderReferenceMap', 'newMoleculeWarning', 'newMoleculeInsertedMessage', 'showDataAudit','showSources','status','fail','legendItem','renderColorLegend','renderLegend','colorMap','reliabilityColor','projectionSVG'];
  vm.runInContext(names.map(realFunction).join('\n'),env.context);
  env.context.bindNewMolecules();
  return {...env,$,requests,downloaded};
}
const scenarios={
  async dialogs(){
    const {context,$,api,document}=await setup();
    assert.equal($('query-csv').closest('form').id,'new-molecules-form');
    assert.equal($('query-smiles-col').hasAttribute('required'),true);
    assert.equal($('options').querySelector('#query-csv'),null);
    assert.equal(document.querySelectorAll('[data-help]').length,77);
    for(const language of ['en','zh']){
      api.setLanguage(language);context.updateNewNavigation();$('open-new-molecules').click();assert.equal($('new-molecules-dialog').open,true);
      assert.equal($('open-new-molecules').textContent,api.t('u8.add'));$('new-molecules-dialog').close();
      context.job.options.representation='coordinates';context.updateNewNavigation();assert.equal($('open-new-molecules').disabled,true);assert.equal($('new-molecules-availability').textContent,api.t('u8.morganOnly'));
      context.job.options.representation='morgan';context.showQueryInput(input);context.updateNewNavigation();assert.equal($('insert-new-molecules').disabled,false);assert.equal($('query-id-col').value,'gen_id');
      assert.equal($('new-unit').disabled,true);$('new-activity-col').value='value';$('new-activity-col').dispatchEvent({type:'change'});assert.equal($('new-unit').disabled,false);
      $('new-unit-col').value='units';$('new-unit-col').dispatchEvent({type:'change'});assert.equal($('new-unit').disabled,true);
      for(const id of ['query-csv','query-smiles-col','query-id-col','new-activity-col','new-unit','new-unit-col','new-relation-col','new-document-col','new-label-col','new-class-filter','new-label-filter']){
        const button=document.querySelector(`[data-help-control=${id}]`);button.click();assert.equal($('parameter-tooltip').textContent,api.t(button.dataset.help));button.click();
      }
    }
    context.job=null;context.updateNewNavigation();assert.equal($('open-new-molecules').disabled,true);
  },
  async workflow(){
    const {context,$,requests}=await setup();const original=context.job;
    $('query-csv').files=[{name:'new.csv'}];await context.loadQueryInput();assert.equal(context.job,original);assert.equal(requests[0].url,'/api/new-molecules/input');
    $('new-activity-col').value='value';$('new-unit').value='nM';$('new-unit-col').value='units';$('new-document-col').value='batch';$('new-relation-col').value='relation';$('new-label-col').value='label';
    await context.insertNewMolecules({preventDefault(){}});const sent=JSON.parse(requests.at(-1).options.body);assert.equal(sent.unit,'p');assert.equal(sent.activity_col,'value');assert.equal(sent.unit_col,'units');assert.equal(sent.document_col,'batch');assert.equal(requests.at(-1).url,'/api/new-molecules');assert.equal(context.job,original);
    context.showQueryInput(input,{value_col:'value',value_unit:'nM',unit_column:'units',document_column:'batch'});assert.equal($('new-activity-col').value,'value');assert.equal($('new-document-col').value,'batch');
    await context.removeQueryInput();assert.equal(context.job,original);assert.equal(context.queryInputInfo,null);assert.equal(context.newMolecules,null);assert.equal(requests.at(-1).url,'/api/new-molecules/clear');
  },
  async list(){
    const {context,$,api,requests}=await setup();
    for(const language of ['en','zh']){
      api.setLanguage(language);$('new-class-filter').value='';$('new-label-filter').value='';$('open-new-list').click();assert.equal($('new-list-dialog').open,true);
      assert.equal($('new-list-table').querySelectorAll('tbody tr').length,4);
      assert.equal($('new-list-counts').children.length,6);
      const headers=$('new-list-table').querySelectorAll('th');assert.equal(headers.length,12);headers[0].querySelector('button').click();assert.equal($('new-list-table').querySelectorAll('tbody tr')[0].dataset.queryId,context.newListSort.direction<0?'N4':'N1');
      $('new-label-filter').value='round 1';$('new-label-filter').dispatchEvent({type:'change'});assert.equal($('new-list-table').querySelectorAll('tbody tr').length,2);
      $('new-class-filter').value='mixed';$('new-class-filter').dispatchEvent({type:'change'});const row=$('new-list-table').querySelector('tbody tr');assert.equal(row.dataset.queryId,'N1');row.click();assert.equal(requests.at(-1).query,'id=N1');assert.equal($('new-list-dialog').open,false);
      $('new-class-filter').value='';$('new-label-filter').value='';context.renderNewList();const invalid=$('new-list-table').querySelectorAll('tbody tr').find(row=>row.dataset.queryId==='N4');invalid.click();assert.equal($('new-list-note').textContent,api.t('u8.invalidNote'));
    }
  },
  async csv(){
    const {context,$,downloaded}=await setup();context.renderNewList(true);$('new-class-filter').value='mixed';$('download-new-list').click();const csv=await downloaded[0].text();assert.match(csv, /"query_id","label","status"/);assert.match(csv,/"N1"/);assert.doesNotMatch(csv,/"N2"|"N3"|"N4"/);
    context.newMolecules.rows=[{...added,label:'round, "one"\nnext'}];$('new-class-filter').value='';context.downloadNewList();assert.match(await downloaded[1].text(),/"round, ""one""\nnext"/);
  },
  async layers(){
    const {context,$,api,node}=await setup();
    for(const language of ['en','zh']){api.setLanguage(language);context.renderLegend('neighbor_evidence');assert.equal($('legend').children.filter(n=>n.dataset.evidenceCategory).length,6);for(const key of ['active','medium','inactive','mixed','insufficient','none'])assert.ok($('legend').textContent.includes(api.t('u8.'+key)));}
    const cells=[reference,{...added,measurement_status:'not_measured'},{...added,grid_col:2}].map(point=>{
      const cell=node('g'),mark=node('rect');cell.dataset.gridRow=point.grid_row;cell.dataset.gridCol=point.grid_col;mark.setAttribute('fill','#abc');cell.append(mark);cell.querySelector=selector=>selector===':scope > rect'?mark:null;return {cell,mark,point};});
    context.points=cells.map(c=>c.point);$('map-stage').querySelectorAll=()=>cells.map(c=>c.cell);
    for(const mode of ['activity','support','reliability','neighbor_evidence']){$('color-mode').value=mode;context.colorMap();
      assert.equal(cells[0].cell.classList.contains('reference-evidence-faded'),mode==='neighbor_evidence');
      if(['support','reliability'].includes(mode))for(const c of cells.slice(1))assert.equal(c.mark.getAttribute('fill'),'#dfe4e9');
      if(mode==='activity'){assert.equal(cells[1].mark.getAttribute('fill'),'#FFFFFF');assert.equal(cells[2].mark.getAttribute('fill'),'#c44e52');}
      if(mode==='neighbor_evidence'){assert.equal(cells[2].mark.style['--evidence-color'],'#2a9d8f');assert.equal(cells[2].mark.style['--evidence-width'],'3');}
    }
    for(const [key,width,dash] of [['mixed','3','6 4'],['insufficient','3','1 4'],['none','1.5','none']]){cells[2].point.neighbor_evidence=key;context.colorMap();assert.equal(cells[2].mark.style['--evidence-width'],width);assert.equal(cells[2].mark.style['--evidence-dash'],dash);}
    const css=fs.readFileSync(path.join(process.argv[2],'app.css'),'utf8');assert.equal(Number(css.match(/\.reference-evidence-faded\s*\{\s*opacity:\s*([.\d]+)/)[1]),.25);assert.match(css,/new-molecule-outline>rect[^}]*#173e34[^}]*1\.5 \* var\(--zoom-inverse/);assert.match(css,/new-evidence-outline>rect[^}]*var\(--evidence-width, 3\) \* var\(--zoom-inverse/);
  },
  async inspector(){
    const {context,$,api}=await setup();
    for(const language of ['en','zh']){api.setLanguage(language);context.showQuerySelection({molecule:added,query_neighbors:[{measured_id:'R1',rank:1,tanimoto:.9}],records:[{document_chembl_id:'R1',raw_activity_value:10,raw_activity_unit:'nM',relation:'='}]});
      assert.equal($('new-placement').textContent,api.t('u8.placement',{id:'R1',distance:'1.0'}));assert.match($('query-summary').textContent,language==='en'?/Disagrees with similar reference molecules/:/与相近的参照分子不一致/);
      assert.equal($('new-records').querySelectorAll('tbody tr').length,1);assert.ok($('query-neighbors').textContent.includes(api.t('u8.disclaimer')));
      context.showQuerySelection({molecule:{...added,reference_support_assessed:false},query_neighbors:[],records:[]});assert.equal($('new-support-note').hidden,false);assert.equal($('new-support-note').textContent,api.t('u8.supportNotAssessed'));
    }
    vm.runInContext(realFunction('inspect'),context);const measuredDetail={molecule:reference,query_match:{query_id:'CTRL',status:'in_reference',reference_id:'R1'},new_molecule_measurement:{own_median:8,reference_median:7,median_difference:1,own_n_records:3},summary:{neighborhood:{},n_records:0},records:[],neighbors:[],sources:[],svg:'',neighborhood_svg:''};
    await context.inspect(null,measuredDetail);assert.match($('query-match-note').textContent,/你的测定：中位 8（3 条记录）；参照分子中位 7；差值 \+1/);
  },
  async audit(){
    const {context,$,api}=await setup();const report={query_molecules:{input_rows:5,mapped:1,invalid_smiles:1,duplicate_query:0,in_reference:1,neighbor_evidence_counts:{active:1},reference_support_assessed:false,excluded_reasons:{non_exact_relation:1},document_pair_offsets:[{document_a:'R1',document_b:'R3',shared_entries:3,median_offset:1.2,sign_agreement:1,systematic_offset:true}]},query_reference_matches:[{query_id:'CTRL',molecule_id:'R1',reference_id:'R1',source_count:1,linked_source_count:1,own_median:8,reference_median:7,median_difference:1}],warnings:[]};
    for(const language of ['en','zh']){api.setLanguage(language);context.showDataAudit(report);const audit=$('audit-body');assert.ok(audit.textContent.includes(api.t('u8.supportNotAssessed')));assert.ok(audit.textContent.includes(api.t('u8.ownMedian')));assert.match(audit.textContent,/R1.*R3/);assert.ok(audit.textContent.includes(api.t('u8.differenceSummary',{median:'+1',minimum:'+1',maximum:'+1'})));assert.equal(audit.querySelectorAll('.new-offset-pairs tbody tr').length,1);}
  },
  async outlines(){
    const {context,$,node}=await setup();
    const svg=node('svg'), translated=node('g'),cell=node('g'),rect=node('rect');
    translated.setAttribute('transform','translate(100,200)');cell.dataset.gridRow='0';cell.dataset.gridCol='0';cell.classList.add('reference-evidence-faded');rect.setAttribute('x','12');rect.setAttribute('fill','green');
    rect.cloneNode=()=>{const copy=node('rect');for(const [key,value] of Object.entries(rect.attributes))copy.setAttribute(key,value);return copy;};
    cell.append(rect);translated.append(cell);svg.append(translated);$('map-stage').append(svg);
    cell.querySelector=selector=>selector===':scope > rect'?rect:null;
    Object.assign(context,{mapCells:[cell],points:[reference],currentCell:{row:0,col:0},neighborRelations:new Map(),highlightedGroup:'',transform:{x:0,y:0}});
    $('measurement-filter').value='all';vm.runInContext(realFunction('highlight'),context);context.highlight();
    const outlines=svg.children.find(item=>item.dataset.mapOutlines);assert.ok(outlines);assert.equal(outlines.parentNode,svg);assert.equal(outlines.children[0].getAttribute('transform'),'translate(100,200)');
    cell.getBoundingClientRect=()=>({left:20,top:30,width:10,height:10});$('map-viewport').getBoundingClientRect=()=>({left:0,top:0,width:100,height:100});context.centerOnCurrentCell();assert.equal(context.transform.x,25);assert.equal(context.transform.y,15);
  },
  async payload(){
    const {context,$,api}=await setup(), payload=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
    Object.assign(context,{job:payload.job,newMolecules:payload.new,points:payload.points});
    for(const language of ['en','zh']){
      api.setLanguage(language);context.renderNewList(true);
      const row=$('new-list-table').querySelectorAll('tbody tr').find(row=>row.dataset.queryId==='N1');
      assert.ok(Number(row.children[5].textContent)>0, '真实 Tanimoto 字段必须显示');
      context.showQueryInput(payload.new.input,payload.new.options);assert.equal($('new-activity-col').value,'value');assert.equal($('new-document-col').value,'batch');
      context.showQuerySelection(payload.detail);assert.equal($('new-records').querySelectorAll('tbody tr').length,2);assert.match($('new-records').textContent,/nM.*uM/);
      assert.ok($('query-summary').textContent.includes(api.t('u8.neighborCountsUnassessed',{similar:payload.detail.molecule.similar_reference_count,classes:payload.detail.molecule.evidence_reference_classes.split(';').filter(Boolean).map(value=>api.t('d.'+value)).join(' / ') || '—'})));
      context.showDataAudit(payload.job.report);assert.ok($('audit-body').textContent.includes(api.t('u8.ownMedian')));assert.equal($('audit-body').querySelectorAll('.new-excluded-records tbody tr').length,1);
    }
  },
  async u8b_formatting(){
    const {context,$,api}=await setup();
    const record={compound:'LAB-1',run:'B1',amount:30,units:'nM',rel:'>',exclusion_reason:'non_exact_relation',private_note:'不应显示'};
    const report={query_molecules:{mapped:1,excluded_reasons:{non_exact_relation:57,invalid_activity_value:0},
      excluded_records:Array.from({length:57},(_,i)=>({...record,compound:`LAB-${i+1}`})),
      excluded_record_columns:{compound_id:'compound',evidence_unit:'run',value:'amount',unit:'units',relation:'rel'},
      no_retained_activity_compounds:['LAB-1','LAB-2'],document_pair_offsets:[
        {document_a:'B1',document_b:'B2',shared_entries:13,median_offset:-.025295,sign_agreement:.615385,systematic_offset:false},
        {document_a:'B1',document_b:'B3',shared_entries:5,median_offset:1.175469,sign_agreement:1,systematic_offset:true}]},
      query_reference_matches:[-.08240689,.038395,.1640176].map((value,i)=>({query_id:`CTRL-${i}`,molecule_id:`R${i}`,reference_id:`R${i}`,own_median:7+value,reference_median:7,median_difference:value})),warnings:[]};
    vm.runInContext(realFunction('inspect'),context);
    for(const language of ['en','zh']){
      api.setLanguage(language);context.showDataAudit(report);
      const audit=$('audit-body');
      assert.equal(audit.querySelectorAll('.new-excluded-reasons tbody tr').length,1);
      const table=audit.querySelector('.new-excluded-records');
      assert.equal(table.querySelectorAll('th').length,6);
      assert.equal(table.querySelectorAll('tbody tr').length,50);
      assert.equal(table.querySelectorAll('th').map(n=>n.textContent).join('|'),['d.Compound ID','d.Evidence units','u8.value','u8.unit','u8.relation','u8.reason'].map(api.t).join('|'));
      assert.equal(table.querySelector('tbody tr').children.slice(0,5).map(n=>n.textContent).join('|'),'LAB-1|B1|30|nM|>');
      assert.ok(!table.textContent.includes('不应显示'));
      assert.ok(audit.textContent.includes(api.t('u8.excludedMore',{count:7})));
      assert.ok(audit.textContent.includes('new_molecules_excluded_records.csv'));
      const matches=audit.querySelectorAll('.query-reference-matches tbody tr');
      assert.equal(matches.map(n=>n.children.at(-1).textContent).join('|'),'−0.082|+0.038|+0.164');
      assert.ok(audit.textContent.includes(language==='zh'?'−0.082 至 +0.164':'−0.082 to +0.164'));
      assert.ok(audit.textContent.includes(language==='zh'?'LAB-1、LAB-2':'LAB-1, LAB-2'));
      const offsets=audit.querySelectorAll('.new-offset-pairs tbody tr');
      assert.equal(offsets[0].children[2].textContent,'−0.03');assert.equal(offsets[0].children[3].textContent,'62%');
      assert.equal(offsets[1].children[2].textContent,'+1.18');assert.equal(offsets[1].children[3].textContent,'100%');
      context.showQuerySelection({molecule:{...added,structure_group:'',placement_distance:6.08276253},query_neighbors:[]});
      assert.equal($('new-placement').textContent,api.t('u8.placement',{id:'R1',distance:'6.1'}));
      assert.ok($('query-summary').textContent.includes(api.t('u8.noStructureGroup',{similarity:'0.6'})));
      assert.equal($('query-summary').querySelectorAll('dt').filter(n=>n.textContent===api.t('u8.neighborEvidence')).length,1);
      assert.equal($('query-summary').querySelectorAll('dt').filter(n=>n.textContent===api.t('u8.similarReferenceMolecules')).length,1);
      assert.equal($('new-map-markers').hidden,false);assert.equal($('new-map-markers').textContent,api.t('u8.mapMarkers'));
      for(const difference of [-.08240689,.1640176]){
        await context.inspect(null,{molecule:reference,summary:{neighborhood:{},n_records:0},records:[],neighbors:[],sources:[],svg:'',neighborhood_svg:'',query_match:{query_id:'CTRL',status:'in_reference',reference_id:'R1'},new_molecule_measurement:{own_median:7+difference,reference_median:7,median_difference:difference,own_n_records:2}});
        assert.ok($('query-match-note').textContent.includes(difference<0?'−0.082':'+0.164'));
      }
    }
    report.query_molecules.excluded_records=[{compound:'LAB',amount:0,exclusion_reason:'nonpositive_activity_value'}];
    report.query_molecules.excluded_record_columns={compound_id:'compound',value:'amount'};
    report.query_molecules.excluded_record_defaults={unit:'p',relation:'='};
    context.showDataAudit(report);const row=$('audit-body').querySelector('.new-excluded-records tbody tr');
    assert.equal(row.children.slice(0,5).map(n=>n.textContent).join('|'),'LAB|—|0|p|=');
    report.query_molecules.excluded_record_columns.unit='units';
    report.query_molecules.excluded_record_columns.relation='rel';
    report.query_molecules.excluded_records[0].units=null;
    report.query_molecules.excluded_records[0].rel=null;
    context.showDataAudit(report);
    const missing=$('audit-body').querySelector('.new-excluded-records tbody tr');
    assert.equal(missing.children.slice(3,5).map(n=>n.textContent).join('|'),'—|—','已选择列中的空值不能显示为默认单位或精确关系');
  },
  async u8c_match_notes(){
    const {context,$,api}=await setup();
    vm.runInContext(realFunction('inspect'),context);
    const detail={molecule:reference,summary:{neighborhood:{},n_records:0},records:[],neighbors:[],sources:[],svg:'',neighborhood_svg:''};
    const duplicate={query_id:'DUPREF_001',first_query_id:'K0000',status:'duplicate_query',reference_id:'R1'};
    const measurement={own_median:8,reference_median:7,median_difference:1,own_n_records:3};
    const expected={
      en:[
        'DUPREF_001 repeats K0000 (same structure). K0000 is not shown separately: it has the same structure as the reference molecule R1. Its evidence and sources are shown below.',
        'K0000 is not shown separately: it has the same structure as the reference molecule R1. Its evidence and sources are shown below. Your measurements: median 8 (3 records); reference molecule median 7; difference +1.',
        'DUPREF_001 repeats K0000 (same structure). K0000 is not shown separately: it has the same structure as the reference molecule R1. Its evidence and sources are shown below. Your measurements: median 8 (3 records); reference molecule median 7; difference +1.'
      ],
      zh:[
        'DUPREF_001 与 K0000 结构相同（重复）。K0000 未单独入图：它与参照分子 R1 结构相同。下面显示该参照分子的证据与来源文献。',
        'K0000 未单独入图：它与参照分子 R1 结构相同。下面显示该参照分子的证据与来源文献。你的测定：中位 8（3 条记录）；参照分子中位 7；差值 +1。',
        'DUPREF_001 与 K0000 结构相同（重复）。K0000 未单独入图：它与参照分子 R1 结构相同。下面显示该参照分子的证据与来源文献。你的测定：中位 8（3 条记录）；参照分子中位 7；差值 +1。'
      ]
    };
    for(const language of ['en','zh']){
      api.setLanguage(language);
      const cases=[
        {query_match:duplicate},
        {query_match:{query_id:'K0000',status:'in_reference',reference_id:'R1'},new_molecule_measurement:measurement},
        {query_match:duplicate,new_molecule_measurement:measurement}
      ];
      for(const [index,match] of cases.entries()){
        await context.inspect(null,{...detail,...match});
        assert.equal($('query-match-note').hidden,false);
        assert.equal($('query-match-note').textContent,expected[language][index]);
      }
    }
  },
  async u8c_neighbor_counts(){
    const {context,$,api}=await setup();
    const cases=[
      {similar:0,supported:0,en:'0 similar reference molecules; 0 with supported classes.',zh:'相近的参照分子 0 个，其中类别可靠的 0 个。'},
      {similar:1,supported:0,classes:'',en:'1 similar reference molecule; 0 with supported classes.',zh:'相近的参照分子 1 个，其中类别可靠的 0 个。'},
      {similar:3,supported:0,en:'3 similar reference molecules; 0 with supported classes.',zh:'相近的参照分子 3 个，其中类别可靠的 0 个。'},
      {similar:1,supported:1,classes:'active',en:'1 similar reference molecule; 1 with supported classes: active.',zh:'相近的参照分子 1 个，其中类别可靠的 1 个：活性。'},
      {similar:3,supported:2,en:'3 similar reference molecules; 2 with supported classes: active / inactive.',zh:'相近的参照分子 3 个，其中类别可靠的 2 个：活性 / 非活性。'},
      {similar:1,supported:0,assessed:false,en:'Similar reference molecules: 1; observed classes: medium.',zh:'相近参照分子 1 个；实测类别 中等活性。'}
    ];
    for(const language of ['en','zh']){
      api.setLanguage(language);
      for(const item of cases){
        context.showQuerySelection({molecule:{...added,similar_reference_count:item.similar,supported_reference_count:item.supported,
          supported_reference_classes:item.classes ?? 'active;inactive',evidence_reference_classes:'medium',reference_support_assessed:item.assessed ?? true},query_neighbors:[]});
        const stats=$('query-summary').querySelector('dl');
        const index=stats.children.findIndex(n=>n.tagName==='DT' && n.textContent===(language==='en'?'Similar reference molecules':'相近的参照分子'));
        assert.ok(index>=0);
        assert.equal(stats.children[index+1].textContent,item[language]);
        assert.equal($('new-support-note').hidden,item.assessed !== false);
      }
    }
  },
  async u8b_markers(){
    const {context,$,node}=await setup();
    const refs=Array.from({length:7},(_,i)=>({...reference,molecule_id:`R${i+1}`,molecule_identity_key:`R${i+1}`,grid_col:i}));
    const selected={...added,grid_row:1,grid_col:1};
    const second={...selected,molecule_id:'N3',query_id:'N3',grid_col:2,nearest_reference_id:'R7'};
    const svg=node('svg'),wrap=node('g');wrap.setAttribute('transform','translate(100,200)');svg.append(wrap);$('map-stage').append(svg);
    const stageQuery=$('map-stage').querySelector.bind($('map-stage'));
    $('map-stage').querySelector=selector=>selector==='[data-map-outlines]' ? svg.children.find(n=>n.dataset.mapOutlines) || null : stageQuery(selector);
    const matrix=(x=0,y=0)=>({x,y,inverse(){return matrix(-x,-y);},multiply(other){return matrix(x+other.x,y+other.y);}});
    svg.getCTM=()=>matrix(10,20);
    svg.createSVGPoint=()=>({x:0,y:0,matrixTransform(m){return {x:this.x+m.x,y:this.y+m.y};}});
    const cells=[...refs,selected,second].map(p=>{
      const cell=node('g'),rect=node('rect');cell.dataset.gridRow=p.grid_row;cell.dataset.gridCol=p.grid_col;
      rect.setAttribute('x',p.grid_col*100);rect.setAttribute('y',p.grid_row*100);rect.setAttribute('width',100);rect.setAttribute('height',100);
      rect.getBBox=()=>({x:p.grid_col*100,y:p.grid_row*100,width:100,height:100});
      rect.getCTM=()=>matrix(p.is_new?10:110,p.is_new?20:220);
      rect.cloneNode=()=>{const copy=node('rect');for(const [key,value] of Object.entries(rect.attributes))copy.setAttribute(key,value);return copy;};
      cell.append(rect);cell.querySelector=selector=>selector===':scope > rect'?rect:null;
      (p.is_new?svg:wrap).append(cell);if(!p.is_new)cell.classList.add('reference-evidence-faded');return cell;
    });
    Object.assign(context,{points:[...refs,selected,second],mapCells:cells,neighborRelations:new Map(),transform:{s:1,x:0,y:0},tileWidth:100,hasStructures:true});
    $('measurement-filter').value='all';$('show-neighbor-links').checked=false;
    vm.runInContext(['highlight','inspect','closeInspector','applyTransform'].map(realFunction).join('\n'),context);
    // 第五名的两项并列使真实近邻表有六项；第七个分子不应标出。
    const detail={molecule:selected,summary:{},query_neighbors:refs.slice(0,6).map((p,i)=>({measured_id:p.molecule_id,rank:Math.min(i+1,5),tanimoto:i<4?.9-i*.05:.6})),records:[],svg:''};
    await context.inspect(null,detail);
    assert.equal(cells.filter(n=>n.classList.contains('neighbor-representation_only')).length,6);
    assert.ok(!cells[6].classList.contains('neighbor-representation_only'));
    let outlines=svg.children.find(n=>n.dataset.mapOutlines),line=outlines.children.find(n=>n.dataset.newReferenceLink);
    assert.equal(outlines.getAttribute('pointer-events'),'none');assert.ok(line);
    assert.equal(outlines.children.filter(n=>n.dataset.outline==='representation_only').length,6);
    assert.equal(['x1','y1','x2','y2'].map(k=>line.getAttribute(k)).join(','),'150,150,150,250');
    const css=fs.readFileSync(path.join(process.argv[2],'app.css'),'utf8');
    assert.match(css,/\[data-new-reference-link\][^}]*stroke: #173e34[^}]*stroke-width: calc\(1\.5 \* var\(--zoom-inverse[^}]*stroke-dasharray:/);
    for(const scale of [.25,1,4]){context.transform.s=scale;context.applyTransform();assert.equal(1.5*Number($('map-stage').style['--zoom-inverse'])*scale,1.5);}
    await context.inspect(null,{...detail,molecule:second,query_neighbors:[{measured_id:'R7',rank:1,tanimoto:.9}]});
    assert.equal(cells.filter(n=>n.classList.contains('neighbor-representation_only')).length,1);
    outlines=svg.children.find(n=>n.dataset.mapOutlines);line=outlines.children.find(n=>n.dataset.newReferenceLink);
    assert.equal(['x1','y1','x2','y2'].map(k=>line.getAttribute(k)).join(','),'250,150,750,250');
    context.closeInspector();assert.equal(svg.children.filter(n=>n.dataset.mapOutlines).length,0);
    assert.equal(cells.filter(n=>n.classList.contains('neighbor-representation_only')).length,0);
    await context.inspect(null,detail);
    await context.inspect(null,{molecule:refs[0],summary:{neighborhood:{},n_records:0},neighbors:[],sources:[],records:[],svg:'',neighborhood_svg:''});
    assert.equal(svg.children.flatMap(n=>n.children).filter(n=>n.dataset.newReferenceLink).length,0);
    assert.equal(cells.filter(n=>n.classList.contains('neighbor-representation_only')).length,0);
  },
  async u8b_summary_exports(){
    const {context,$,api,requests}=await setup();
    vm.runInContext(realFunction('stat'),context);
    const fn=realFunction('showResult'),start=fn.indexOf('  const report = job.report;'),end=fn.indexOf("  $('summary').hidden = false;");
    vm.runInContext(`function showSummary(){${fn.slice(start,end)}}`,context);
    for(const language of ['en','zh']){
      api.setLanguage(language);context.job.molecules=934;context.job.report={query_molecules:{mapped:45}};context.showSummary();
      assert.equal($('summary').children[1].children[0].textContent,'889');assert.equal($('summary').children[1].children[1].textContent,api.t('u8.referenceMolecules'));
      assert.equal($('summary').children[2].children[0].textContent,'45');
      context.job.molecules=889;context.job.report={};context.showSummary();
      assert.equal($('summary').children.length,4);assert.equal($('summary').children[1].children[1].textContent,api.t('d.Mapped molecules'));
    }
    const begin=fn.indexOf("  $('downloads').replaceChildren();"),finish=fn.indexOf("  $('saved-location').textContent",begin);
    vm.runInContext(`function showExports(){${fn.slice(begin,finish)}}`,context);
    context.newMolecules.files=['new_molecules_excluded_records.csv'];context.download=(url,name)=>requests.push({url,name});context.showExports();
    const button=$('downloads').children.find(n=>n.textContent==='new_molecules_excluded_records.csv');assert.ok(button);button.click();
    assert.equal(requests.at(-1).url,'/api/file?name=new_molecules_excluded_records.csv');
  },
  async projection(){const {context}=await setup();assert.equal(context.projectionSVG().children.length,1);},
};
scenarios[process.argv[3]]().catch(error=>{console.error(error);process.exitCode=1;});
