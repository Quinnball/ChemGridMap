'use strict';
const assert=require('node:assert/strict'),vm=require('node:vm'),fs=require('node:fs'),path=require('node:path');
const {runtime,realFunction,dictionaries}=require('./u4_web_contract.cjs');
const reference={molecule_id:'R1',grid_row:0,grid_col:0,activity_class:'active'};
const measured={query_id:'M001',molecule_id:'M001',is_new:true,status:'mapped',canonical_smiles:'CCF',grid_row:0,grid_col:1,measurement_status:'measured',activity_class:'inactive',measured_class:'inactive',activity_pchembl:5,measured_value:10000,measured_unit:'nM',class_source:'value',result_source:'manual',neighbor_evidence:'active',neighbor_comparison:'disagrees',label:'A',batch:'B1',note:'my note',recorded_at:'2026-11-12T10:00:00+00:00'};
async function setup(){
 const env=await runtime(),requests=[],downloads=[],$=id=>env.document.getElementById(id);
 const job={state:'done',options:{representation:'morgan',input_format:'molecule',lower_threshold:6,upper_threshold:7},files:['map.svg'],output_dir:'/local/run',molecules:2,report:{query_molecules:{mapped:1,feedback_counts:{recorded:1,agrees:0,disagrees:1,no_comparison:0}}}};
 Object.assign(env.context,{$,t:env.api.t,translateServer:env.api.translateServer,busy:false,job,points:[reference,measured],view:'grid',newMolecules:{rows:[measured],directory:'/local/run/new_molecules/manual',files:['new_molecules/manual/neighbor_evidence_feedback.csv'],replacement:{total:1,manual:1,recorded:0,directory:'/local/run/new_molecules/manual'}},newListSort:{key:'query_id',direction:1},queryInputInfo:{},inputInfo:null,structureGroups:[],highlightedGroup:'',mapCells:[],currentCell:null,
 colors:{active:'#2a9d8f',medium:'#e6a23c',inactive:'#c44e52',unlabelled:'#b8c2cc'},fadedColors:{},reliabilityBins:[[0,'#cde1f5']],
 setBusy(value){env.context.busy=value;},showResult:async()=>{},updateNavigation(){},closeInspector(){},highlight(){},translateSVGText(){},applyTransform(){},setView(value){env.context.view=value;},centerOnCurrentCell(){},
 inspect:async(query)=>requests.push({inspect:query}),download:async(url,name)=>downloads.push({url,name}),
 getJSON:async(url,options)=>{requests.push({url,options});return {status:'ok',canonical_smiles:'CCF',svg:'<svg/>',query_id:'M001',row:measured,point:measured,job};}
 });
 env.window.confirm=()=>true;
 const names=['bindNewMolecules','downloadNewList','table','measuredLabel','showQuerySelection','refreshValidation','updateNewNavigation','newMoleculeOptions','insertNewMolecules','newListColumns','visibleNewMolecules','newListValue','renderNewList','format','countLabel','percentage','referenceMapDetails','renderReferenceMap','showDataAudit','newMoleculeWarning','legendItem','renderColorLegend','renderLegend','colorMap','reliabilityColor','status','fail'];
 vm.runInContext(names.map(realFunction).join('\n'),env.context);
 env.context.bindNewMolecules();
 return {...env,$,requests,downloads,job};
}
const scenarios={
 async payload(){
  const {context,$,api}=await setup(),payload=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
  Object.assign(context,{job:payload.job,newMolecules:payload.new,points:payload.points});
  for(const language of ['en','zh']){
   api.setLanguage(language);context.renderNewList(true);
   assert.equal($('new-list-feedback-counts').textContent,api.t('u9b.counts',payload.new.report.feedback_counts));
   for(const detail of payload.details){
    context.renderMeasuredResult(detail);const result=detail.measured_result;
    assert.equal($('measured-result-actions').hidden,!result.editable);
    assert.ok($('measured-result-summary').textContent.includes(api.t('u9b.source.'+result.result_source)));
    if(result.result_source==='manual'){$('edit-result').click();assert.equal($('result-unit').value,'nM');assert.equal($('result-value').value,'15');}
    if(result.result_source==='recorded'){$('edit-result').click();assert.equal($('result-mode').value,'class');assert.equal($('result-class').value,'inactive');assert.equal($('result-value').value,'');}
   }
  }
 },
 async errors(){
  const {api}=await setup();
  const appSource=fs.readFileSync(path.join(process.argv[2],'..','app.py'),'utf8');
  const source=appSource.slice(appSource.indexOf('def normalize_measurement('),appSource.indexOf('\nclass AppSession('));
  assert.match(source,/class ManualMoleculesMixin/);
  for(const [,message] of source.matchAll(/raise ValueError\("([^"]+)"\)/g)){api.setLanguage('en');assert.equal(api.translateServer(message),message);api.setLanguage('zh');assert.notEqual(api.translateServer(message),message,message);}
 },
 async controls(){
  const {context,$,api,document}=await setup();
  assert.equal(document.querySelectorAll('[data-help]').length,77);
  for(const language of ['en','zh']){
   api.setLanguage(language);
   for(const tab of ['upload','one']){$('new-'+tab+'-tab').click();assert.equal($('new-'+tab+'-panel').hidden,false);assert.equal($('new-'+tab+'-tab').getAttribute('aria-selected'),'true');assert.equal($('new-'+(tab==='one'?'upload':'one')+'-panel').hidden,true);}
   assert.equal($('new-one-tab').textContent,api.t('u9b.oneTab'));
   for(const prefix of ['single','result'])for(const kind of ['mode','value','unit','class','batch','note']){
    const help=document.querySelector(`[data-help-control=${prefix+'-'+kind}]`);assert.ok(help);help.click();assert.equal($('parameter-tooltip').textContent,api.t(help.dataset.help));help.click();
   }
   const help=document.querySelector('[data-help-control=measured-result]');assert.match(api.t(help.dataset.help),language==='zh'?/不是预测.*活性悬崖.*不会改变参照图.*重启软件/:/not a prediction.*activity cliffs.*does not change.*restarting/);
  }
 },
 async preview_cache(){
  const {context,$,api,requests}=await setup();
  for(const language of ['en','zh']){
   api.setLanguage(language);$('single-smiles').value='CCF';
   for(const old of ['duplicate_new','in_reference']){
    context.previewSingleMolecule.last={smiles:'CCF',status:old,reference_id:'R1',query_id:'M001'};
    context.setNewMoleculeTab('upload');context.setNewMoleculeTab('one');await new Promise(resolve=>setImmediate(resolve));
    assert.equal(requests.at(-1).url,'/api/new-molecules/preview');assert.equal(context.previewSingleMolecule.last.status,'ok');
    context.previewSingleMolecule.last={smiles:'CCF',status:old};$('open-new-molecules').click();await new Promise(resolve=>setImmediate(resolve));
    assert.equal(context.previewSingleMolecule.last.status,'ok');assert.equal($('single-preview-status').textContent,api.t('u9b.previewOk',{smiles:'CCF'}));
   }
  }
 },
 async preview(){
  const {context,$,api,requests}=await setup();
  for(const language of ['en','zh']){
   api.setLanguage(language);
   for(const [status,key,id] of [['ok','previewOk',''],['invalid','previewInvalid',''],['in_reference','previewReference','R1'],['duplicate_new','previewDuplicate','M001']]){
    const result={status,canonical_smiles:'CCF',svg:status==='invalid'?'':'<svg/>',...(status==='in_reference'?{reference_id:id}:{query_id:id})};
    context.getJSON=async(url,options)=>{requests.push({url,options});return result;};$('single-smiles').value='CCF';await context.previewSingleMolecule();
    assert.equal($('single-preview-status').textContent,api.t('u9b.'+key,{smiles:'CCF',id}));assert.equal($('single-preview-view').hidden,!id);
    if(id){$('new-molecules-dialog').showModal();$('single-preview-view').click();assert.equal(requests.at(-1).inspect,status==='in_reference'?'row=0&col=0':'id='+id);assert.equal($('new-molecules-dialog').open,false);}
   }
  }
  context.points=[{...measured,query_id:'R1',molecule_id:'R1',grid_row:3,grid_col:4},reference];
  context.renderSinglePreview({status:'in_reference',reference_id:'R1'});$('single-preview-view').click();assert.equal(requests.at(-1).inspect,'row=0&col=0');
  context.renderSinglePreview({status:'duplicate_new',query_id:'R1'});$('single-preview-view').click();assert.equal(requests.at(-1).inspect,'id=R1');
  let delay,callback,cancelled=0;context.setTimeout=(fn,ms)=>{delay=ms;callback=fn;return 9;};context.clearTimeout=()=>cancelled++;
  $('single-smiles').value='CCCl';$('single-smiles').dispatchEvent({type:'input'});assert.equal(delay,400);assert.equal(cancelled,1);assert.equal(context.previewSingleMolecule.last,null);assert.equal(callback,context.previewSingleMolecule);
  // 一个旧请求晚返回时，不得覆盖后来输入的预览。
  let resolve;context.getJSON=()=>new Promise(done=>{resolve=done;});const pending=context.previewSingleMolecule();$('single-smiles').value='CCBr';context.scheduleSinglePreview();resolve({status:'ok',canonical_smiles:'CCCl',svg:'OLD'});await pending;assert.notEqual($('single-preview-svg').innerHTML,'OLD');
 },
 async validation(){
  const {context,$,api,requests}=await setup();
  for(const language of ['en','zh']){
   api.setLanguage(language);
   for(const [form,id] of [['single-molecule-form','single-smiles'],['single-molecule-form','single-value'],['measured-result-form','result-value']]){
    $(id).validity={valueMissing:true,valid:false};$(id).disabled=false;$(form).dispatchEvent({type:'invalid',target:$(id)});assert.equal($(id).validationMessage,language==='zh'?api.t('d.Please fill out this field.'):'');
    $(id).validity={valid:true};$(form).dispatchEvent({type:'input',target:$(id)});assert.equal($(id).validationMessage,'');
   }
   for(const prefix of ['single','result']){
    $(prefix+'-mode').value='value';context.updateMeasurementForm(prefix);assert.equal($(prefix+'-numeric').hidden,false);assert.equal($(prefix+'-qualitative').hidden,true);
    for(const value of ['0','-1','NaN','Infinity','> 10','0x10','abc','']){$(prefix+'-value').value=value;assert.throws(()=>context.readMeasurement(prefix),{message:api.t('u9b.invalidValue')});}
    $(prefix+'-value').value='12.5';$(prefix+'-unit').value='nM';$(prefix+'-batch').value='batch';$(prefix+'-note').value='note';assert.equal(context.readMeasurement(prefix).value,12.5);
    $(prefix+'-unit').value='watts';assert.throws(()=>context.readMeasurement(prefix),{message:api.t('u9b.invalidUnit')});
    $(prefix+'-mode').value='class';$(prefix+'-class').value='medium';context.updateMeasurementForm(prefix);const result=context.readMeasurement(prefix);assert.equal(result.measured_class,'medium');assert.equal(result.value,undefined);assert.equal($(prefix+'-value').disabled,true);
   }
   $('single-mode').value='none';assert.equal(context.readMeasurement('single'),null);
   $('single-mode').value='value';$('single-value').value='>10';const count=requests.length;await context.addSingleMolecule({preventDefault(){}});assert.equal(requests.length,count);assert.equal($('single-feedback').textContent,api.t('u9b.invalidValue'));
  }
 },
 async add(){
  const {context,$,api,requests}=await setup();let refreshed=0;
  context.refreshNewMoleculeMutation=async()=>{refreshed++;};
  for(const language of ['en','zh']){
   api.setLanguage(language);$('new-molecules-dialog').showModal();$('single-smiles').value='CCF';$('single-id').value='';$('single-label').value='A';$('single-mode').value='value';$('single-value').value='10000';$('single-unit').value='nM';$('single-note').value='my note';
   context.previewSingleMolecule.last={smiles:'CCF',status:'ok'};await context.addSingleMolecule({preventDefault(){}});
   const request=requests.findLast(r=>r.url==='/api/new-molecules/add');const body=JSON.parse(request.options.body);assert.equal(body.measurement.value,10000);assert.equal(body.measurement.unit,'nM');assert.equal(body.query_id,'');assert.equal(body.note,'my note');
   assert.equal($('new-molecules-dialog').open,true);assert.equal($('single-smiles').value,'');assert.equal($('single-id').value,'');assert.equal($('single-label').value,'');assert.equal($('single-mode').value,'none');assert.equal($('single-feedback').textContent,api.t('u9b.added',{id:'M001',category:api.t('u8.active')}));
  }
  assert.equal(refreshed,2);
 },
 async measured(){
  const {context,$,api,requests}=await setup();context.refreshNewMoleculeMutation=async()=>{};
  for(const language of ['en','zh']){
   api.setLanguage(language);
   context.renderMeasuredResult({molecule:{...measured,measurement_status:'not_measured',result_source:''}});assert.equal($('record-result').hidden,false);$('record-result').click();assert.equal($('measured-result-form').hidden,false);
   for(const source of ['manual','recorded','upload']){
    context.renderMeasuredResult({molecule:{...measured,result_source:source},measured_result:{editable:source!=='upload'}});
    assert.equal($('record-result').hidden,true);assert.equal($('measured-result-actions').hidden,source==='upload');assert.ok($('measured-result-summary').textContent.includes(api.t('u9b.disagrees')));assert.match($('measured-result-summary').textContent,/B1.*my note.*2026-11-12/);
    if(source==='upload'){assert.ok($('measured-result-summary').textContent.includes(api.t('u9b.fromUpload')));assert.equal($('measured-result-form').hidden,true);}
    else {$('edit-result').click();assert.equal($('result-value').value,'10000');assert.equal($('result-unit').value,'nM');}
   }
   context.showQuerySelection({molecule:{...measured,result_source:'recorded',neighbor_evidence:'mixed',neighbor_comparison:'no_comparison'},query_neighbors:[],records:[]});assert.ok($('query-summary').textContent.includes(api.t('u9b.no_comparison')));assert.ok(!$('query-summary').textContent.includes('u8.no_comparison'));assert.ok($('measured-result-summary').textContent.includes(api.t('u9b.noComparisonReason',{category:api.t('u8.mixed')})));
   $('edit-result').click();$('result-mode').value='class';$('result-class').value='active';await context.saveMeasuredResult({preventDefault(){}});assert.equal(requests.findLast(r=>r.url)?.url,'/api/new-molecules/result');assert.equal(JSON.parse(requests.findLast(r=>r.url).options.body).measurement.measured_class,'active');
   await context.saveMeasuredResult(null,true);assert.equal(requests.findLast(r=>r.url)?.url,'/api/new-molecules/result/clear');assert.equal(requests.at(-1).inspect,'id=M001');
  }
 },
 async list_audit(){
  const {context,$,api,downloads}=await setup();
  context.newMolecules.rows=[measured,{...measured,query_id:'M002',neighbor_comparison:'agrees',result_source:'recorded'},{...measured,query_id:'M003',neighbor_comparison:'no_comparison',neighbor_evidence:'mixed',result_source:'upload'},{...measured,query_id:'M004',measurement_status:'not_measured',neighbor_comparison:'',result_source:''}];
  for(const language of ['en','zh']){
   api.setLanguage(language);$('new-mismatch-filter').checked=false;context.renderNewList(true);assert.equal($('new-list-table').querySelectorAll('tbody tr').length,4);
   assert.equal($('new-list-feedback-counts').textContent,api.t('u9b.counts',{recorded:3,agrees:1,disagrees:1,no_comparison:1}));assert.ok($('new-list-table').textContent.includes(api.t('u9b.resultSource')));assert.match($('new-list-table').textContent,/my note/);
   $('new-mismatch-filter').checked=true;$('new-mismatch-filter').dispatchEvent({type:'change'});assert.equal($('new-list-table').querySelectorAll('tbody tr').length,1);assert.equal($('new-list-table').querySelector('tbody tr').dataset.queryId,'M001');
   context.showDataAudit(context.job.report);assert.ok($('audit-body').textContent.includes('neighbor_evidence_feedback.csv'));assert.ok($('audit-body').textContent.includes(api.t('u9b.counts',{recorded:1,agrees:0,disagrees:1,no_comparison:0})));
   context.newMolecules.files=['new_molecules/manual/neighbor_evidence_feedback.csv'];$('download-neighbor-feedback').click();assert.match(downloads.at(-1).url,/neighbor_evidence_feedback.csv/);
   context.newMolecules.files=[];context.renderNewList();assert.equal($('download-neighbor-feedback').disabled,false);$('download-neighbor-feedback').click();assert.equal(downloads.at(-1).url,'/api/new-molecules/feedback');
  }
 },
 async replacement(){
  const {context,$,api,window,requests}=await setup();
  for(const language of ['en','zh']){
   api.setLanguage(language);let seen;window.confirm=text=>{seen=text;return false;};await context.insertNewMolecules({preventDefault(){}});assert.equal(requests.length,0);assert.equal(seen,api.t('u9b.replaceConfirm',context.newMolecules.replacement));
   window.confirm=()=>true;assert.equal(context.confirmNewMoleculeReplacement(),true);
   const original=context.newMolecules.replacement;context.newMolecules.replacement={total:10,manual:0,recorded:0};window.confirm=()=>{throw Error('ordinary upload should not ask');};assert.equal(context.confirmNewMoleculeReplacement(),true);context.newMolecules.replacement=original;
  }
 },
 async markers(){
  const {context,$,api,node}=await setup();
  const cell=node('g'),rect=node('rect');rect.setAttribute('width','100');rect.setAttribute('height','100');rect.setAttribute('x','100');rect.setAttribute('y','200');rect.setAttribute('fill','#c44e52');cell.dataset.gridRow=0;cell.dataset.gridCol=1;cell.append(rect);
  cell.querySelector=selector=>selector===':scope > rect'?rect:cell.children.find(item=>selector.includes('data-neighbor-mismatch')&&item.dataset.neighborMismatch)||null;
  $('map-stage').querySelectorAll=()=>[cell];
  for(const language of ['en','zh'])for(const mode of ['activity','support','reliability','structure','neighbor_evidence']){
   api.setLanguage(language);$('color-mode').value=mode;context.colorMap();const markers=cell.children.filter(item=>item.dataset.neighborMismatch);assert.equal(markers.length,1);assert.match(markers[0].getAttribute('d'),/^M111,202/);assert.ok($('legend').textContent.includes(api.t('u9b.mismatchLegend')));
  }
  const backendMarker=cell.children.find(item=>item.dataset.neighborMismatch);backendMarker.setAttribute('data-rendered-by','server');context.colorMap();assert.equal(cell.children.filter(item=>item.dataset.neighborMismatch).length,1);assert.equal(cell.children.find(item=>item.dataset.neighborMismatch),backendMarker);assert.equal(backendMarker.getAttribute('fill'),'#6f1d6b');assert.equal($('legend').querySelector('.neighbor-mismatch-legend').firstElementChild.style.background,'#6f1d6b');
  context.points=[reference,{...measured,neighbor_comparison:'agrees'}];context.colorMap();assert.equal(cell.children.filter(item=>item.dataset.neighborMismatch).length,0);
 },
 async mutation_highlight(){
  const {context,$,job}=await setup();let cleared=0;
  context.clearSourceHighlight=()=>{cleared++;context.highlight.source=null;};
  const grouping={groups:[{structure_group:'G1',n_query_members:3,n_new_members:3}],rule:'coverage',protocol:{method:'butina'}};
  context.getJSON=async url=>url==='/api/new-molecules'?{rows:[measured],files:[]}:url==='/api/groups'?grouping:[reference,measured];
  context.showResult=async languageOnly=>{assert.equal(languageOnly,true);assert.equal(context.showResult.grouping,grouping);};
  context.highlight.source={result:{scope:'new',map_molecule_count:2},cells:new Set(['0,1'])};
  await context.refreshNewMoleculeMutation({job,point:measured});assert.equal(cleared,1);assert.equal(context.highlight.source,null);assert.equal(context.showResult.grouping.groups[0].n_new_members,3);
  const referenceHighlight={result:{scope:'reference',map_molecule_count:1},cells:new Set(['0,0'])};context.highlight.source=referenceHighlight;
  await context.refreshNewMoleculeMutation({job,point:measured});assert.equal(cleared,1);assert.equal(context.highlight.source,referenceHighlight);
 },
 async cropped_patch(){
  const {context,$,node,document}=await setup();const svg=node('svg'),ref=node('g'),incoming=node('g'),source=node('svg'),sourceRef=node('g');
  ref.setAttribute('data-grid-row','0');ref.setAttribute('data-grid-col','0');incoming.setAttribute('data-grid-row','0');incoming.setAttribute('data-grid-col','1');incoming.setAttribute('data-new-molecule','true');
  sourceRef.setAttribute('data-reference-map','true');sourceRef.setAttribute('transform','translate(200,100)');source.append(sourceRef,incoming);svg.append(ref);$('map-stage').append(svg);
  svg.viewBox={baseVal:{x:0,y:0,width:620,height:300}};source.documentElement={viewBox:{baseVal:{width:1000,height:500}}};
  document.importNode=value=>value;context.DOMParser=class{parseFromString(){return source;}};
  context.mapCells=[ref];context.colorMap=()=>{};context.transform={s:2,x:30,y:40};context.patchNewMoleculeCell('<svg/>',measured);
  assert.equal(ref.parentNode,svg);assert.equal(svg.children[0],ref);assert.equal(incoming.getAttribute('transform'),'translate(-200,-100)');
  assert.equal(svg.getAttribute('viewBox'),'-200 -100 1000 500');assert.equal(context.transform.x,-370);assert.equal(context.transform.y,-160);
  // viewBox 原点变化由外层平移抵消，参照图在屏幕的位置保持不变。
  assert.equal(context.transform.x+200*context.transform.s,30);assert.equal(context.transform.y+100*context.transform.s,40);
 },
 async patch(){
  const {context,$,node,document}=await setup();const svg=node('svg'),ref=node('g'),other=node('g'),incoming=node('g'),source=node('svg');
  ref.setAttribute('data-grid-row','0');ref.setAttribute('data-grid-col','0');other.setAttribute('data-grid-row','1');other.setAttribute('data-grid-col','1');incoming.setAttribute('data-grid-row','0');incoming.setAttribute('data-grid-col','1');incoming.setAttribute('data-new-molecule','true');svg.append(ref,other);$('map-stage').append(svg);source.append(incoming);
  const find=svg.querySelector.bind(svg);svg.querySelector=selector=>selector==='[data-reference-map]'?ref:find(selector);
  document.importNode=node=>node;context.DOMParser=class{parseFromString(){return source;}};
  let painted;context.colorMap=cells=>painted=cells;context.mapCells=[ref,other];context.patchNewMoleculeCell('<svg/>',measured);
  assert.equal(ref.parentNode,svg);assert.equal(other.parentNode,svg);assert.equal(svg.children[0],ref);assert.equal(svg.children[1],other);assert.equal(painted.length,1);assert.equal(painted[0],incoming);assert.equal(context.gridSVG,'<svg/>');
  const next=node('g');next.setAttribute('data-grid-row','0');next.setAttribute('data-grid-col','1');next.setAttribute('data-new-molecule','true');incoming.replaceWith=function(value){const at=svg.childNodes.indexOf(this);svg.childNodes[at]=value;value.parentNode=svg;this.parentNode=null;};source.replaceChildren(next);context.patchNewMoleculeCell('<svg>updated</svg>',measured);assert.equal(ref.parentNode,svg);assert.equal(other.parentNode,svg);assert.equal(svg.children[2],next);assert.equal(painted[0],next);
 }
};
(async()=>{const name=process.argv[3];assert.ok(name in scenarios);await scenarios[name]();console.log('U9b '+name+': passed (English and Chinese where applicable)');})().catch(error=>{console.error(error);process.exitCode=1;});
