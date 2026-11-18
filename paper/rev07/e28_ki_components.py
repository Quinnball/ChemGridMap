import os
import pandas as pd, numpy as np, glob
from pathlib import Path
from scipy import stats
CL=("inactive","medium","active")
REF=Path(os.environ.get("E28_REF", Path(__file__).resolve().parent / "e28_ref")) / "primary"
docs={}
for d in REF.glob("*_Ki"):
    ld=pd.read_csv(d/"later_documents.csv"); r=pd.read_csv(d/"reference.csv")
    s=set(";".join(ld.documents.dropna()).split(";"))|set(";".join(r.documents_all.dropna()).split(";"))
    docs[d.name]=s-{""}
names=sorted(docs); parent={n:n for n in names}
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
for i,a in enumerate(names):
    for b in names[i+1:]:
        if docs[a]&docs[b]: parent[find(a)]=find(b)
comp={n:find(n) for n in names}
groups={}
for n,c in comp.items(): groups.setdefault(c,[]).append(n)
print("Ki datasets",len(names),"components (any shared document):",sorted(len(v) for v in groups.values()))
c=pd.read_csv("compounds_annotated_primary.csv.gz"); c=c[c.endpoint.eq("Ki")&c.measured_class.isin(CL)]
c["component"]=c.dataset.map(comp)
rows=[]
for g,t in c.groupby("component"):
    cc=t[t.ev.isin(CL)]
    if len(cc)==0: continue
    pm=cc.measured_class.value_counts(normalize=True); pc=cc.ev.value_counts(normalize=True)
    n1=cc.measured_class.eq(cc.ev).mean()-sum(pm.get(x,0)*pc.get(x,0) for x in CL)
    e6=t[t.tanimoto>=0.6]; ag=e6.measured_class.eq(e6.nearest_reference_class)
    n2=ag[e6.group.eq("consistent")].mean()-ag[e6.group.eq("unsettled")].mean()
    far=t[t.group.eq("none")]; n3=ag[e6.group.eq("consistent")].mean()-far.measured_class.eq(far.nearest_reference_class).mean()
    rows.append(dict(component=g,datasets=t.dataset.nunique(),consistent=len(cc),n1=n1,n2=n2,n3=n3))
r=pd.DataFrame(rows); print(r.round(3).to_string())
for col in ["n1","n2","n3"]:
    v=r[col].dropna().values; print(col,"positive",(v>0).sum(),"of",len(v), "exact one-sided P (sign-rank)", stats.wilcoxon(v,alternative="greater").pvalue if len(v)>=2 else None)
r.to_csv("ki_components_primary.csv",index=False)
