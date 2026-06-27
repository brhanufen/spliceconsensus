"""Failure-mode analysis: which splice-disrupting variants (SDVs) do the tools miss, and where?
Distance to the nearest splice site is computed from the MFASS minigene construct
(rel_position within [upstream intron | exon | downstream intron]); this matches an
independent genomic-coordinate distance at r=0.997. Recall is compared at a common
operating point (threshold = 90th percentile of neutral scores, ~10% false-positive rate)."""
import csv, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES=f"{ROOT}/results"
rows=[]
with open(f"{ROOT}/data/mfass/snv_data_clean.txt") as f:
    for d in csv.DictReader(f,delimiter="\t"):
        if d.get("category")!="mutant": continue
        sl=d.get("strong_lof")
        try:
            t0=int(d["rel_position"])-1; i1=int(d["intron1_len"]); ex=int(d["exon_len"])
        except: continue
        dss=min(abs(t0-i1), abs(t0-(i1+ex)))          # bp to nearest splice-site junction
        region="exon" if i1<=t0<i1+ex else "intron"
        rows.append({"id":d["id"],"dist_ss":dss,"region":region,
                     "sdv":(1 if sl=="TRUE" else 0 if sl=="FALSE" else np.nan)})
df=pd.DataFrame(rows)
tools=[]
for t in ["pangolin","spliceai","splicetx"]:
    fp=f"{RES}/scores_{t}.csv"
    if os.path.exists(fp): df=df.merge(pd.read_csv(fp).rename(columns={"score":t}),on="id",how="left"); tools.append(t)
df=df.merge(pd.read_csv(f"{ROOT}/data/mfass_labels.csv")[["id","spanr"]],on="id",how="left"); tools.append("spanr")
b=df[df["sdv"].isin([0,1])].dropna(subset=["dist_ss"]).copy(); b["sdv"]=b["sdv"].astype(int)
print(f"evaluable: {len(b)} ({int(b['sdv'].sum())} SDV)")
thr={t: np.nanpercentile(b.loc[b['sdv']==0,t],90) for t in tools}

bins=[(0,2,"0-2"),(3,5,"3-5"),(6,10,"6-10"),(11,20,"11-20"),(21,99,"21-50+")]
def bn(d):
    for lo,hi,nm in bins:
        if lo<=d<=hi: return nm
b["dbin"]=b["dist_ss"].apply(bn); sdv=b[b["sdv"]==1]; order=[nm for _,_,nm in bins]
print(f"\n{'dist(bp)':10}{'#SDV':>6} "+" ".join(f"{t:>9}" for t in tools))
recall={t:[] for t in tools}; counts=[]
for nm in order:
    s=sdv[sdv["dbin"]==nm]; counts.append(len(s)); line=f"{nm:10}{len(s):>6} "
    for t in tools:
        r=(s[t]>=thr[t]).mean() if len(s) else np.nan; recall[t].append(r); line+=f"{r:>9.2f}"
    print(line)
near=sdv[sdv["dist_ss"]<=3]; far=sdv[sdv["dist_ss"]>10]
print(f"\nSDVs <=3bp from splice site: {len(near)} ({100*len(near)/len(sdv):.0f}%) | >10bp: {len(far)} ({100*len(far)/len(sdv):.0f}%)")
for t in tools: print(f"  {t:9}: recall near={near[t].ge(thr[t]).mean():.2f} far={far[t].ge(thr[t]).mean():.2f}")
missed=sdv[~sdv[tools].ge(pd.Series(thr)).any(axis=1)]
print(f"\nSDVs missed by ALL {len(tools)} tools: {len(missed)}/{len(sdv)} ({100*len(missed)/len(sdv):.0f}%)")
print(f"  median dist of missed: {missed['dist_ss'].median():.0f}bp | %>10bp: {100*(missed['dist_ss']>10).mean():.0f}% | %exon: {100*(missed['region']=='exon').mean():.0f}% ({int((missed['region']=='exon').sum())}/{len(missed)})")

plt.rcParams.update({"font.size":13})
fig,ax=plt.subplots(figsize=(8,5.2)); x=np.arange(len(order))
nм={"pangolin":"Pangolin","spliceai":"SpliceAI","splicetx":"SpliceTransformer","spanr":"SPANR"}
co={"pangolin":"#1f77b4","spliceai":"#ff7f0e","splicetx":"#2ca02c","spanr":"#d62728"}
for t in tools: ax.plot(x,recall[t],"o-",lw=2.4,color=co[t],label=nм[t])
ax.set_xticks(x); ax.set_xticklabels([f"{nm}\n(n={c})" for nm,c in zip(order,counts)])
ax.set_xlabel("Distance of SDV to nearest splice site (bp)"); ax.set_ylabel("Recall at ~10% false-positive rate")
ax.set_title("Detection declines with distance from the splice site")
ax.legend(); ax.set_ylim(0,1.02); ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{RES}/fig_failuremode.png",dpi=160)
pd.DataFrame({"dist_bin":order,"n_sdv":counts,**{t:recall[t] for t in tools}}).to_csv(f"{RES}/failuremode_recall.csv",index=False)
print("\nwrote fig_failuremode.png + failuremode_recall.csv")
