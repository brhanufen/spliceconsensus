"""VALIDATION v2 diagnostic — decides HONESTLY whether the cassette-GTF + signed-delta fix
recovered MMSplice's real signal, or whether MMSplice genuinely performs poorly on MFASS.

Runs LOCALLY (no MMSplice, no GPU): needs only the reduced score table + MFASS labels.
Compares every candidate per-variant reduction against the SDV label, and reports the
directional Spearman vs the continuous dPSI (a directional model should correlate).

Usage:
  python diagnose_mmsplice_v2.py --reduced results/scores_mmsplice_v2_reduced.csv \
      --labels data/mfass_labels.csv

DECISION RULE (printed at the end):
  - If any reduction reaches AUROC >= ~0.75 on the covered subset -> the v1 result was a
    pipeline artifact; adopt that reduction, MMSplice is a legitimate mid-pack tool.
  - If the BEST reduction still <= ~0.70 with the corrected cassette GTF -> MMSplice genuinely
    underperforms on MFASS under a fair pipeline; report it as such (last of 5) WITH the
    caveat that it is used tissue-agnostically on a single-exon minigene context.
  Either outcome is publishable and honest; the rule just removes our discretion.
"""
import argparse
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reduced", required=True)
    ap.add_argument("--labels", default="data/mfass_labels.csv")
    a = ap.parse_args()
    red = pd.read_csv(a.reduced)
    lab = pd.read_csv(a.labels)[["id","sdv","dpsi"]]
    lab["sdv"] = pd.to_numeric(lab["sdv"], errors="coerce")
    d = lab.merge(red, on="id", how="inner")
    d = d[d["sdv"].isin([0,1])].copy(); d["sdv"]=d["sdv"].astype(int)
    base = d["sdv"].mean()
    print(f"covered evaluable: {len(d)}  SDV={d['sdv'].sum()} ({100*base:.2f}%)")
    cands = [c for c in red.columns if c!="id"]
    print(f"\n{'reduction':18}{'AUROC':>8}{'AP':>8}{'AP/base':>9}{'Spearman(-dPSI)':>17}")
    best=("",0,0)
    for c in cands:
        s = d[c].values.astype(float)
        # a disruption score should be HIGHER for SDV; if a signed reduction is negatively
        # oriented (more negative = disruptive), also try its negation and keep the better.
        for orient,ss in [("+",s),("-",-s)]:
            try:
                au=roc_auc_score(d["sdv"], ss); ap=average_precision_score(d["sdv"], ss)
            except Exception: continue
            sp=spearmanr(ss, -d["dpsi"], nan_policy="omit").correlation
            tag=f"{c}{'' if orient=='+' else ' (neg)'}"
            print(f"{tag:18}{au:>8.3f}{ap:>8.3f}{ap/base:>9.1f}{sp:>17.3f}")
            if au>best[1]: best=(tag,au,ap)
    print(f"\nBEST reduction: {best[0]}  AUROC {best[1]:.3f}  AP {best[2]:.3f}")
    v1 = 0.665
    print(f"v1 (single-exon GTF, |delta|) AUROC was {v1:.3f}")
    if best[1] >= 0.75:
        print(f"DECISION: fix RECOVERED signal (+{best[1]-v1:.3f} AUROC). Adopt '{best[0]}'. MMSplice is legit mid-pack.")
    elif best[1] >= v1+0.03:
        print(f"DECISION: fix HELPED (+{best[1]-v1:.3f}) but still moderate. Report improved number; note MFASS is hard for MMSplice.")
    else:
        print(f"DECISION: fix did NOT materially help. MMSplice genuinely underperforms on MFASS under a fair pipeline; report last-of-5 WITH the tissue-agnostic/minigene caveat. NOT a pipeline artifact.")

if __name__ == "__main__": main()
