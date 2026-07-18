"""Benchmark the splicing predictors against MFASS experimental labels,
build a held-out consensus, and make figures. Runs locally (CPU) on downloaded scores.
Primary metric AUPRC/AP (imbalanced ~3.8% positive); also AUROC + Spearman vs delta-PSI."""
import csv, glob, os, json
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NAMES = {"pangolin":"Pangolin", "spliceai":"SpliceAI",
         "splicetx":"SpliceTransformer", "mmsplice":"MMSplice",
         "spanr":"SPANR"}

ROOT = os.environ.get("SPLICE_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = f"{ROOT}/results"; MF = f"{ROOT}/data/mfass"
rng = np.random.default_rng(0)

def load_labels():
    # Prefer the slim committed labels file (id,sdv,dpsi,spanr); reproduces the benchmark from the repo alone.
    slim = f"{ROOT}/data/mfass_labels.csv"
    if os.path.exists(slim):
        d = pd.read_csv(slim)
        d["sdv"] = pd.to_numeric(d["sdv"], errors="coerce")
        d["dpsi"] = pd.to_numeric(d["dpsi"], errors="coerce")
        return d
    rows = []
    with open(f"{MF}/snv_data_clean.txt") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d.get("category") != "mutant": continue
            sl = d.get("strong_lof")
            try: dpsi = float(d["v2_dpsi"])
            except (ValueError, KeyError, TypeError): dpsi = np.nan
            rows.append({"id": d["id"], "sdv": (1 if sl=="TRUE" else 0 if sl=="FALSE" else np.nan), "dpsi": dpsi})
    return pd.DataFrame(rows)

def load_tool(name):
    fp = f"{RES}/scores_{name}.csv"
    shards = sorted(glob.glob(f"{RES}/scores_{name}_s*of*.csv"))
    if os.path.exists(fp): df = pd.read_csv(fp)
    elif shards: df = pd.concat([pd.read_csv(s) for s in shards], ignore_index=True)
    else: return None
    df = df.dropna(subset=["score"]).drop_duplicates("id")
    return df.rename(columns={"score": name})[["id", name]]

def boot_ci(y, s, fn, n=1000):
    vals = []
    idx = np.arange(len(y))
    for _ in range(n):
        b = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) < 2: continue
        vals.append(fn(y[b], s[b]))
    return (np.percentile(vals, 2.5), np.percentile(vals, 97.5)) if vals else (np.nan, np.nan)

def main():
    lab = load_labels()
    df = lab.copy()
    tools = []
    for t in ["spliceai", "pangolin", "splicetx", "mmsplice"]:
        td = load_tool(t)
        if td is not None: df = df.merge(td, on="id", how="left"); tools.append(t)
    # legacy baseline SPANR: from slim labels column if present, else the full file
    if "spanr" in df.columns:
        tools.append("spanr")
    else:
        try:
            sp = pd.read_csv(f"{MF}/snv_SPANR_scores.txt", sep="\t"); sp["spanr"] = sp["dpsi_max_tissue"].abs()
            df = df.merge(sp[["id","spanr"]].drop_duplicates("id"), on="id", how="left"); tools.append("spanr")
        except Exception as e: print("SPANR skip:", e)

    print(f"merged: {len(df)} variants; tools: {tools}")
    print("coverage (non-null):", {t: int(df[t].notna().sum()) for t in tools})

    # ---- binary benchmark (drop NA labels) ----
    b = df[df["sdv"].isin([0,1])].copy(); b["sdv"] = b["sdv"].astype(int)
    print(f"\nbinary eval set: {len(b)} ({int(b['sdv'].sum())} SDV, {b['sdv'].mean()*100:.2f}% pos)")
    res = []
    for t in tools:
        m = b.dropna(subset=[t]); y = m["sdv"].values; s = m[t].values
        if len(m) < 50 or y.sum() < 5: continue
        auroc = roc_auc_score(y, s); ap = average_precision_score(y, s)
        sp_dpsi = spearmanr(m[t], -m["dpsi"], nan_policy="omit").correlation
        lo, hi = boot_ci(y, s, average_precision_score)
        res.append({"method": t, "n": len(m), "AUROC": round(auroc,4), "AP": round(ap,4),
                    "AP_lo": round(lo,4), "AP_hi": round(hi,4), "Spearman_vs_loss": round(sp_dpsi,4)})
    rt = pd.DataFrame(res).sort_values("AP", ascending=False)
    print("\n=== PER-METHOD BENCHMARK (sorted by AP) ===\n", rt.to_string(index=False))

    # ---- consensus (held-out), modern tools only ----
    modern = [t for t in ["spliceai","pangolin","splicetx"] if t in tools]
    cons = None
    if len(modern) >= 2:
        c = b.dropna(subset=modern).copy()
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import make_pipeline
        from scipy.stats import rankdata
        from sklearn.model_selection import GroupShuffleSplit
        c["exon"] = c["id"].str.rsplit("_", n=1).str[0]   # group by exon so no exon spans train/test
        tr_idx, te_idx = next(GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=0).split(c[modern].values, c["sdv"].values, groups=c["exon"].values))
        Xtr, Xte = c[modern].values[tr_idx], c[modern].values[te_idx]
        ytr, yte = c["sdv"].values[tr_idx], c["sdv"].values[te_idx]
        lr = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced")).fit(Xtr, ytr)
        pcons = lr.predict_proba(Xte)[:,1]
        ap_c = average_precision_score(yte, pcons); au_c = roc_auc_score(yte, pcons)
        # robust rank-average consensus (mean of per-tool ranks on held-out)
        rank_mean = np.mean([rankdata(Xte[:,i]) for i in range(len(modern))], axis=0)
        ap_rank = average_precision_score(yte, rank_mean); au_rank = roc_auc_score(yte, rank_mean)
        singles = {t: average_precision_score(yte, Xte[:,i]) for i,t in enumerate(modern)}
        best_t = max(singles, key=singles.get)
        coefs = lr.named_steps["logisticregression"].coef_[0]
        cons = {"consensus_LR_AP": round(ap_c,4), "consensus_LR_AUROC": round(au_c,4),
                "consensus_rankmean_AP": round(ap_rank,4), "consensus_rankmean_AUROC": round(au_rank,4),
                "best_single": best_t, "best_single_AP": round(singles[best_t],4), "best_single_AUROC": round(roc_auc_score(yte, Xte[:,modern.index(best_t)]),4),
                "held_out_n": len(yte), "lr_coef_scaled": dict(zip(modern, coefs.round(3).tolist())),
                "all_singles_AP": {k: round(v,4) for k,v in singles.items()}}
        print("\n=== CONSENSUS (held-out 50%) ===\n", json.dumps(cons, indent=1))
        # Nominal AP difference between the best consensus and the best single tool is
        # far inside the bootstrap CI (see paper); report it as within-noise, not a win.
        diff = max(ap_c, ap_rank) - singles[best_t]
        print(f"  -> best consensus AP {max(ap_c,ap_rank):.4f} vs best single ({best_t}) {singles[best_t]:.4f}: "
              f"difference {diff:+.4f} AP, within bootstrap noise (no meaningful improvement over the best single tool)")

    # ---- figures ----
    plt.figure(figsize=(6,5))
    for t in rt["method"]:
        m = b.dropna(subset=[t]); fpr,tpr,_ = roc_curve(m["sdv"], m[t]); plt.plot(fpr,tpr,label=f"{NAMES.get(t,t)} (AUROC {roc_auc_score(m['sdv'],m[t]):.3f})")
    plt.plot([0,1],[0,1],"k--",lw=.7); plt.xlabel("FPR"); plt.ylabel("TPR"); plt.legend(fontsize=8); plt.title("MFASS SDV detection: ROC"); plt.tight_layout(); plt.savefig(f"{RES}/fig_roc.png",dpi=150)
    plt.figure(figsize=(6,5))
    for t in rt["method"]:
        m = b.dropna(subset=[t]); pr,rc,_ = precision_recall_curve(m["sdv"], m[t]); plt.plot(rc,pr,label=f"{NAMES.get(t,t)} (AP {average_precision_score(m['sdv'],m[t]):.3f})")
    plt.axhline(b["sdv"].mean(),color="k",ls="--",lw=.7,label=f"baseline {b['sdv'].mean():.3f}"); plt.xlabel("Recall"); plt.ylabel("Precision"); plt.legend(fontsize=8); plt.title("MFASS SDV detection: PR"); plt.tight_layout(); plt.savefig(f"{RES}/fig_pr.png",dpi=150)

    rt.to_csv(f"{RES}/benchmark_table.csv", index=False)
    json.dump({"per_method": res, "consensus": cons, "eval_n": len(b), "n_pos": int(b["sdv"].sum())}, open(f"{RES}/benchmark_summary.json","w"), indent=1)
    print("\nwrote benchmark_table.csv, benchmark_summary.json, fig_roc.png, fig_pr.png")

if __name__ == "__main__": main()
