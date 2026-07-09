"""VALIDATION v2 scorer: keep the SIGN of delta_logit_psi and emit every candidate
reduction, so we can test which (if any) recovers MMSplice's real signal. No |.| baked in.

Emits results/scores_mmsplice_v2_raw.csv with columns:
    id, delta_logit_psi_signed, abs_delta, ref_acceptor_diff?, ...   (whatever predict_all_table gives)
plus a per-variant reduced table for each candidate score. The DECISION is made by
diagnose_mmsplice_v2.py against the MFASS labels, not here.

Run on HCC (same env as v1):
  python score_mmsplice_v2.py --gtf mfass_exons_v2.gtf --vcf mfass.sorted.vcf.gz --fasta ref/hg38.fa \
      --out_raw results/scores_mmsplice_v2_raw.csv
"""
import argparse, os
import numpy as np, pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True); ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out_raw", default="results/scores_mmsplice_v2_raw.csv")
    a = ap.parse_args()
    from mmsplice.vcf_dataloader import SplicingVCFDataloader
    from mmsplice import MMSplice, predict_all_table
    from mmsplice.utils import max_varEff

    dl = SplicingVCFDataloader(a.gtf, a.fasta, a.vcf, tissue_specific=False)
    model = MMSplice()
    df = predict_all_table(model, dl, pathogenicity=False, splicing_efficiency=False)
    print("predict_all_table columns:", list(df.columns))
    print("rows (variant x exon):", len(df))
    # Keep the SIGNED delta_logit_psi. Reduce per variant TWO ways for the diagnostic:
    #   (a) max |delta|  (v1 behavior)           (b) min signed delta (most-skipping)
    idcol = "ID" if "ID" in df.columns else df.columns[0]
    df = df.rename(columns={idcol: "id"})
    df["abs_delta"] = df["delta_logit_psi"].abs()
    # per-variant reductions
    g = df.groupby("id")
    red = pd.DataFrame({
        "id": list(g.groups.keys()),
        "max_abs_delta":  g["delta_logit_psi"].apply(lambda s: s.abs().max()).values,    # v1
        "min_signed":     g["delta_logit_psi"].min().values,                              # most negative = most skipping
        "max_signed":     g["delta_logit_psi"].max().values,
        "signed_at_maxabs": g["delta_logit_psi"].apply(lambda s: s.iloc[s.abs().values.argmax()]).values,
    })
    os.makedirs(os.path.dirname(a.out_raw) or ".", exist_ok=True)
    df.to_csv(a.out_raw, index=False)
    red.to_csv(a.out_raw.replace("_raw.csv","_reduced.csv"), index=False)
    print(f"wrote {a.out_raw} ({len(df)} rows) and reduced table ({len(red)} variants)")

if __name__ == "__main__": main()
