"""Score MFASS variants with MMSplice -> results/scores_mmsplice.csv (columns: id,score).

Uses MMSplice's standard VCF interface (deltaLogitPSI). The disruption-magnitude score
is |delta_logit_psi| per variant, matching how the other tools are reduced to a single
magnitude in this benchmark (SpliceAI max|delta|, Pangolin max gain/loss magnitude, etc.).
When a variant maps to multiple exon predictions, the max-magnitude effect is kept
(mmsplice.utils.max_varEff), the tool's own recommended per-variant reduction.

Run on HCC (GPU optional; MMSplice is small, CPU is fine):
    conda create -n mmsplice python=3.7 -y && conda activate mmsplice
    pip install mmsplice cyvcf2 kipoiseq        # full stack (needs a Linux box; builds cyvcf2)
    python build_mfass_vcf.py --mfass <MFASS.txt> --hg38 <hg38.fa> --out_vcf mfass.vcf --out_gtf mfass_exons.gtf
    bcftools sort mfass.vcf -Oz -o mfass.sorted.vcf.gz && tabix -p vcf mfass.sorted.vcf.gz
    python score_mmsplice.py --gtf mfass_exons.gtf --vcf mfass.sorted.vcf.gz --fasta <hg38.fa> --out results/scores_mmsplice.csv
"""
import argparse, os
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--vcf", required=True, help="sorted, bgzipped, tabix-indexed")
    ap.add_argument("--fasta", required=True, help="hg38 (must match VCF contig names, chr*)")
    ap.add_argument("--out", default="results/scores_mmsplice.csv")
    a = ap.parse_args()

    from mmsplice.vcf_dataloader import SplicingVCFDataloader
    from mmsplice import MMSplice, predict_all_table
    from mmsplice.utils import max_varEff

    dl = SplicingVCFDataloader(a.gtf, a.fasta, a.vcf, tissue_specific=False)
    model = MMSplice()
    # predict_all_table returns per (variant, exon) rows incl. 'delta_logit_psi' and 'ID'
    df = predict_all_table(model, dl, pathogenicity=False, splicing_efficiency=False)
    # reduce to one row per variant by max |delta_logit_psi| (tool-recommended reduction)
    df = max_varEff(df)  # adds/uses 'delta_logit_psi'; keeps max-effect exon per variant
    # 'ID' column holds the VCF ID we set to the MFASS variant id
    idcol = "ID" if "ID" in df.columns else ("variant" if "variant" in df.columns else df.columns[0])
    out = pd.DataFrame({
        "id":   df[idcol].astype(str),
        "score": df["delta_logit_psi"].abs().values,   # disruption magnitude
    }).dropna().drop_duplicates("id")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out}: {len(out)} variants scored")
    print(out.head())

if __name__ == "__main__":
    main()
