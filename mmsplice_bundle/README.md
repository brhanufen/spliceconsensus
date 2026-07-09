# MMSplice-on-MFASS bundle — run-ready for HCC (GPU box)

Adds MMSplice as a 5th predictor to the SpliceConsensus benchmark, using MMSplice's
standard VCF interface. Everything here is prepared and validated; it needs a Linux box
with the MFASS table + hg38 FASTA (the same inputs score_one_tool.py already uses on HCC).

## Why this runs on HCC and not in the prep sandbox
MMSplice's stack (TensorFlow, cyvcf2) has no macOS/py3.8 wheels and cyvcf2 builds from
source — fine on HCC Linux, impossible in the prep sandbox. The scripts below were
written against MMSplice's documented API and the coordinate/strand logic was validated
against the MFASS table (99.96% exon-region concordance; see "Validation" below).

## Inputs needed on HCC (paths you already have from the Proto scoring run)
- MFASS table:  data/mfass/snv_data_clean.txt   (32k rows; the repo has it)
- hg38 FASTA:   ref/hg38.fa                      (same file score_one_tool.py used)

## Steps
```bash
# 0. environment (Linux)
conda create -n mmsplice python=3.7 -y && conda activate mmsplice
pip install mmsplice cyvcf2 kipoiseq pyfaidx pandas

# 1. build the VCF + exon GTF from MFASS's hg38 coordinates (checks ref against FASTA,
#    flips minus-strand alleles — same normalize() logic as score_one_tool.py)
python build_mfass_vcf.py \
    --mfass data/mfass/snv_data_clean.txt \
    --hg38  ref/hg38.fa \
    --out_vcf mfass.vcf --out_gtf mfass_exons.gtf

# 2. sort + bgzip + index (needs bcftools + tabix; both on HCC modules)
bcftools sort mfass.vcf -Oz -o mfass.sorted.vcf.gz
tabix -p vcf mfass.sorted.vcf.gz

# 3. score with MMSplice -> results/scores_mmsplice.csv  (id,score = |delta_logit_psi|)
python score_mmsplice.py \
    --gtf mfass_exons.gtf --vcf mfass.sorted.vcf.gz --fasta ref/hg38.fa \
    --out results/scores_mmsplice.csv

# 4. fold into the benchmark (see INTEGRATE.md for the exact 5 edits)
python src/analyze.py         # per-method table + figures, now incl. mmsplice
python src/failure_mode.py    # blind-spot recall, now "missed by all 5"
```

## Files
- build_mfass_vcf.py  — MFASS table -> VCF + exon GTF (self-contained, no GENCODE needed)
- score_mmsplice.py   — MMSplice deltaLogitPSI scoring -> scores_mmsplice.csv
- INTEGRATE.md        — the exact edits to analyze.py + failure_mode.py, and what to report
- README.md           — this file

## Validation done during prep (in the sandbox, no FASTA required)
- MFASS construct span == intron1+exon+intron2 for 28,962/28,972 (100.0%).
- Every variant's hg38 position falls within its construct span (28,972/28,972).
- Reconstructed genomic exon boundaries agree with MFASS's construct-based exon/intron
  region label for 28,961/28,972 (99.96%) — so the GTF puts each variant in the right exon.
- Strand handling copied verbatim from score_one_tool.py normalize() (paper's r=0.997).

## The scientific question this answers (for the paper)
MMSplice is the model designed for modular exonic/intronic effects — the class of variant
in the paper's "exon-interior blind spot". After step 4, check the failure-mode plot:
- If MMSplice ALSO degrades in the exon interior -> the blind-spot claim gets STRONGER.
- If it recovers exon-interior SDVs the splice-site models miss -> report it as the partial
  exception and soften "all predictors" to "all splice-site-based predictors".
Either way the paper is better with MMSplice in than with it excluded.
