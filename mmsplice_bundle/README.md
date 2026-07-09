# MMSplice-on-MFASS bundle

Scores MMSplice as the fourth benchmarked predictor in the SpliceConsensus MFASS
benchmark, from MMSplice's standard VCF interface plus a cassette-exon GTF. The scripts
are self-contained; they need a Linux box with the MFASS table and an hg38 FASTA (the
same inputs `src/score_one_tool.py` uses), because MMSplice's stack (TensorFlow, cyvcf2)
has no macOS/py3.8 wheels and cyvcf2 builds from source.

## Canonical pipeline (produces the paper's MMSplice scores: AUROC 0.758, AP 0.256)

MMSplice needs each MFASS exon to sit as an **internal** exon so both of its splice
modules fire. A single-exon GTF leaves roughly half the variants (the deep-intronic
ones) unplaced and gives a misleadingly low score; the cassette GTF below places each
MFASS exon between two dummy flanking exons and restores full 28,972-variant coverage.

```bash
# 0. environment (Linux)
conda create -n mmsplice python=3.7 -y && conda activate mmsplice
pip install mmsplice cyvcf2 kipoiseq pyfaidx pandas

# 1. VCF from MFASS hg38 coordinates (checks ref against FASTA, flips minus-strand
#    alleles, same normalize() logic as score_one_tool.py)
python build_mfass_vcf.py \
    --mfass data/mfass/snv_data_clean.txt \
    --hg38  ref/hg38.fa \
    --out_vcf mfass.vcf --out_gtf /dev/null   # the single-exon GTF here is NOT used

# 2. cassette-exon GTF: each MFASS exon as the internal exon of a 3-exon transcript
python validation_v2/build_mfass_gtf_v2.py \
    --mfass data/mfass/snv_data_clean.txt \
    --out_gtf mfass_exons_v2.gtf --flank 300 --dummy 50

# 3. sort + bgzip + index (bcftools + tabix; both on HCC modules)
bcftools sort mfass.vcf -Oz -o mfass.sorted.vcf.gz
tabix -p vcf mfass.sorted.vcf.gz

# 4. score with MMSplice -> results/scores_mmsplice_v2_reduced.csv
#    (emits four candidate reductions; the paper uses max_abs_delta = |delta logit Psi|)
python validation_v2/score_mmsplice_v2.py \
    --gtf mfass_exons_v2.gtf --vcf mfass.sorted.vcf.gz --fasta ref/hg38.fa \
    --out_raw results/scores_mmsplice_v2_raw.csv

# 5. take the max_abs_delta column as the per-variant score -> results/scores_mmsplice.csv
#    (id, score), then fold into the benchmark (see INTEGRATE below):
python src/analyze.py         # per-method table + figures, incl. mmsplice
python src/failure_mode.py    # blind-spot recall, "missed by all 5"
```

## Reduction choice
`score_mmsplice_v2.py` emits four per-variant reductions of the signed delta-logit-Psi.
`max_abs_delta` (the magnitude) gives the best discrimination (AUROC 0.758); the signed
reductions all score lower (~0.64). The benchmark uses `max_abs_delta`.

## Files
- `build_mfass_vcf.py`: MFASS table -> VCF (self-contained, no GENCODE needed)
- `validation_v2/build_mfass_gtf_v2.py`: cassette-exon GTF builder (the GTF the paper uses)
- `validation_v2/score_mmsplice_v2.py`: MMSplice scoring; emits raw + reduced tables
- `validation_v2/diagnose_mmsplice_v2.py`: local sanity check on scores vs labels
- `INTEGRATE.md`: the exact edits that add mmsplice to analyze.py + failure_mode.py

## Validation (checked during prep, no FASTA required)
- MFASS construct span == intron1 + exon + intron2 for 28,962/28,972 (100.0%).
- Every variant's hg38 position falls within its construct span (28,972/28,972).
- Reconstructed genomic exon boundaries agree with MFASS's construct-based region label
  for 28,961/28,972 (99.96%), so the GTF places each variant in the correct exon.
- Strand handling copied from score_one_tool.py normalize() (paper's r=0.997 concordance).

## What MMSplice adds to the paper
MMSplice is the one model built for modular exonic/intronic effects rather than
splice-site recognition. On this benchmark it ranks fourth (AUROC 0.758, AP 0.256),
above the legacy SPANR model, and shows the **same** distance-dependent decline into the
exon interior as the splice-site models. That makes the paper's exon-interior blind-spot
finding a shared limitation of current predictors rather than an artifact of
splice-site-centric design.
