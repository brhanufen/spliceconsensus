#!/bin/bash
# Download the full MFASS dataset (for re-scoring / baselines). The benchmark itself
# reproduces from the committed data/mfass_labels.csv without this.
set -e
mkdir -p "$(dirname "$0")/mfass"; cd "$(dirname "$0")/mfass"
BASE=https://raw.githubusercontent.com/KosuriLab/MFASS/master/processed_data/snv
for f in snv_data_clean.txt sdv_simple_list.txt snp_positions_hg38.bed snv_SPANR_scores.txt; do
  echo "downloading $f"; curl -sL "$BASE/$f" -o "$f"
done
echo "done. For re-scoring you also need hg38.fa (UCSC goldenPath) + the Proto tool stack."
