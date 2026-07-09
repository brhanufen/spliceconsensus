# MMSplice validation v2 — is the weak v1 result real, or a pipeline artifact?

## Why this exists
The first MMSplice run scored AUROC 0.665 / AP 0.093 — near-random for a tool designed
for exactly MFASS's variant class. Before that number goes in the paper (as a claim about
a real published tool), we must rule out two pipeline suspects:

1. **Single-exon GTF (leading suspect).** v1 made each MFASS exon its OWN single-exon
   transcript. MMSplice's acceptor/donor modules need the exon to be INTERNAL (flanked by
   other exons) to score both splice sites; a lone exon can make it skip its splice modules
   entirely -> degenerate scores. build_mfass_gtf_v2.py fixes this: each MFASS exon E becomes
   the internal exon of a 3-exon cassette transcript, with dummy exons 300 bp away so the
   introns bounding E are REAL genomic sequence.

2. **|delta| magnitude reduction.** v1 used max|delta_logit_psi|. MMSplice's disruption
   signal is DIRECTIONAL (negative = more skipping). score_mmsplice_v2.py keeps the sign and
   emits every candidate reduction; diagnose_mmsplice_v2.py picks the best against the labels.

## Steps (HCC, same env as v1)
```bash
# 1. corrected GTF (VCF unchanged — reuse mfass.sorted.vcf.gz from v1)
python build_mfass_gtf_v2.py --mfass data/mfass/snv_data_clean.txt \
    --out_gtf mfass_exons_v2.gtf --flank 300 --dummy 50

# 2. re-score keeping the sign (~30 s, same Slurm setup as v1)
python score_mmsplice_v2.py --gtf mfass_exons_v2.gtf --vcf mfass.sorted.vcf.gz \
    --fasta ref/hg38.fa --out_raw results/scores_mmsplice_v2_raw.csv
#    -> also writes results/scores_mmsplice_v2_reduced.csv

# 3. decide (runs LOCALLY too — only needs the reduced table + labels)
python diagnose_mmsplice_v2.py --reduced results/scores_mmsplice_v2_reduced.csv \
    --labels data/mfass_labels.csv
```

## Second, independent check (do if convenient)
MMSplice's own papers reported an MFASS metric. If you can find their published auPRC/auROC
on MFASS, compare it to our best v2 reduction on the SAME variants. If we land near their
number, our pipeline is faithful. If we're far below it, the pipeline is still under-
extracting and we should NOT report MMSplice as last-of-5.

## The decision is rule-based, not discretionary (integrity)
diagnose_mmsplice_v2.py prints one of three verdicts by AUROC:
- >= 0.75  -> v1 was a pipeline artifact; adopt the winning reduction; MMSplice is legit mid-pack.
- +0.03 to <0.75 -> fix helped but MFASS is genuinely hard for MMSplice; report improved number.
- < v1+0.03 -> fix did NOT help; MMSplice genuinely underperforms under a FAIR pipeline;
               report last-of-5 WITH the tissue-agnostic/minigene caveat. This is now defensible.

## What to send back to the analysis session
- results/scores_mmsplice_v2_reduced.csv (the per-variant reductions)
- the diagnose_mmsplice_v2.py printout (which verdict, and the AUROC/AP of the best reduction)
- MMSplice's published MFASS number if you found one
Only after this do we lock MMSplice's numbers into the manuscript — in whichever direction
the data actually points. No number goes in until the plausibility check is passed or the
underperformance is confirmed to be real.
