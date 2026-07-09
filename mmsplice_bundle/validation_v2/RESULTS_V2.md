# MMSplice validation v2 — result (run 2026-07-05, diagnosed 2026-07-07)

## Verdict: v1 was substantially a PIPELINE ARTIFACT. MMSplice is legit mid-pack, not last-of-5.

Rule-based decision (diagnose_mmsplice_v2.py, threshold AUROC >= 0.75):
> DECISION: fix RECOVERED signal (+0.093 AUROC). Adopt 'max_abs_delta'. MMSplice is legit mid-pack.

## What the fix did
| | v1 (single-exon GTF) | v2 (cassette GTF) |
|---|---|---|
| variants scored | 14,595 (50.3% of 27,733) | **28,972 (100%)** |
| covered evaluable | 13,960 | **27,733** |
| AUROC | 0.665 | **0.758** |
| AP (auPRC) | 0.093 | **0.256** |

The single-exon GTF was the real culprit: with each MFASS exon isolated, MMSplice skipped
its acceptor/donor modules and (a) failed to score ~half the variants and (b) produced
near-random scores on the rest. Making each exon INTERNAL in a 3-exon cassette (real genomic
flanking introns) let both splice modules fire on all variants.

## The |delta| reduction was NOT the problem
Best reduction is `max_abs_delta` — the SAME sign-agnostic |delta_logit_psi| v1 used. The
signed reductions (min_signed, etc.) did WORSE (AUROC ~0.64-0.65). So directionality was a
red herring; the GTF context was the whole story.

Full reduction sweep (covered evaluable = 27,733; base rate 3.79%):
```
reduction            AUROC      AP  AP/base  Spearman(-dPSI)
max_abs_delta        0.758   0.256      6.8            0.090   <- best
signed_at_maxabs(neg)0.648   0.213      5.6            0.130
min_signed (neg)     0.649   0.213      5.6            0.130
max_signed (neg)     0.643   0.210      5.5            0.130
```

## Independent check vs published MMSplice-on-MFASS
Published MMSplice auPRC on the SAME MFASS set (27,733 variants / 1,050 SDVs): ~0.361
(Yin et al., Brief. Bioinform. 2022, bbac334; MMSplice+SpliceAI "MMAIpsi" combo reaches 0.472).
- v1 AP 0.093 was ~4x too low (26% of published) -> clearly broken.
- v2 AP 0.256 is 71% of published -> same ballpark, no longer near-random.
- Residual gap is expected: we score on a synthetic single-cassette minigene context with a
  tissue-agnostic model and no reference-PSI/psi term, not the variant's real multi-exon
  transcript. So v2 modestly under-extracts but is faithful enough to place MMSplice correctly.

## Where MMSplice lands in the 5-tool benchmark (using v2 max_abs_delta, all on 27,733)
| method | AUROC | AP |
|---|---|---|
| pangolin | 0.888 | 0.421 |
| spliceai | 0.819 | 0.321 |
| splicetx | 0.786 | 0.317 |
| **mmsplice (v2)** | **0.758** | **0.256** |
| spanr | 0.748 | 0.228 |

MMSplice moves from LAST (v1) to **4th of 5** — above the SPANR legacy baseline, below the
three modern deep-learning tools. Report it as mid-pack, not last-of-5.

## Recommended next step (not yet done — awaiting sign-off)
Replace `results/scores_mmsplice.csv` with the v2 max_abs_delta scores (id-keyed) and re-run
`analyze.py` + `failure_mode.py`, so the manuscript's per-method table and the exon-interior
blind-spot analysis use the fair (cassette-GTF) MMSplice numbers. The blind-spot claim should
be re-checked against v2 — MMSplice now scores all variants, so its exon-interior recall is a
fair test rather than a coverage artifact.

## Files
- results/scores_mmsplice_v2_reduced.csv  (MFASS-id keyed; 4 candidate reductions per variant)
- results/scores_mmsplice_v2_reduced_mfassid.csv (same, on HCC)
- results/scores_mmsplice_v2_raw.csv       (HCC; per variant x exon, all MMSplice module diffs)
- results/diagnose_mmsplice_v2.txt         (the diagnostic printout)
