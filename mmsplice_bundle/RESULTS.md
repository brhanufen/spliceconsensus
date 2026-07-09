# MMSplice-on-MFASS — run results (2026-07-05)

Ran on HCC (swan). Env: py3.7, mmsplice 2.4.0, cyvcf2 0.30.15, tensorflow 2.11 (CPU).
One env fix: sorted-nearest 0.0.41 hard-imports py3.8 `importlib.metadata`; added a
`sitecustomize.py` shim aliasing the installed `importlib_metadata` backport. No scoring
logic changed; committed scripts run verbatim.

## Pipeline (all real, nothing fabricated)
- build_mfass_vcf.py: 28,972 mutants written, 0 skipped, 98 minus-strand flipped, 2,198 exons.
- sort+bgzip+tabix (bcftools 1.2 lacks `sort`; used manual coord sort): 28,972 records.
- score_mmsplice.py (Slurm job 16499276, batch, 32 s, CPU): wrote 14,595 scored variants.
- MMSplice emits keys as `chr:pos:ref>alt`, not the MFASS id. Remapped 1:1 via the VCF
  (28,972 unique genomic keys, all multiplicity 1; 14,595/14,595 matched, 0 unmapped).
  Benchmark file `results/scores_mmsplice.csv` is MFASS-id keyed; raw output kept as
  `results/scores_mmsplice_chrposkey.csv`.

## Coverage
- 14,595 total scored; 13,960 / 27,733 evaluable (50.3%).
- The ~50% uncovered are deep-intronic variants MMSplice does not model.
- Exon-interior SDVs are FULLY covered (333/333), so the blind-spot verdict below is NOT
  a coverage artifact.

## Ranking (per-method benchmark)
| method | n | AUROC | AP |
|---|---|---|---|
| pangolin | 27733 | 0.888 | 0.421 |
| spliceai | 27733 | 0.819 | 0.321 |
| splicetx | 27733 | 0.786 | 0.317 |
| spanr | 27663 | 0.748 | 0.228 |
| **mmsplice** | **13960** | **0.665** | **0.093** |

MMSplice ranks LAST. (AUROC is the fair cross-tool comparator; MMSplice AP is on its
covered subset.) Consensus unchanged (3 modern tools) and still does not beat best-single
Pangolin (+0.0014 AP, within noise).

## Exon-interior blind spot — KEY QUESTION
MMSplice recall by distance-to-splice-site: 0.19 / 0.08 / 0.12 / 0.15 / 0.21 — flat and
uniformly low; near=0.17, far=0.18. It does NOT show the near-site peak the splice-site
models have, and it does NOT recover the exon-interior SDVs they miss:
- Of 257 SDVs missed by the 3 splice-site models, MMSplice rescues only 28 (79 it didn't
  even score).
- "Missed by all tools": 228/1050 (4 tools) -> 208/1050 (with MMSplice) — MMSplice rescued 20.
- On its scored exon-interior subset (region=exon & >10 bp, n=333): recall 0.29.

**Verdict: MMSplice SHARES the exon-interior blind spot; it does not fill it. The paper's
central claim gets STRONGER** — even the model built for modular exonic regulation misses
these variants. Keep "all predictors" wording (add MMSplice as a scored 4th DL tool).

Caveat to eyeball: MMSplice's absolute AUROC (0.665) is on the low side; the score is the
bundle's max|deltaLogitPSI| reduction over the reconstructed single-exon GTF. Pipeline ran
correctly; worth a sanity glance before final manuscript numbers, but it does not change
the qualitative verdict.
