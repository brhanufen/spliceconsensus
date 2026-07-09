# Folding MMSplice into the existing benchmark

After `score_mmsplice.py` writes `results/scores_mmsplice.csv`, MMSplice becomes a
5th tool. The two analysis scripts pick tools up by name; add "mmsplice" in each place.

## src/analyze.py  (3 edits)
1. Line ~56, add mmsplice to the loaded tools:
     for t in ["spliceai", "pangolin", "splicetx"]:
   ->
     for t in ["spliceai", "pangolin", "splicetx", "mmsplice"]:

   (mmsplice is loaded by the same load_tool(); it will appear in the per-method
    benchmark table and figures automatically.)

2. (Optional) include mmsplice in the consensus. Line ~87:
     modern = [t for t in ["spliceai","pangolin","splicetx"] if t in tools]
   ->
     modern = [t for t in ["spliceai","pangolin","splicetx","mmsplice"] if t in tools]

   Decide deliberately: the paper's consensus story is "modern deep-learning tools";
   MMSplice is modular/linear. Cleanest is to KEEP the consensus as the 3 modern
   tools and report MMSplice only in the per-tool ranking + failure-mode. If you add
   it, re-check the "consensus does not beat best single" sentence still holds.

## src/failure_mode.py  (2 edits)
1. Line ~25, add mmsplice to the load loop:
     for t in ["pangolin","spliceai","splicetx"]:
   ->
     for t in ["pangolin","spliceai","splicetx","mmsplice"]:

2. Lines ~52-53, give it a label + colour in the figure dicts:
     nм={"pangolin":"Pangolin","spliceai":"SpliceAI","splicetx":"SpliceTransformer","spanr":"SPANR"}
     co={"pangolin":"#1f77b4","spliceai":"#ff7f0e","splicetx":"#2ca02c","spanr":"#d62728"}
   -> add:
     "mmsplice":"MMSplice"      to nм
     "mmsplice":"#9467bd"       to co   (purple, unused so far)

## What to report
- MMSplice's row in the per-method table (AUROC, AP, CI): where does it rank?
- MMSplice's line in the failure-mode plot, KEY QUESTION for the paper:
  does it share the exon-interior blind spot, or does it recover exon-interior SDVs
  the splice-site models miss? Recompute the "missed by ALL tools" count with
  MMSplice included (it becomes "missed by all 5").
- If MMSplice shares the blind spot -> paper's central claim gets STRONGER
  ("even the model built for exonic regulation misses these").
  If it partially fills the blind spot -> soften "all predictors" to "all splice-site
  models", and report MMSplice as the partial exception. Either outcome is publishable
  and better than excluding it.

The manuscript edits these changes support (five-tool table, blind-spot paragraph, scoring description) are already applied in the committed `docs/SpliceConsensus_Report.tex`.
