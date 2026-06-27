"""Labeled correctness tests (synthetic arrays here ONLY — not part of any reported result)."""
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

def test_metrics_sanity():
    y = np.array([0,0,0,1,1])
    assert roc_auc_score(y, np.array([0.1,0.2,0.3,0.9,0.8])) == 1.0      # perfect ranking
    assert abs(average_precision_score(y, y.astype(float)) - 1.0) < 1e-9  # perfect

def test_strand_normalization():
    COMP = str.maketrans("ACGTacgt","TGCAtgca")
    rc = lambda s: s.translate(COMP)[::-1]
    # a minus-strand allele (ref 'A' where + genome is 'T') normalizes to +strand ref 'T'
    base, ref, alt = "T", "A", "G"
    assert base == rc(ref)                 # detected as minus-strand
    assert (rc(ref), rc(alt)) == ("T","C")  # flipped to +strand

if __name__ == "__main__":
    test_metrics_sanity(); test_strand_normalization(); print("smoke tests pass")
