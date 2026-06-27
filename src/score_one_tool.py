"""Score all MFASS mutants with ONE tool -> scores_<tool>.csv.
Decoupled per tool (one tool's failure can't lose another's), resumable, chunked.
SpliceTransformer is sub-batched to fit GPU memory."""
import csv, os, time, argparse
from pyfaidx import Fasta
import numpy as np

import os
HG38 = os.environ.get("HG38", "ref/hg38.fa")
MFASS = os.environ.get("MFASS", "data/mfass/snv_data_clean.txt")
RES = os.environ.get("RES", "results")
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]
def normalize(x, fetch):
    """Return (+strand ref, +strand alt) or None. Auto-flips minus-strand alleles."""
    p0 = x["pos1"]-1
    try: base = fetch(x["chrom"], p0, p0+1)
    except Exception: return None
    if len(base)!=1: return None
    if base==x["ref"]: return x["ref"], x["alt"]
    if base==rc(x["ref"]): return rc(x["ref"]), rc(x["alt"])  # alleles on minus strand -> +strand
    return None  # true ref mismatch (coord/indel) -> skip

def load_mutants():
    rows = []
    with open(MFASS) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d.get("category") != "mutant": continue
            try: pos = int(d["snp_position_hg38_1based"])
            except (ValueError, KeyError, TypeError): continue
            rows.append({"id": d["id"], "chrom": d["chr"], "pos1": pos, "ref": d["ref_allele"],
                         "alt": d["alt_allele"], "strand": d["strand"]})
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=["spliceai","pangolin","splicetx"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=0)   # variants per write-flush
    ap.add_argument("--sub", type=int, default=16)    # splicetx variants per model call
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshards", type=int, default=1)
    args = ap.parse_args()
    chunk = args.chunk or {"spliceai":500,"pangolin":300,"splicetx":64}[args.tool]
    fa = Fasta(HG38)
    def fetch(c,s0,e0): return str(fa[c][s0:e0]).upper()
    rows = load_mutants()
    if args.limit: rows = rows[:args.limit]
    if args.nshards > 1: rows = rows[args.shard::args.nshards]
    out = f"{RES}/scores_{args.tool}_s{args.shard}of{args.nshards}.csv" if args.nshards > 1 else f"{RES}/scores_{args.tool}.csv"
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            for d in csv.DictReader(f): done.add(d["id"])
    todo = [x for x in rows if x["id"] not in done]
    print(f"[{args.tool}] total {len(rows)} done {len(done)} todo {len(todo)}", flush=True)
    new = not os.path.exists(out)
    fh = open(out,"a",newline=""); w = csv.DictWriter(fh, fieldnames=["id","score"])
    if new: w.writeheader(); fh.flush()
    t0 = time.time(); n = 0

    if args.tool == "spliceai":
        from proto_tools import run_spliceai_score, SpliceAIScoreInput, SpliceAIScoreConfig, SpliceAIVariant
        cfg = SpliceAIScoreConfig(reference_fasta=HG38, annotation="grch38", device="cuda")
        def score(ch):
            vs, keys = [], []
            for x in ch:
                nr = normalize(x, fetch)
                if nr is None: continue
                vs.append(SpliceAIVariant(chromosome=x["chrom"], position=x["pos1"], ref=nr[0], alt=nr[1])); keys.append(x["id"])
            r = {}
            if vs:
                o = run_spliceai_score(SpliceAIScoreInput(variants=vs), config=cfg)
                for k,res in zip(keys,o.results):
                    v = res.metrics.get("max_delta_score"); r[k] = v if isinstance(v,(int,float)) else 0.0
            return r
    elif args.tool == "pangolin":
        from proto_tools import run_pangolin_score_variants, PangolinScoreVariantsInput, PangolinScoreVariantsConfig, PangolinVariant
        cfg = PangolinScoreVariantsConfig(device="cuda"); F = 5050
        def score(ch):
            pv, keys = [], []
            for x in ch:
                nr = normalize(x, fetch)
                if nr is None: continue
                p0 = x["pos1"]-1; ws = p0-F; seq = fetch(x["chrom"], ws, p0+F+1)
                if len(seq) != 2*F+1 or seq[F] != nr[0]: continue
                pv.append(PangolinVariant(sequence=seq, variant_position=p0-ws, reference_bases=nr[0], alternate_bases=nr[1], strand=x["strand"])); keys.append(x["id"])
            r = {}
            if pv:
                o = run_pangolin_score_variants(PangolinScoreVariantsInput(variants=pv), config=cfg)
                for k,res in zip(keys,o.results):
                    mg = res.metrics.get("max_gain") or 0; ml = res.metrics.get("max_loss") or 0; r[k] = max(mg,-ml)
            return r
    else:
        from proto_tools import run_splice_transformer, SpliceTransformerInput, SpliceTransformerConfig
        cfg = SpliceTransformerConfig(device="cuda"); T, C = 1000, 4000
        def win(c,p0,al,strand):
            s = p0-(C+T//2); full = list(fetch(c,s,p0+(C+T//2))); 
            if al: full[p0-s] = al
            full = "".join(full)
            if strand == "-": full = rc(full)
            return full[:C], full[C:C+T], full[C+T:C+T+C]
        def score(ch):
            r = {}
            for j in range(0, len(ch), args.sub):
                sub = ch[j:j+args.sub]; L,Tt,R,keys = [],[],[],[]
                for x in sub:
                    nr = normalize(x, fetch)
                    if nr is None: continue
                    p0 = x["pos1"]-1
                    ok = True; seqs = []
                    for al in (None, nr[1]):
                        l,t,rr = win(x["chrom"],p0,al,x["strand"])
                        if len(l)!=C or len(t)!=T or len(rr)!=C: ok=False; break
                        seqs.append((l,t,rr))
                    if not ok: continue
                    for l,t,rr in seqs: L.append(l); Tt.append(t); R.append(rr)
                    keys.append(x["id"])
                if not keys: continue
                o = run_splice_transformer(SpliceTransformerInput(target_seqs=Tt,left_contexts=L,right_contexts=R), config=cfg)
                pred = np.array(o.prediction)
                for i,k in enumerate(keys):
                    rp,ap = pred[2*i],pred[2*i+1]
                    r[k] = float(max(np.abs(ap[:,1]-rp[:,1]).max(), np.abs(ap[:,2]-rp[:,2]).max()))
            return r

    for i in range(0, len(todo), chunk):
        ch = todo[i:i+chunk]
        try:
            sc = score(ch)
            for x in ch: w.writerow({"id": x["id"], "score": sc.get(x["id"], "")})
            fh.flush(); n += len(ch)
            el = time.time()-t0
            print(f"[{args.tool}] {n}/{len(todo)} {el:.0f}s {el/max(n,1):.3f}s/var eta {(len(todo)-n)*el/max(n,1)/60:.0f}min", flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); print(f"[{args.tool}] CHUNK {i} ERR {repr(e)[:160]}", flush=True)
    fh.close(); print(f"{args.tool.upper()}_DONE", flush=True)

if __name__ == "__main__": main()
