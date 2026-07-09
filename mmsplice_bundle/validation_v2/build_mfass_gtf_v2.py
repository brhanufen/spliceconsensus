"""VALIDATION v2: build a GTF where each MFASS exon is an INTERNAL (cassette) exon,
so MMSplice computes BOTH acceptor and donor modules with real genomic flanking intron
sequence. The v1 GTF made each exon a single-exon transcript, which can make MMSplice
skip its splice modules (first exon has no acceptor, last has no donor, lone exon has
neither) -> degenerate near-random scores. This is the leading suspect for the weak v1 result.

Construction (standard cassette trick, fully honest):
  For each MFASS exon E on its strand, emit one transcript with three exons:
     [dummy upstream exon] -- intron -- [E, the real MFASS exon] -- intron -- [dummy downstream exon]
  The dummy exons are placed FLANK bp into the real genomic sequence beyond E's boundaries,
  so the introns bounding E are REAL genomic sequence (the same context the other tools use).
  Only E carries scored variants; the dummies exist solely to define E's two splice sites.

Outputs mfass_exons_v2.gtf. VCF is unchanged (reuse mfass.sorted.vcf.gz from v1).
Run on HCC:
  python build_mfass_gtf_v2.py --mfass data/mfass/snv_data_clean.txt --out_gtf mfass_exons_v2.gtf --flank 300 --dummy 50
"""
import csv, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mfass", required=True)
    ap.add_argument("--out_gtf", default="mfass_exons_v2.gtf")
    ap.add_argument("--flank", type=int, default=300, help="bp from E boundary to the dummy exon (intron length)")
    ap.add_argument("--dummy", type=int, default=50, help="length of each dummy flanking exon")
    a = ap.parse_args()

    exons = {}  # exon_id -> (chrom, start0, end0, strand)
    with open(a.mfass) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d.get("category") != "mutant": continue
            try:
                s0=int(d["start_hg38_0based"]); i1=int(d["intron1_len"])
                exl=int(d["exon_len"]); i2=int(d["intron2_len"])
            except (KeyError,ValueError,TypeError): continue
            off = i1 if d["strand"]=="+" else i2
            ex_start0 = s0+off; ex_end0 = s0+off+exl
            exons.setdefault(d["id"].rsplit("_",1)[0], (d["chr"], ex_start0, ex_end0, d["strand"]))

    F, DUM = a.flank, a.dummy
    with open(a.out_gtf, "w") as g:
        for exon_id,(chrom,s0,e0,strand) in exons.items():
            # dummy exons flank E by F bp of (real genomic) intron on each side
            up_s, up_e   = s0 - F - DUM, s0 - F          # upstream dummy (0-based)
            dn_s, dn_e   = e0 + F,       e0 + F + DUM     # downstream dummy
            tid=f"{exon_id}_cass"; gid=f"{exon_id}_g"
            g.write(f'{chrom}\tMFASS\tgene\t{up_s+1}\t{dn_e}\t.\t{strand}\t.\tgene_id "{gid}"; gene_name "{gid}";\n')
            g.write(f'{chrom}\tMFASS\ttranscript\t{up_s+1}\t{dn_e}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}";\n')
            # three exons; MFASS exon E is the internal one -> both splice sites defined
            for (xs,xe,tag) in [(up_s,up_e,"up"),(s0,e0,"E"),(dn_s,dn_e,"dn")]:
                g.write(f'{chrom}\tMFASS\texon\t{xs+1}\t{xe}\t.\t{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; exon_id "{exon_id}_{tag}";\n')
    print(f"wrote {a.out_gtf}: {len(exons)} cassette transcripts (MFASS exon internal, flank={F}bp, dummy={DUM}bp)")
    print("NOTE: variants sit only in the internal exon E and its immediate flanks; dummy exons are unscored scaffolding.")

if __name__ == "__main__": main()
