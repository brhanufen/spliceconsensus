"""Build a MMSplice-ready VCF + exon GTF for the MFASS variants, directly from the
MFASS table's hg38 coordinates. No external GENCODE download needed: the exon
annotation is reconstructed from MFASS's own construct coordinates, guaranteeing
that every variant's exon is present and correctly bounded.

Strand handling reuses the exact logic from score_one_tool.py: each variant's ref
allele is checked against the real hg38 FASTA base and minus-strand alleles are
flipped to the plus strand (the paper's r=0.997 coordinate concordance comes from this).

Run on HCC (needs hg38.fa + pyfaidx):
    python build_mfass_vcf.py --mfass data/mfass/snv_data_clean.txt --hg38 ref/hg38.fa \
        --out_vcf mfass.vcf --out_gtf mfass_exons.gtf
Outputs a bgzip/sort-ready VCF (sort+bgzip+tabix commands printed at the end) and a GTF.
"""
import csv, argparse, sys
from pyfaidx import Fasta

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mfass", required=True)
    ap.add_argument("--hg38", required=True)
    ap.add_argument("--out_vcf", default="mfass.vcf")
    ap.add_argument("--out_gtf", default="mfass_exons.gtf")
    a = ap.parse_args()
    fa = Fasta(a.hg38)
    def base(chrom, pos1):  # 1-based -> base
        return str(fa[chrom][pos1-1:pos1]).upper()

    variants = []      # (chrom, pos1, id, ref+, alt+, exon_id)
    exons = {}         # exon_id -> (chrom, start0, end0, strand)
    n=skipped=flipped=0
    with open(a.mfass) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d.get("category") != "mutant": continue
            n+=1
            try:
                pos1 = int(d["snp_position_hg38_1based"]); chrom=d["chr"]
                s0 = int(d["start_hg38_0based"]); e0=int(d["end_hg38_0based"])
                i1=int(d["intron1_len"]); exl=int(d["exon_len"]); i2=int(d["intron2_len"])
            except (KeyError,ValueError,TypeError):
                skipped+=1; continue
            ref,alt,strand = d["ref_allele"], d["alt_allele"], d["strand"]
            b = base(chrom,pos1)
            if b==ref: rp,ap_=ref,alt
            elif b==rc(ref): rp,ap_=rc(ref),rc(alt); flipped+=1  # minus-strand alleles -> +strand
            else: skipped+=1; continue                          # true mismatch/indel -> skip
            # exon genomic coords: construct is plus-oriented; exon offset from the plus-strand
            # 5' end is intron1_len on +strand genes, intron2_len on -strand genes.
            off = i1 if strand=="+" else i2
            ex_start0 = s0+off; ex_end0 = s0+off+exl
            exon_id = d["id"].rsplit("_",1)[0]
            exons.setdefault(exon_id,(chrom,ex_start0,ex_end0,strand))
            variants.append((chrom,pos1,d["id"],rp,ap_,exon_id))

    # write VCF (unsorted; sort/bgzip/tabix afterward)
    with open(a.out_vcf,"w") as v:
        v.write("##fileformat=VCFv4.2\n")
        for c in sorted({x[0] for x in variants}):
            v.write(f"##contig=<ID={c}>\n")
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom,pos1,vid,ref,alt,_ in variants:
            v.write(f"{chrom}\t{pos1}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")

    # write GTF (one gene+transcript+exon per MFASS exon; MMSplice needs exon+transcript+gene rows)
    with open(a.out_gtf,"w") as g:
        for exon_id,(chrom,s0,e0,strand) in exons.items():
            attr=f'gene_id "{exon_id}"; transcript_id "{exon_id}"; exon_id "{exon_id}"; gene_name "{exon_id}";'
            for feat in ("gene","transcript","exon"):
                g.write(f"{chrom}\tMFASS\t{feat}\t{s0+1}\t{e0}\t.\t{strand}\t.\t{attr}\n")

    print(f"mutants read: {n}  written: {len(variants)}  skipped(mismatch/indel): {skipped}  minus-strand flipped: {flipped}")
    print(f"distinct exons in GTF: {len(exons)}")
    print("\nNext (on HCC):")
    print(f"  bcftools sort {a.out_vcf} -Oz -o mfass.sorted.vcf.gz && tabix -p vcf mfass.sorted.vcf.gz")
    print(f"  (or: bgzip {a.out_vcf} && tabix -p vcf {a.out_vcf}.gz  after sort)")

if __name__=="__main__": main()
