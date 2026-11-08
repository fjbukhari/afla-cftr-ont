#!/usr/bin/env python3
"""CFTR poly-TG / poly-T genotype, c.1210-34TG[n]T[n] (legacy IVS8 TG(m)T(n)); education and research use only,
not for clinical use.

For every read that spans the whole tract, the sequence between two fixed reference flanks (16 bases each side,
taken from the reference FASTA) is read as (TG)m(T)n. The two most common (m, n) are the alleles.

  samtools view -F 0x904 sample.bam chr7:117548590-117548650 | \\
      afla_cftr_polyt.py --sam - --fasta hg38.fa --sample S1 --platform illumina --out S1.polyt.json

Why it matters (teaching summary, see the report for wording): 5T lowers correct splicing of exon 10 (legacy
exon 9) and is a variant of varying clinical consequence, more so with TG12 or TG13; 7T and 9T are normal.
With R117H, 5T in cis is associated with CF and 7T with milder or no disease; cis/trans needs parental
samples because the amplicons are too short to join R117H (exon 4) and IVS8 in one read.

Nanopore reads miscount long homopolymers (T runs and, less often, TG repeats); ONT results are labelled lower
confidence and the full count distribution is reported so it can be judged.
Standard library only.
"""
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import afla_cftr_hgvs as hg  # noqa: E402

try:
    import edlib
except ImportError:
    edlib = None

TG_C, T_END_C = "1210-34", "1210-6"  # first base of the TG repeat, last base of the T run (GRCh38: TG11 T7)
FLANK = 16
TRACT_RE = re.compile(r"^((?:TG)*)(T*)$")


def locus(tx, fasta):
    a, b = tx.c_pos_to_g(TG_C), tx.c_pos_to_g(T_END_C)
    left = fasta.fetch(hg.CHROM, a - FLANK, a - 1)
    right = fasta.fetch(hg.CHROM, b + 1, b + FLANK)
    ref_tract = fasta.fetch(hg.CHROM, a, b)
    return {"start": a, "end": b, "left": left, "right": right, "ref": ref_tract}


def _find(q, s, k):
    if edlib is not None:
        r = edlib.align(q, s, mode="HW", task="locations", k=k)
        if r["editDistance"] < 0:
            return None
        start, end = r["locations"][0]  # edlib gives no start when the match reaches the very beginning of s
        return (0 if start is None else start), end
    i = s.find(q)
    return (i, i + len(q) - 1) if i >= 0 else None


def read_tract(seq, loc, k=2):
    """Returns ('TG11T7' parts) as (m, n, exact) or None if the read does not span the tract."""
    rc = seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    for s in (seq, rc):
        L = _find(loc["left"], s, k)
        if not L:
            continue
        R = _find(loc["right"], s[L[1] + 1:], k)
        if not R:
            continue
        mid = s[L[1] + 1: L[1] + 1 + R[0]]
        m = TRACT_RE.match(mid)
        if m:
            return len(m.group(1)) // 2, len(m.group(2)), True
        # imperfect tract (sequencing error inside): count the TG units at the start and the T run at the end
        tg = re.match(r"^((?:TG)*)", mid).group(1)
        t = re.search(r"(T*)$", mid).group(1)
        return len(tg) // 2, len(t), False
    return None


def call(counts, min_reads=10, min_frac=0.2):
    """counts: Counter {(m, n): reads}. Returns dict with genotype and support."""
    n = sum(counts.values())
    res = {"spanning_reads": n, "alleles": [], "genotype": "", "status": ""}
    if n < min_reads:
        res["status"] = f"not called: only {n} reads span the tract (need {min_reads})"
        return res
    top = counts.most_common()
    a1, c1 = top[0]
    second = next(((a, c) for a, c in top[1:] if c >= min_frac * n), None)
    if second:
        alleles = sorted([a1, second[0]], key=lambda x: (x[1], x[0]))
        res["alleles"] = [{"tg": a[0], "t": a[1], "reads": counts[a]} for a in alleles]
        res["status"] = "heterozygous"
    else:
        res["alleles"] = [{"tg": a1[0], "t": a1[1], "reads": c1}] * 2
        res["status"] = "homozygous"
    res["genotype"] = "/".join(f"TG{a['tg']}-{a['t']}T" for a in res["alleles"])
    res["poly_t"] = "/".join(f"{a['t']}T" for a in res["alleles"])
    res["poly_tg"] = "/".join(f"TG{a['tg']}" for a in res["alleles"])
    res["hgvs"] = " / ".join(f"c.1210-34TG[{a['tg']}]T[{a['t']}]" for a in res["alleles"])
    return res


def interpret(res, platform):
    notes = []
    if not res.get("alleles"):
        return notes
    ts = [a["t"] for a in res["alleles"]]
    if 5 in ts:
        tg5 = [a["tg"] for a in res["alleles"] if a["t"] == 5]
        risk = " with TG12 or TG13 (higher risk of disease)" if any(t >= 12 for t in tg5) else ""
        notes.append(f"5T present{risk}: 5T reduces correct splicing of exon 10 (legacy exon 9). It is a variant "
                     "of varying clinical consequence (CFTR-related disorders such as CBAVD; CF only with another "
                     "CF-causing variant in trans, or R117H in cis).")
    if any(t not in (5, 7, 9) for t in ts):
        notes.append("Unusual T count: check the read distribution" +
                     (" (Nanopore often miscounts T runs)" if platform == "ont" else "") + ".")
    if platform == "ont":
        notes.append("Nanopore result: lower confidence for repeat counts; confirm 5T by an orthogonal method.")
    return notes


def iter_sam(handle):
    for line in handle:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) > 9 and not int(f[1]) & 0x904:
            yield f[9]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sam", help="SAM records (from samtools view), '-' for stdin")
    ap.add_argument("--fastq", help="or reads as FASTQ (.gz)")
    ap.add_argument("--fasta", required=True, help="GRCh38 reference FASTA (with .fai)")
    ap.add_argument("--exons", default=str(hg.DEFAULT_EXONS))
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--platform", choices=["illumina", "ont"], default="illumina")
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--out", help="JSON output (default: print)")
    a = ap.parse_args(argv)
    tx = hg.Transcript.load(a.exons)
    loc = locus(tx, hg.Fasta(a.fasta))
    m = TRACT_RE.match(loc["ref"])
    if not m:
        sys.exit(f"reference tract {loc['ref']} is not (TG)n(T)n: wrong reference or exon table?")
    if a.sam:
        seqs = iter_sam(sys.stdin if a.sam == "-" else open(a.sam))
    elif a.fastq:
        import gzip
        op = gzip.open if a.fastq.endswith(".gz") else open
        seqs = (line.strip() for i, line in enumerate(op(a.fastq, "rt")) if i % 4 == 1)
    else:
        sys.exit("give --sam or --fastq")
    counts, imperfect = Counter(), 0
    for s in seqs:
        r = read_tract(s, loc, k=1 if a.platform == "illumina" else 3)
        if r:
            counts[(r[0], r[1])] += 1
            imperfect += not r[2]
    res = call(counts, min_reads=a.min_reads)
    res.update(sample=a.sample, platform=a.platform, imperfect_reads=imperfect,
               reference=f"TG{len(m.group(1)) // 2}-{len(m.group(2))}T",
               locus=f"{hg.CHROM}:{loc['start']}-{loc['end']}",
               distribution=[{"tg": k[0], "t": k[1], "reads": v} for k, v in counts.most_common(12)],
               notes=interpret(res, a.platform))
    if a.platform == "ont" and res.get("alleles"):
        res["confidence"] = "lower (Nanopore)"
    out = json.dumps(res, indent=1)
    if a.out:
        Path(a.out).write_text(out)
    else:
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
