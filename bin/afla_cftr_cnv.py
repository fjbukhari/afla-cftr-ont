#!/usr/bin/env python3
"""CFTR exon copy number from amplicon depth, compared with the other samples of the same run
(education and research use only, not for clinical use).

  afla_cftr_cnv.py --depth S1=S1.regions.bed.gz S2=S2.regions.bed.gz ... --out cnv.json [--tsv cnv.tsv]

Depth files: mosdepth --by <targets.bed> output (<prefix>.regions.bed.gz: chrom, start, end, name, mean depth).
For each sample, each target's depth is divided by the sample's median target depth; that is compared with the
median of the same target in the other samples (leave-one-out). A ratio near 0.5 suggests a heterozygous
deletion, near 0 a homozygous deletion, near 1.5 a duplication. Only exon targets are called; at least 4 samples
are needed (3 to compare with). Calls on a single exon are labelled for confirmation (MLPA or another method):
a variant under a PCR primer can make one amplicon drop out and look like a deletion.
Standard library only.
"""
import argparse
import gzip
import json
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import afla_cftr_hgvs as hg  # noqa: E402

DEL_HOM, DEL_HET, DUP = 0.2, 0.7, 1.35
MIN_SAMPLES = 4
EXON_RE = re.compile(r"CFTR_exon(\d+)")


def read_depth(path):
    out = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                out[p[3]] = float(p[4])
    return out


def state(r):
    if r is None:
        return "not assessed"
    if r < DEL_HOM:
        return "homozygous deletion?"
    if r < DEL_HET:
        return "heterozygous deletion?"
    if r > DUP:
        return "duplication?"
    return "normal"


def known_deletions(path=hg.RESOURCES / "CFTR_large_deletions.tsv"):
    out = []
    if Path(path).exists():
        for r in hg.read_tsv(path):
            a, _, b = r["exons"].partition("-")
            out.append((int(a), int(b or a), r["legacy"], r["cdna"]))
    return out


def analyse(depths, min_depth=30, deletions=hg.RESOURCES / "CFTR_large_deletions.tsv"):
    """depths: {sample: {target: mean depth}}. Returns {sample: result dict}."""
    samples = sorted(depths)
    names = sorted({t for d in depths.values() for t in d if EXON_RE.match(t)},
                   key=lambda t: int(EXON_RE.match(t).group(1)))
    norm = {}
    for s in samples:
        vals = [depths[s].get(t, 0.0) for t in names]
        med = statistics.median(vals) if vals else 0
        norm[s] = {t: (depths[s].get(t, 0.0) / med if med > 0 else None) for t in names}
    res = {}
    known = known_deletions(deletions)
    for s in samples:
        others = [o for o in samples if o != s]
        rows = []
        for t in names:
            exon = int(EXON_RE.match(t).group(1))
            ref = [norm[o][t] for o in others if norm[o][t] is not None]
            ratio = None
            note = ""
            if len(others) + 1 < MIN_SAMPLES:
                note = f"needs at least {MIN_SAMPLES} samples in the run"
            elif depths[s].get(t, 0.0) < min_depth and statistics.median(
                    [depths[o].get(t, 0.0) for o in others]) < min_depth:
                note = f"depth below {min_depth} in most samples"
            elif ref and statistics.median(ref) > 0 and norm[s][t] is not None:
                ratio = norm[s][t] / statistics.median(ref)
                spread = statistics.pstdev(ref) / statistics.mean(ref) if len(ref) > 1 and statistics.mean(ref) else 0
                if spread > 0.3:
                    note = f"variable in this run (CV {spread:.0%}): less reliable"
            rows.append({"target": t, "exon": exon, "legacy": hg.LEGACY_EXON[exon], "depth": round(depths[s].get(t, 0.0), 1),
                         "ratio": round(ratio, 2) if ratio is not None else None, "state": state(ratio), "note": note})
        # merge neighbouring abnormal exons into segments
        segs, cur = [], None
        for r in rows:
            st = r["state"]
            if st in ("normal", "not assessed"):
                cur = None
                continue
            kind = "deletion" if "deletion" in st else "duplication"
            if cur and cur["kind"] == kind and r["exon"] == cur["exons"][-1] + 1:
                cur["exons"].append(r["exon"])
                cur["ratios"].append(r["ratio"])
            else:
                cur = {"kind": kind, "exons": [r["exon"]], "ratios": [r["ratio"]]}
                segs.append(cur)
        calls = []
        for g in segs:
            a, b = g["exons"][0], g["exons"][-1]
            mean_r = sum(g["ratios"]) / len(g["ratios"])
            zyg = "homozygous" if g["kind"] == "deletion" and mean_r < DEL_HOM else (
                "heterozygous" if g["kind"] == "deletion" else "")
            label = f"exon{'s' if b > a else ''} {a}" + (f"-{b}" if b > a else "")
            leg = f"legacy {hg.LEGACY_EXON[a]}" + (f"-{hg.LEGACY_EXON[b]}" if b > a else "")
            match = [k for k in known if k[0] == a and k[1] == b and g["kind"] == "deletion"]
            calls.append({
                "kind": g["kind"], "zygosity": zyg, "exons": label, "legacy_exons": leg,
                "mean_ratio": round(mean_r, 2), "known": (f"{match[0][2]} ({match[0][3]})" if match else ""),
                "confirm": "single exon: confirm by MLPA or another method (allele drop-out under a primer can "
                           "look like a deletion)" if a == b else "confirm by MLPA or another method",
            })
        assessed = sum(1 for r in rows if r["ratio"] is not None)
        res[s] = {"sample": s, "exons_assessed": assessed, "exons_total": len(rows), "calls": calls, "exons": rows,
                  "status": ("not assessed: " + rows[0]["note"]) if rows and assessed == 0 else "assessed"}
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth", nargs="+", required=True, help="sample=path/to/sample.regions.bed.gz")
    ap.add_argument("--min-depth", type=int, default=30)
    ap.add_argument("--deletions", default=str(hg.RESOURCES / "CFTR_large_deletions.tsv"),
                    help="known multi-exon deletions (inside a container, pass this file as an input)")
    ap.add_argument("--out", required=True, help="JSON with one entry per sample")
    ap.add_argument("--tsv", help="also a flat table")
    a = ap.parse_args(argv)
    depths = {}
    for item in a.depth:
        s, _, p = item.partition("=")
        if not p:
            p, s = s, Path(s).name.split(".")[0]
        depths[s] = read_depth(p)
    res = analyse(depths, a.min_depth, a.deletions)
    Path(a.out).write_text(json.dumps(res, indent=1))
    if a.tsv:
        with open(a.tsv, "w") as f:
            f.write("sample\ttarget\texon\tlegacy_exon\tdepth\tratio\tstate\tnote\n")
            for s, r in res.items():
                for e in r["exons"]:
                    f.write(f"{s}\t{e['target']}\t{e['exon']}\t{e['legacy']}\t{e['depth']}\t"
                            f"{'' if e['ratio'] is None else e['ratio']}\t{e['state']}\t{e['note']}\n")
    for s, r in res.items():
        calls = "; ".join(f"{c['zygosity']} {c['kind']} {c['exons']} ({c['legacy_exons']}) ratio {c['mean_ratio']}"
                          + (f" = {c['known']}" if c["known"] else "") for c in r["calls"]) or "none"
        print(f"{s}: {r['status']}; calls: {calls}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
