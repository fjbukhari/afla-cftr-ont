#!/usr/bin/env sh
# Usage: pad_merge_bed.sh <in.bed> <ref.fa.fai> <padding>
# Removes header lines, matches chromosome naming to the reference (adds or removes "chr"),
# drops contigs not in the reference, pads each region, clamps to chromosome ends, sorts and merges overlaps.
set -eu
bed=$1; fai=$2; pad=$3
awk -v OFS='\t' -v pad="$pad" '
  NR==FNR { len[$1]=$2; next }
  /^(#|track|browser)/ || NF<3 { next }
  {
    c=$1
    if (!(c in len)) { if (("chr" c) in len) c="chr" c; else if (sub(/^chr/,"",c) && (c in len)) {} else next }
    if (c=="M" && ("MT" in len)) c="MT"
    s=$2-pad; if (s<0) s=0
    e=$3+pad; if (e>len[c]) e=len[c]
    print c, s, e
  }' "$fai" "$bed" \
| sort -k1,1 -k2,2n \
| awk -v OFS='\t' '
  $1!=c || $2>e { if (c!="") print c,s,e; c=$1; s=$2; e=$3; next }
  $3>e { e=$3 }
  END { if (c!="") print c,s,e }'
