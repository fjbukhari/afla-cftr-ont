// Reference genome and target-region preparation.

// Makes the .fai (position index) and .dict (chromosome list) when they are not next to the FASTA.
process PREPARE_REF {
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    input:
    path ref
    output:
    path "${ref}.fai", emit: fai
    path "*.dict", emit: dict
    script:
    """
    samtools faidx ${ref}
    samtools dict ${ref} -o ${ref.baseName}.dict
    """
}

// One-time bwa index of the reference (about 1 hour for the human genome).
process BWA_INDEX {
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    publishDir "${params.out_dir}/reference", mode: "copy"
    input:
    path ref
    output:
    path "${ref}.{amb,ann,bwt,pac,sa}"
    script:
    """
    bwa index ${ref}
    """
}

// Cleans a user BED (headers, chr naming), pads each region and merges overlaps.
process PREPARE_BED {
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    publishDir "${params.out_dir}/targets", mode: "copy"
    input:
    path bed
    path fai
    output:
    path "targets.padded.bed", emit: padded
    path "targets.named.bed", emit: named
    script:
    // targets.named.bed: the original target regions, unpadded, with their names (usually the gene) kept, so
    // coverage can be reported per gene for the clinical report. Contig names are matched to the reference.
    """
    pad_merge_bed.sh ${bed} ${fai} ${params.bed_padding} > targets.padded.bed
    awk -v OFS='\\t' 'NR==FNR { len[\$1]=1; next }
      /^(#|track|browser)/ || NF<3 { next }
      { c=\$1; if (!(c in len)) { if (("chr" c) in len) c="chr" c; else if (sub(/^chr/,"",c) && (c in len)) {} else next }
        print c, \$2, \$3, (NF>=4 && \$4!="" ? \$4 : c ":" \$2 "-" \$3) }' ${fai} ${bed} \\
      | sort -k1,1 -k2,2n > targets.named.bed
    """
}

// Without a BED file: use the regions the reads actually cover (at least N reads deep).
process AUTO_TARGETS {
    tag "${sample}"
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    publishDir "${params.out_dir}/targets", mode: "copy"
    input:
    tuple val(sample), path(aln), path(idx)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.auto_targets.bed")
    script:
    """
    samtools depth --reference ${ref_name} -Q 20 -q 20 ${aln} \
      | awk -v OFS='\\t' -v min=${params.auto_targets_min_depth} '
          \$3>=min { if(\$1==c && \$2==e+1){e=\$2} else { if(c!="") print c,s-1,e; c=\$1; s=\$2; e=\$2 } }
          END { if(c!="") print c,s-1,e }' > covered.bed
    pad_merge_bed.sh covered.bed ${ref_name}.fai ${params.bed_padding} > ${sample}.auto_targets.bed
    """
}
