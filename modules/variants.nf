// Variant calling, normalisation/quality filter and annotation.

// GATK HaplotypeCaller (CPU): the widely used standard germline caller.
// Amplicon mode: by default HaplotypeCaller keeps only 50 reads per start position. In amplicon data
// thousands of clean reads share one start, so it would discard ~97% of them while keeping every read that
// starts elsewhere (mostly error-rich ones): noise then looks like 5-20% variants. Down-sampling is switched
// off, and soft-clipped (adapter/probe) bases are kept out of local assembly.
process CALL_GATK {
    tag "${sample}"
    label "big"
    container "broadinstitute/gatk:4.6.2.0"
    input:
    tuple val(sample), path(aln), path(idx), path(bed)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.raw.vcf.gz"), emit: vcf
    script:
    def amp = (params.amplicon.toString() == 'true') ? "--max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases" : ""
    """
    gatk --java-options "-Xmx16g" HaplotypeCaller -R ${ref_name} -I ${aln} -L ${bed} ${amp} \
      --native-pair-hmm-threads ${task.cpus} -O ${sample}.raw.vcf.gz
    """
}

// DeepVariant (Google): deep-learning caller; uses the GPU image when use_gpu is on.
process CALL_DEEPVARIANT {
    tag "${sample}"
    label "big"
    label "gpu"
    container "google/deepvariant:${params.deepvariant_version}${(params.use_gpu.toString() == 'true') ? '-gpu' : ''}"
    input:
    tuple val(sample), path(aln), path(idx), path(bed)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.raw.vcf.gz"), emit: vcf
    script:
    """
    run_deepvariant --model_type=${params.deepvariant_model} --ref=${ref_name} --reads=${aln} \
      --regions=${bed} --output_vcf=${sample}.raw.vcf.gz --num_shards=${task.cpus} \
      --intermediate_results_dir=dv_tmp
    rm -rf dv_tmp
    """
}

// Splits multi-allelic sites, left-aligns indels, labels low-quality calls (they are kept, marked LowQual).
process NORMALISE_FILTER {
    tag "${sample}"
    container "staphb/bcftools:1.23.1"
    publishDir { "${params.out_dir}/${sample}/variants" }, mode: "copy"
    input:
    tuple val(sample), path(vcf)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.tbi"), emit: vcf
    script:
    """
    bcftools norm -m -both -f ${ref_name} --check-ref s ${vcf} -Ou \
      | bcftools filter -m + -s LowQual \
          -e 'QUAL<${params.min_qual} || FMT/DP<${params.min_dp} || FMT/GQ<${params.min_gq}' \
          -Oz -o ${sample}.vcf.gz
    bcftools index -t ${sample}.vcf.gz
    """
}

// VCF given as input: sort/compress/index it so it enters the same path as called variants.
process PREPARE_INPUT_VCF {
    tag "${sample}"
    container "staphb/bcftools:1.23.1"
    input:
    tuple val(sample), path(vcf)
    output:
    tuple val(sample), path("${sample}.input.vcf.gz"), path("${sample}.input.vcf.gz.tbi"), emit: vcf
    script:
    """
    bcftools sort ${vcf} -Oz -o ${sample}.input.vcf.gz
    bcftools index -t ${sample}.input.vcf.gz
    """
}

// Ensembl VEP, offline cache. The Ensembl cache names chromosomes 1,2..MT; the VCF is renamed for VEP
// and renamed back afterwards, so the output keeps the reference's own naming.
process ANNOTATE_VEP {
    tag "${sample}"
    label "big"
    container "ensemblorg/ensembl-vep:release_115.2"
    publishDir { "${params.out_dir}/${sample}/variants" }, mode: "copy"
    input:
    tuple val(sample), path(vcf), path(tbi)
    path vep_cache
    path extra_files
    val extra_args
    output:
    tuple val(sample), path("${sample}.annotated.vcf.gz"), emit: vcf
    path "${sample}.vep_summary.html"
    script:
    """
    # rename chr1 -> 1 (and chrM -> MT) only if this VCF uses "chr" names
    CHR=\$(zcat ${vcf} | awk '/^##contig=<ID=chr/ {f=1} END {print f+0}')
    zcat ${vcf} | awk -v chr=\$CHR 'BEGIN{OFS="\\t"}
        !chr { print; next }
        /^##contig=<ID=chr/ { sub(/ID=chrM,/,"ID=MT,"); sub(/ID=chr/,"ID="); print; next }
        /^#/ { print; next }
        { if (\$1=="chrM") \$1="MT"; else sub(/^chr/,"",\$1); print }' > for_vep.vcf
    vep --input_file for_vep.vcf --format vcf --output_file vep_out.vcf --vcf --force_overwrite \
      --stats_file ${sample}.vep_summary.html \
      --offline --cache --dir_cache ${vep_cache} --cache_version ${params.vep_cache_version} \
      --assembly GRCh38 --species homo_sapiens --fork ${task.cpus} --buffer_size 5000 \
      --pick --mane --symbol --canonical --biotype --numbers --variant_class --sift b --polyphen b \
      --af --af_gnomade --af_gnomadg --max_af --check_existing --pubmed --gene_phenotype \
      ${extra_args}
    # VEP writes no file when there are no variants
    [ -e vep_out.vcf ] || cp for_vep.vcf vep_out.vcf
    awk -v chr=\$CHR 'BEGIN{OFS="\\t"}
        !chr { print; next }
        /^##contig=<ID=MT,/ { sub(/ID=MT,/,"ID=chrM,"); print; next }
        /^##contig=<ID=/ { sub(/ID=/,"ID=chr"); print; next }
        /^#/ { print; next }
        { if (\$1=="MT") \$1="chrM"; else \$1="chr" \$1; print }' vep_out.vcf | bgzip > ${sample}.annotated.vcf.gz
    """
}
