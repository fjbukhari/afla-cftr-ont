// CFTR-specific steps: poly-T/TG, exon copy number and the one-page report.
// The helper scripts need resource files passed as inputs: inside a container only bin/ and the task folder are visible.

// ---------------------------------------------------------------- Nanopore
process MINIMAP2_INDEX {
    label "big"
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
    input:
    path fa
    output:
    path "${fa.baseName}.map-ont.mmi"
    script:
    """
    minimap2 -x map-ont -t ${task.cpus} -d ${fa.baseName}.map-ont.mmi ${fa}
    """
}

// The minimap2 image has no samtools, so aligning and sorting are two steps.
process ALIGN_ONT {
    tag "${sample}"
    label "big"
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
    input:
    tuple val(sample), path(reads)
    path index
    output:
    tuple val(sample), path("${sample}.sam.gz")
    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:ONT' ${index} ${reads} | gzip -1 > ${sample}.sam.gz
    """
}

process SORT_ONT {
    tag "${sample}"
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    cpus 4
    publishDir { "${params.out_dir}/${sample}/alignment" }, mode: "copy"
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
    """
    samtools sort -@ ${task.cpus} -m 1G -T tmp_sort -o ${sample}.bam ${sam}
    samtools index -@ ${task.cpus} ${sample}.bam
    """
}

// bcftools with its Nanopore profile: default caller, small image, runs on any computer.
process CALL_BCFTOOLS_ONT {
    tag "${sample}"
    container "staphb/bcftools:1.23.1"
    cpus 4
    input:
    tuple val(sample), path(aln), path(idx), path(bed)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.raw.vcf.gz")
    script:
    """
    bcftools mpileup -X ${params.bcftools_profile} -f ${ref_name} -R ${bed} -d 8000 -L 8000 -a FORMAT/AD,FORMAT/DP \\
      --threads ${task.cpus} -Ou ${aln} \\
      | bcftools call -m -v -f GQ --threads ${task.cpus} -Oz -o ${sample}.raw.vcf.gz
    """
}

// Clair3 settings for amplicons: all candidate sites go through the full-alignment model and phasing.
process CALL_CLAIR3 {
    tag "${sample}"
    label "big"
    container "quay.io/biocontainers/clair3:1.2.0--py310h779eee5_0"
    input:
    tuple val(sample), path(aln), path(idx), path(bed)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.raw.vcf.gz")
    script:
    """
    MODELS=\$(dirname \$(which run_clair3.sh))/models
    [ -d "\$MODELS/${params.clair3_model}" ] || { echo "Clair3 model ${params.clair3_model} not found. Available: \$(ls \$MODELS)" >&2; exit 1; }
    run_clair3.sh --bam_fn=${aln} --ref_fn=${ref_name} --bed_fn=${bed} --threads=${task.cpus} --platform=ont \
      --model_path=\$MODELS/${params.clair3_model} --output=clair3_out --sample_name=${sample} \
      --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --include_all_ctgs --enable_variant_calling_at_sequence_head_and_tail
    mv clair3_out/merge_output.vcf.gz ${sample}.raw.vcf.gz
    rm -rf clair3_out
    """
}

// Read-based phasing: two variants in one read set can be placed on the same or on opposite chromosomes.
process WHATSHAP {
    tag "${sample}"
    container "quay.io/biocontainers/whatshap:2.8--py312h9b2995c_1"
    cpus 2
    publishDir { "${params.out_dir}/${sample}/variants" }, mode: "copy"
    input:
    tuple val(sample), path(vcf), path(tbi), path(aln), path(idx)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path("${sample}.phased.vcf.gz"), path("${sample}.phased.vcf.gz.tbi")
    script:
    """
    whatshap phase --ignore-read-groups --reference ${ref_name} -o ${sample}.phased.vcf.gz ${vcf} ${aln}
    python -c "import pysam; pysam.tabix_index('${sample}.phased.vcf.gz', preset='vcf', force=True, keep_original=True)"
    """
}

// ---------------------------------------------------------------- both platforms
// Reads over the poly-T/TG tract (c.1210-34), as SAM text for the caller in the next step.
process CFTR_POLYT_READS {
    tag "${sample}"
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    input:
    tuple val(sample), path(aln), path(idx)
    output:
    tuple val(sample), path("${sample}.polyt.sam")
    script:
    """
    samtools view -F 0x904 ${aln} chr7:117548590-117548650 > ${sample}.polyt.sam
    """
}

process CFTR_POLYT {
    tag "${sample}"
    container "quay.io/biocontainers/python:3.12"
    publishDir { "${params.out_dir}/${sample}/cftr" }, mode: "copy"
    input:
    tuple val(sample), path(sam)
    path ref_files
    val ref_name
    path exons
    output:
    tuple val(sample), path("${sample}.polyt.json")
    script:
    """
    afla_cftr_polyt.py --sam ${sam} --fasta ${ref_name} --exons ${exons} --sample ${sample} \\
      --platform ${params.platform} --out ${sample}.polyt.json
    """
}

// Exon copy number: every sample of the run compared with the others.
process CFTR_CNV {
    container "quay.io/biocontainers/python:3.12"
    publishDir "${params.out_dir}/cftr", mode: "copy"
    input:
    path regions
    path exons
    path deletions
    output:
    path "cftr_cnv.json", emit: json
    path "cftr_cnv.tsv"
    script:
    def fs = regions instanceof List ? regions : [regions]
    def args = fs.collect { f -> "${f.name.replace('.targets.regions.bed.gz', '')}=${f.name}" }.join(' ')
    """
    afla_cftr_cnv.py --depth ${args} --deletions ${deletions} --min-depth ${params.cnv_min_depth} \\
      --out cftr_cnv.json --tsv cftr_cnv.tsv
    """
}

process CFTR_REPORT {
    tag "${sample}"
    container "quay.io/biocontainers/python:3.12"
    cache false     // seconds to run; always rebuilt so report changes apply
    publishDir { "${params.out_dir}/${sample}" }, mode: "copy"
    input:
    tuple val(sample), path(vcf), path(thresholds), path(regions), path(polyt), path(cnv)
    path ref_files
    val ref_name
    path exons
    path cftr2
    val extra_args
    output:
    path "${sample}.cftr.html"
    path "${sample}.cftr.json"
    script:
    """
    afla_cftr.py --sample ${sample} --vcf ${vcf} --polyt ${polyt} --cnv ${cnv} --thresholds ${thresholds} \\
      --regions ${regions} --fasta ${ref_name} --cftr2 ${cftr2} --exons ${exons} --platform ${params.platform} \\
      --min-depth ${params.report_min_depth} --out ${sample}.cftr.html --json ${sample}.cftr.json ${extra_args}
    """
}
