// Read clean-up, alignment and coverage.

// fastp: removes adapters and low-quality ends, writes a QC report. Several lanes are joined first.
// Besides auto-detection it is given a list of known Illumina adapters: auto-detection misses some
// (e.g. the TruSight Amplicon read-1 tail), and untrimmed adapter creates false indels.
process FASTP {
    tag "${sample}"
    container "staphb/fastp:1.3.3"
    cpus 8
    publishDir { "${params.out_dir}/${sample}/qc" }, mode: "copy", pattern: "*.{html,json}"
    input:
    tuple val(sample), path(r1), path(r2)
    path adapters
    output:
    tuple val(sample), path("${sample}.R1.trim.fq.gz"), path("${sample}.R2.trim.fq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html"
    script:
    def r1s = r1 instanceof List ? r1 : [r1]
    def r2s = r2 instanceof List ? r2 : [r2]
    def in1 = r1s.size() > 1 ? "lanes_R1.fq.gz" : r1s[0]
    def in2 = r2s.size() > 1 ? "lanes_R2.fq.gz" : r2s[0]
    def join = r1s.size() > 1 ? "cat ${r1s.join(' ')} > lanes_R1.fq.gz; cat ${r2s.join(' ')} > lanes_R2.fq.gz" : ""
    """
    ${join}
    fastp -i ${in1} -I ${in2} \
      -o ${sample}.R1.trim.fq.gz -O ${sample}.R2.trim.fq.gz \
      --detect_adapter_for_pe --adapter_fasta ${adapters} -w ${task.cpus} \
      -j ${sample}.fastp.json -h ${sample}.fastp.html -R "${sample} fastp report"
    """
}

// bwa mem -> (duplicate marking) -> sorted, indexed CRAM or BAM. Streamed: no large temporary SAM file.
// Amplicon panels with probe/primer positions: the probe bases (which come from the oligo, not the patient)
// are soft-clipped with samtools ampliconclip, and the callers are told to ignore soft-clipped bases.
process ALIGN {
    tag "${sample}"
    label "big"
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    publishDir { "${params.out_dir}/${sample}/alignment" }, mode: "copy"
    input:
    tuple val(sample), path(r1), path(r2)
    path ref_files
    val ref_name
    path primers
    output:
    tuple val(sample), path("${sample}.${params.output_format}"), path("${sample}.${params.output_format}.{crai,bai}"), emit: aln
    path "${sample}.flagstat.txt", emit: flagstat
    path "${sample}.markdup.txt", optional: true, emit: markdup
    path "${sample}.primer_clip.txt", optional: true, emit: clipstats
    script:
    def r1s = r1 instanceof List ? r1 : [r1]
    def r2s = r2 instanceof List ? r2 : [r2]
    def fmt = params.output_format == "bam" ? "BAM" : "CRAM"
    def t = task.cpus
    def rg = "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA"
    def sort_out = "samtools sort -@ ${t} -m 1G -T tmp_sort --reference ${ref_name} -O ${fmt} -o ${sample}.${params.output_format} -"
    // clipping moves read starts, so sort again afterwards; reads shorter than 30 bp after clipping are dropped
    def clip = primers.name.startsWith("NO_") ? "" :
        "samtools sort -@ ${t} -m 1G -T tmp_pre -u - | samtools ampliconclip --both-ends --filter-len 30 -b ${primers} -f ${sample}.primer_clip.txt -u - | "
    def dedup = (params.mark_duplicates.toString() == 'true' && params.amplicon.toString() != 'true') ?
        "samtools fixmate -m - - | samtools sort -@ ${t} -m 1G -T tmp_sort - | samtools markdup -@ ${t} -s -f ${sample}.markdup.txt --reference ${ref_name} -O ${fmt} - ${sample}.${params.output_format}" :
        "${clip}${sort_out}"
    """
    bwa mem -t ${t} -Y -R '${rg}' ${ref_name} <(cat ${r1s.join(' ')}) <(cat ${r2s.join(' ')}) \
      | ${dedup}
    samtools index -@ ${t} ${sample}.${params.output_format}
    samtools flagstat -@ ${t} ${sample}.${params.output_format} > ${sample}.flagstat.txt
    """
}

// Existing BAM/CRAM given as input: make sure it has an index, and collect alignment statistics.
process INDEX_INPUT_ALN {
    tag "${sample}"
    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
    cpus 4
    input:
    tuple val(sample), path(aln)
    path ref_files
    val ref_name
    output:
    tuple val(sample), path(aln), path("${aln}.{crai,bai}"), emit: aln
    path "${sample}.flagstat.txt", emit: flagstat
    script:
    def ext = aln.name.endsWith(".cram") ? "crai" : "bai"
    """
    samtools index -@ ${task.cpus} ${aln} ${aln}.${ext}
    samtools flagstat -@ ${task.cpus} --reference ${ref_name} ${aln} > ${sample}.flagstat.txt
    """
}

// mosdepth: depth over the target regions, and the share of target bases reaching 10x/20x/30x/50x/100x.
process COVERAGE {
    tag "${sample}"
    container "quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1"
    cpus 4
    publishDir { "${params.out_dir}/${sample}/qc" }, mode: "copy"
    input:
    tuple val(sample), path(aln), path(idx), path(bed)
    path ref_files
    val ref_name
    path named_bed
    output:
    tuple val(sample), path("${sample}.mosdepth.summary.txt"), path("${sample}.thresholds.bed.gz"), path("${sample}.regions.bed.gz"), emit: files
    tuple val(sample), path("${sample}.targets.thresholds.bed.gz"), path("${sample}.targets.regions.bed.gz"), optional: true, emit: targets
    script:
    // second pass on the original, unpadded, named targets: per-gene coverage for the clinical report
    def named = named_bed.name.startsWith("NO_") ? "" :
        "mosdepth -t ${task.cpus} -n --fasta ${ref_name} --by ${named_bed} --thresholds 20,50,100,200,500 ${sample}.targets ${aln}"
    """
    mosdepth -t ${task.cpus} -n --fasta ${ref_name} --by ${bed} --thresholds 10,20,30,50,100 ${sample} ${aln}
    ${named}
    """
}
