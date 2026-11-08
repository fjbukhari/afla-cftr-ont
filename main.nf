#!/usr/bin/env nextflow
// AFLA CFTR workflow: small variants, poly-T/TG, exon copy number and a one-page report per sample, for targeted CFTR
// sequencing on Illumina or Nanopore. Works with any CFTR panel. Becomes main.nf in the shared packages.
// Education and research use only, not for clinical use.

include { BWA_INDEX; PREPARE_REF; PREPARE_BED; PREPARE_BED as PREPARE_REPORT_BED } from './modules/reference'
include { FASTP; ALIGN; COVERAGE } from './modules/reads'
include { CALL_GATK; NORMALISE_FILTER; ANNOTATE_VEP } from './modules/variants'
include { CFTR_POLYT_READS; CFTR_POLYT; CFTR_CNV; CFTR_REPORT; MINIMAP2_INDEX; ALIGN_ONT; SORT_ONT; CALL_BCFTOOLS_ONT; CALL_CLAIR3; WHATSHAP } from './modules/cftr'

// "S1_S8_L001" -> "S1": drop the Illumina sample-sheet number and lane, so lanes of one sample merge
def sampleFromFastq(String id) {
    return id.replaceAll(/(_S\d+)?(_L\d{3})?$/, '')
}

def flag(Object x) {
    if (x instanceof Boolean) return x
    return x != null && x.toString().trim().toLowerCase() in ['true', 'yes', '1', 'on']
}

workflow {
    def ont = params.platform == 'ont'
    if (!(params.platform in ['illumina', 'ont'])) {
        error "platform must be 'illumina' or 'ont'."
    }
    if (!params.fastq || !params.ref) {
        error "Give the FASTQ folder and the GRCh38 reference FASTA."
    }
    if (flag(params.run_annotation) && !params.vep_cache) {
        error "Annotation needs the Ensembl VEP cache folder (vep_cache). Set it, or untick 'Annotate variants'."
    }
    if (params.cftr2_list && !file(params.cftr2_list).exists()) {
        error "CFTR2 list not found: ${params.cftr2_list}"
    }

    // ---------------------------------------------------------------- reference
    def fa = file(params.ref, checkIfExists: true)
    def ref_name = fa.name
    def fai = file("${fa}.fai")
    def dict = file("${fa.parent}/${fa.baseName}.dict")
    def core
    if (fai.exists() && dict.exists()) {
        core = channel.value([fa, fai, dict])
    } else {
        PREPARE_REF(fa)
        core = PREPARE_REF.out.fai.combine(PREPARE_REF.out.dict).map { f, d -> [fa, f, d] }
    }
    def ref_files
    if (ont) {
        ref_files = core.collect()
    } else {
        def idx = ['amb', 'ann', 'bwt', 'pac', 'sa'].collect { e -> file("${fa}.${e}") }
        def bwa_idx
        if (idx.every { f -> f.exists() }) {
            bwa_idx = channel.value(idx)
        } else if (flag(params.build_bwa_index)) {
            bwa_idx = BWA_INDEX(fa).collect()
        } else {
            error "No bwa index found next to ${fa}. Tick 'Build bwa index' (about 1 hour, done once) or choose an indexed reference."
        }
        ref_files = core.combine(bwa_idx).collect()
    }
    def res = "${projectDir}/cftr/resources"
    def exons = file("${res}/CFTR_NM_000492.4_exons.tsv", checkIfExists: true)
    def deletions = file("${res}/CFTR_large_deletions.tsv", checkIfExists: true)
    def report_targets = file("${res}/CFTR_targets_GRCh38.bed", checkIfExists: true)
    def calling_targets = params.bed ? file(params.bed, checkIfExists: true) : report_targets

    // ---------------------------------------------------------------- reads -> alignment
    def d = params.fastq
    def aln
    if (ont) {
        // <folder>/<sample>/*.fastq.gz (demultiplexer or MinKNOW barcode folders) or <folder>/<sample>.fastq.gz
        def root = file(d, checkIfExists: true)
        def reads = channel
            .fromPath(["${d}/**.fastq.gz", "${d}/**.fq.gz", "${d}/**.fastq", "${d}/**.fq"])
            .map { f -> [f.parent == root ? f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '') : f.parent.name, f] }
            .filter { s, _f -> !(s.toLowerCase() in ['unclassified', 'undetermined']) }
            .groupTuple()
            .ifEmpty { error "No .fastq or .fastq.gz files found in ${d}." }
        def mmi = file("${fa.parent}/${fa.baseName}.map-ont.mmi")
        def index = mmi.exists() ? channel.value(mmi) : MINIMAP2_INDEX(fa).collect()
        SORT_ONT(ALIGN_ONT(reads, index))
        aln = SORT_ONT.out
    } else {
        def pairs = channel
            .fromFilePairs(["${d}/**_R{1,2}_001.f*q.gz", "${d}/**_R{1,2}.f*q.gz", "${d}/**_{1,2}.f*q.gz"], size: 2, flat: true)
            .map { id, r1, r2 -> [sampleFromFastq(id), r1, r2] }
            .filter { s, _r1, _r2 -> s.toLowerCase() != 'undetermined' }
            .groupTuple()
            .ifEmpty { error "No paired .fastq.gz files (*_R1/_R2 or *_1/_2) found in ${d}." }
        FASTP(pairs, file("${projectDir}/assets/illumina_adapters.fasta.txt"))
        def primers = params.primer_bed ? file(params.primer_bed, checkIfExists: true) : file("${projectDir}/assets/NO_RESOURCE")
        ALIGN(FASTP.out.reads, ref_files, ref_name, primers)
        aln = ALIGN.out.aln
    }

    // ---------------------------------------------------------------- coverage and small variants
    // Variants are called over the chosen panel; the report always judges coverage against the same CFTR target list.
    def fai_ch = ref_files.map { fs -> fs.find { f -> f.name.endsWith('.fai') } }
    PREPARE_BED(calling_targets, fai_ch)
    PREPARE_REPORT_BED(report_targets, fai_ch)
    def aln_bed = aln.combine(PREPARE_BED.out.padded)
    COVERAGE(aln_bed, ref_files, ref_name, PREPARE_REPORT_BED.out.named)
    def clair3 = ont && params.ont_caller == 'clair3'
    def raw = !ont ? CALL_GATK(aln_bed, ref_files, ref_name) :
        (clair3 ? CALL_CLAIR3(aln_bed, ref_files, ref_name) : CALL_BCFTOOLS_ONT(aln_bed, ref_files, ref_name))
    NORMALISE_FILTER(raw, ref_files, ref_name)
    def called = NORMALISE_FILTER.out.vcf
    def phasing = ont && flag(params.run_phasing)
    if (phasing) {
        called = WHATSHAP(called.join(aln), ref_files, ref_name)
    }

    def annotated = called.map { s, v, _t -> [s, v] }
    if (flag(params.run_annotation)) {
        def extra = [file("${projectDir}/assets/NO_QC")]
        def args = ''
        if (params.clinvar_vcf) {
            def cv = file(params.clinvar_vcf, checkIfExists: true)
            extra = [cv, file("${cv}.tbi", checkIfExists: true)]
            args = "--custom file=${cv.name},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN"
        }
        ANNOTATE_VEP(called, file(params.vep_cache, checkIfExists: true), extra, args)
        annotated = ANNOTATE_VEP.out.vcf
    }

    // ---------------------------------------------------------------- CFTR-specific steps and report
    CFTR_POLYT(CFTR_POLYT_READS(aln), ref_files, ref_name, exons)
    CFTR_CNV(COVERAGE.out.targets.map { _s, _t, r -> r }.collect(), exons, deletions)

    def cftr2 = params.cftr2_list ? file(params.cftr2_list) : file("${projectDir}/assets/NO_RESOURCE")
    def caller = !ont ? 'GATK HaplotypeCaller' :
        (clair3 ? "Clair3 (${params.clair3_model})" : "bcftools (${params.bcftools_profile})") + (phasing ? ', WhatsHap' : '')
    def stamp = new Date().format('yyyy-MM-dd')
    def extra_args = "--caller '${caller}' --version ${workflow.manifest.version} --date ${stamp}" +
        ((params.primer_bed && !ont) ? ' --primers-clipped' : '')
    def per_sample = annotated
        .join(COVERAGE.out.targets)
        .join(CFTR_POLYT.out)
        .combine(CFTR_CNV.out.json)
    CFTR_REPORT(per_sample, ref_files, ref_name, exons, cftr2, extra_args)
}
