nextflow.enable.dsl = 2

include { MISER } from './modules/local/miser'
include { ISOQUANT } from './modules/local/isoquant'
include { PREPARE_MISER_REFERENCE } from './modules/local/prepare_miser_reference'
include { PREPARE_OARFISH_REFERENCE } from './modules/local/prepare_oarfish_reference'
include { BAM_TO_FASTQ } from './modules/local/bam_to_fastq'
include { OARFISH } from './modules/local/oarfish'
include { MERGE_OARFISH } from './modules/local/merge_oarfish'
include { SWISH } from './modules/local/swish'

def metadataPath = params.metadata ?: params.samplesheet

if( !params.samplesheet ) {
    error "Missing required parameter: --samplesheet"
}
if( !params.reference_fasta ) {
    error "Missing required parameter: --reference_fasta"
}
if( !params.reference_gtf ) {
    error "Missing required parameter: --reference_gtf"
}

workflow {
    main:

    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
    ch_metadata = Channel.fromPath(metadataPath, checkIfExists: true)
    ch_reference_fasta = Channel.fromPath(params.reference_fasta, checkIfExists: true)
    ch_reference_gtf = Channel.fromPath(params.reference_gtf, checkIfExists: true)
    ch_reference_bed12 = params.reference_bed12 ? Channel.fromPath(params.reference_bed12, checkIfExists: true) : null

    ch_samples = ch_samplesheet
        .splitCsv(header: true)
        .map { row ->
            def required = ['sample', 'condition', 'bam']
            required.each { key ->
                if( !row[key] ) {
                    error "Samplesheet is missing required column '${key}' or has an empty value"
                }
            }
            def meta = [sample: row.sample as String, condition: row.condition as String]
            tuple(meta, file(row.bam, checkIfExists: true))
        }

    miser_ref = ch_reference_bed12 ? null : PREPARE_MISER_REFERENCE(ch_reference_gtf)
    miser_annotation = ch_reference_bed12 ?: miser_ref.bed12

    miser_out = MISER(ch_samples, ch_reference_fasta, miser_annotation)

    ch_miser_bams = miser_out.corrected_bam.map { meta, bam -> bam }

    isoquant_out = ISOQUANT(ch_miser_bams.collect(), ch_reference_fasta, ch_reference_gtf)

    oarfish_ref = PREPARE_OARFISH_REFERENCE(isoquant_out.isoquant_gtf, ch_reference_fasta)

    fastq_out = BAM_TO_FASTQ(miser_out.corrected_bam)

    oarfish_out = OARFISH(fastq_out.fastq, oarfish_ref.transcriptome_fasta)

    ch_oarfish_files = oarfish_out.oarfish_quant.flatMap { meta, quant, metaInfo, infreps ->
        [quant, metaInfo, infreps]
    }

    merged_oarfish = MERGE_OARFISH(ch_oarfish_files.collect(), oarfish_ref.tx2gene, ch_samplesheet)

    swish_out = SWISH(
        merged_oarfish.count_matrix,
        merged_oarfish.infreps_manifest,
        oarfish_ref.tx2gene,
        ch_metadata
    )

    emit:
    miser_bams = miser_out.corrected_bam
    isoquant_gtf = isoquant_out.isoquant_gtf
    isoquant_counts = isoquant_out.isoquant_counts
    isoquant_exon_counts = isoquant_out.isoquant_exon_counts
    transcriptome_fasta = oarfish_ref.transcriptome_fasta
    tx2gene = oarfish_ref.tx2gene
    miser_bed12 = miser_annotation
    oarfish_quant_dirs = oarfish_out.oarfish_quant
    oarfish_counts = merged_oarfish.count_matrix
    oarfish_usage = merged_oarfish.isoform_usage
    oarfish_infreps_manifest = merged_oarfish.infreps_manifest
    swish_dte = swish_out.dte_results
    swish_dtu = swish_out.dtu_results
}
