process MISER {
    label 'process_high'
    tag "${meta.sample}"
    publishDir "${params.outdir}/01_miser", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path reference_fasta
    path reference_bed12

    output:
    tuple val(meta), path("${meta.sample}.miser.bam"), emit: corrected_bam

    script:
    def extra = params.miser_args ?: ''
    """
    ${params.miser_exe} \
      --input-bam ${bam} \
      --reference ${reference_fasta} \
      --annotation ${reference_bed12} \
      --threads ${task.cpus} \
      --output ${meta.sample}.miser.bam \
      ${extra}
    """

    stub:
    """
    touch ${meta.sample}.miser.bam
    """
}
