process BAM_TO_FASTQ {
    label 'process_medium'
    tag "${meta.sample}"
    publishDir "${params.outdir}/04_fastq", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.fastq.gz"), emit: fastq

    script:
    """
    ${params.samtools_exe} fastq ${bam} | gzip -c > ${meta.sample}.fastq.gz
    """

    stub:
    """
    cat > ${meta.sample}.fastq <<'EOF'
    @read1
    ACGT
    +
    IIII
    EOF
    gzip -c ${meta.sample}.fastq > ${meta.sample}.fastq.gz
    """
}
