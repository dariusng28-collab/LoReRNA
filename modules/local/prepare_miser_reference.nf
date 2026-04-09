process PREPARE_MISER_REFERENCE {
    label 'process_low'
    publishDir "${params.outdir}/00_reference", mode: params.publish_mode

    input:
    path reference_gtf

    output:
    path("reference.bed12"), emit: bed12

    script:
    """
    ${params.python_exe} ${projectDir}/bin/gtf_to_bed12.py \
      --gtf ${reference_gtf} \
      --output reference.bed12
    """

    stub:
    """
    cat > reference.bed12 <<'EOF'
    chr1	0	100	TX1	0	+	0	100	0	1	100	0
    EOF
    """
}
