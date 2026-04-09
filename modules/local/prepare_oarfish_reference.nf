process PREPARE_OARFISH_REFERENCE {
    label 'process_medium'
    publishDir "${params.outdir}/03_oarfish_reference", mode: params.publish_mode

    input:
    path isoquant_gtf
    path reference_fasta

    output:
    path("isoquant.transcriptome.fa"), emit: transcriptome_fasta
    path("tx2gene.tsv"), emit: tx2gene

    script:
    """
    ${params.gffread_exe} ${isoquant_gtf} -g ${reference_fasta} -w isoquant.transcriptome.fa
    ${params.python_exe} ${projectDir}/bin/gtf_to_tx2gene.py --gtf ${isoquant_gtf} --output tx2gene.tsv
    """

    stub:
    """
    cat > isoquant.transcriptome.fa <<'EOF'
    >TX1
    ACGT
    EOF
    cat > tx2gene.tsv <<'EOF'
    transcript_id\tgene_id\tgene_name
    TX1\tGENE1\tGENE1
    EOF
    """
}
