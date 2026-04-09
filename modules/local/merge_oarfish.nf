process MERGE_OARFISH {
    label 'process_low'
    publishDir "${params.outdir}/06_oarfish_merged", mode: params.publish_mode

    input:
    path oarfish_outputs
    path tx2gene
    path samplesheet

    output:
    path("oarfish.transcript_counts.tsv"), emit: count_matrix
    path("oarfish.isoform_usage.tsv"), emit: isoform_usage
    path("oarfish.infreps_manifest.tsv"), emit: infreps_manifest

    script:
    def outputArgs = oarfish_outputs.collect { it.getName() }.join(' ')
    """
    ${params.python_exe} ${projectDir}/bin/oarfish_quant_to_matrices.py \
      --samplesheet ${samplesheet} \
      --tx2gene ${tx2gene} \
      --outdir . \
      ${outputArgs}
    """

    stub:
    """
    cat > oarfish.transcript_counts.tsv <<'EOF'
    transcript_id\tgene_id\tgene_name\tsampleA\tsampleB
    TX1\tGENE1\tGENE1\t10\t12
    EOF
    cat > oarfish.isoform_usage.tsv <<'EOF'
    transcript_id\tgene_id\tgene_name\tsampleA\tsampleB
    TX1\tGENE1\tGENE1\t1.0\t1.0
    EOF
    cat > oarfish.infreps_manifest.tsv <<'EOF'
    sample\tquant\tinfreps\tmeta_info
    sampleA\tsampleA.quant\tsampleA.infreps.pq\tsampleA.meta_info.json
    sampleB\tsampleB.quant\tsampleB.infreps.pq\tsampleB.meta_info.json
    EOF
    """
}
