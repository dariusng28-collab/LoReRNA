process SWISH {
    label 'process_medium'
    publishDir "${params.outdir}/07_swish", mode: params.publish_mode

    input:
    path count_matrix
    path infreps_manifest
    path tx2gene
    path metadata

    output:
    path("swish.dte.tsv"), emit: dte_results
    path("swish.dtu.tsv"), emit: dtu_results

    script:
    def pairArg = params.swish_pair_col ? "--pair-col ${params.swish_pair_col}" : ''
    def batchArg = params.swish_batch_col ? "--batch-col ${params.swish_batch_col}" : ''
    """
    ${params.rscript_exe} ${projectDir}/bin/run_swish.R \
      --counts ${count_matrix} \
      --manifest ${infreps_manifest} \
      --tx2gene ${tx2gene} \
      --metadata ${metadata} \
      --condition-col ${params.swish_condition_col} \
      ${pairArg} \
      ${batchArg} \
      --min-count ${params.swish_min_count} \
      --min-n ${params.swish_min_n} \
      --nperms ${params.swish_nperms} \
      --alpha ${params.swish_alpha} \
      --outdir .
    """

    stub:
    """
    cat > swish.dte.tsv <<'EOF'
    transcript_id\tgene_id\tlog2FC\tqvalue
    TX1\tGENE1\t1.0\t0.01
    EOF
    cat > swish.dtu.tsv <<'EOF'
    transcript_id\tgene_id\tlog2FC\tqvalue
    TX1\tGENE1\t0.5\t0.03
    EOF
    """
}
