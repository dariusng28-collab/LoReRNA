process ISOQUANT {
    label 'process_high'
    publishDir "${params.outdir}/02_isoquant", mode: params.publish_mode

    input:
    path miser_bams
    path reference_fasta
    path reference_gtf

    output:
    path("${params.isoquant_prefix}.transcript_models.gtf"), emit: isoquant_gtf
    path("combined_transcript_counts.tsv"), emit: isoquant_counts
    path("${params.isoquant_prefix}.exon_grouped_counts.tsv"), optional: true, emit: isoquant_exon_counts

    script:
    def bamArgs = miser_bams.collect { it.getName() }.join(' ')
    def extra = params.isoquant_args ?: ''
    """
    ${params.isoquant_exe} \
      --reference ${reference_fasta} \
      --genedb ${reference_gtf} \
      --bam ${bamArgs} \
      --data_type ${params.isoquant_data_type} \
      --threads ${task.cpus} \
      --prefix ${params.isoquant_prefix} \
      --output isoquant_run \
      ${extra}

    cp isoquant_run/${params.isoquant_prefix}.transcript_models.gtf ./
    cp isoquant_run/combined_transcript_counts.tsv ./

    if [[ -f isoquant_run/${params.isoquant_prefix}.exon_grouped_counts.tsv ]]; then
      cp isoquant_run/${params.isoquant_prefix}.exon_grouped_counts.tsv ./
    fi
    """

    stub:
    """
    cat > ${params.isoquant_prefix}.transcript_models.gtf <<'EOF'
    chr1\tIsoQuant\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";
    EOF
    cat > combined_transcript_counts.tsv <<'EOF'
    feature_id\tsampleA\tsampleB
    TX1\t10\t12
    EOF
    cat > ${params.isoquant_prefix}.exon_grouped_counts.tsv <<'EOF'
    feature_id\tsampleA\tsampleB
    exon1\t5\t6
    EOF
    """
}
