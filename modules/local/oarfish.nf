process OARFISH {
    label 'process_high'
    tag "${meta.sample}"
    publishDir "${params.outdir}/05_oarfish", mode: params.publish_mode

    input:
    tuple val(meta), path(reads)
    path transcriptome_fasta

    output:
    tuple val(meta),
        path("${meta.sample}.quant"),
        path("${meta.sample}.meta_info.json"),
        path("${meta.sample}.infreps.pq"),
        emit: oarfish_quant

    script:
    def extra = params.oarfish_args ?: ''
    def bootstraps = params.oarfish_num_bootstraps as Integer
    def bootstrapArg = bootstraps > 0 ? "--num-bootstraps ${bootstraps}" : ''
    """
    ${params.oarfish_exe} \
      --reads ${reads} \
      --reference ${transcriptome_fasta} \
      --seq-tech ${params.seq_tech} \
      --threads ${task.cpus} \
      --filter-group ${params.oarfish_filter_group} \
      ${bootstrapArg} \
      --output ${meta.sample} \
      ${extra}
    """

    stub:
    """
    cat > ${meta.sample}.quant <<'EOF'
    transcript_id\tlength\tnum_reads
    TX1\t1000\t10
    EOF
    cat > ${meta.sample}.meta_info.json <<'EOF'
    {"sample":"${meta.sample}","bootstraps":${params.oarfish_num_bootstraps}}
    EOF
    python - <<'PY'
    import pathlib
    pathlib.Path("${meta.sample}.infreps.pq").touch()
    PY
    """
}
