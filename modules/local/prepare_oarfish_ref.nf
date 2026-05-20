// ============================================================
// PREPARE_OARFISH_REFERENCE
//
// Replicates oarfish_prep.sh Steps 1–5 exactly:
//   Step 1  — gffcompare structural merge across all per-sample
//             transcript_models.gtf files from IsoQuant
//             Deduplicates by intron chain (not transcript_id)
//             Retains all structurally distinct isoforms including
//             MisER-rescued exon variants
//   Step 2  — gffread: merged GTF + genome FASTA → transcript FASTA
//             Validates FASTA sequence count == GTF transcript count
//   Step 3  — tx2gene: extract TCONS_ID → gene_name mapping
//             Falls back to XLOC_ gene_id for novel transcripts
//             with no gene_name annotation
//   Step 4  — Validate FASTA sequence count == tx2gene row count
//             Prevents silent tximport failures from missing mappings
//   Step 5  — BAM cleaning: strip duplicate ts tags from all MisER
//             sorted BAMs and filter secondary/supplementary records
//             (samtools view -F 2304 -x ts)
//             MisER adds ts tag on top of existing minimap2 ts tag.
//             oarfish's strict BAM parser rejects duplicate tags:
//             InvalidData(DuplicateTag(Tag("ts")))
//             -F 2304 = exclude secondary (256) + supplementary (2048)
//
// Why this module collects all GTFs before running:
//   All samples must be quantified against the same reference.
//   MisER-rescued isoforms may appear in only one condition.
//   Without merging, reads get forced into wrong isoforms and
//   counts become incomparable across conditions.
//   .collect() holds the channel until all GTFs are available.
//
// Changelog v1.1.1 (fix #7):
//   /tmp removed — all temp files now written to ${PWD}/tmp_work
//   within the Nextflow process work directory to avoid ramfs
//   exhaustion and cross-job collisions on shared HPC nodes.
//
//   GTF glob changed from 'ls -1 *.gtf' (which would silently
//   include any annotation GTF staged into the work dir) to
//   'ls -1 *.transcript_models.gtf' — a pattern that matches only
//   the IsoQuant output files staged by the preceding step.
//
// Inputs:
//   gtf_files     — collected channel of all transcript_models.gtf
//   miser_bams    — collected channel of (meta, miser.sorted.bam) tuples
//   genome_fasta  — GRCh38.primary_assembly.genome.fa
//
// Outputs:
//   merged_gtf    — gffcmp.combined.gtf
//   merged_fa     — merged_expressed_transcripts.fa
//   tx2gene       — tx2gene.tsv (TCONS_ID → gene_name)
//   clean_bams    — tuple (meta, sample.clean.bam) per sample
//   stats         — gffcmp.stats (gffcompare merge statistics)
// ============================================================

process PREPARE_OARFISH_REFERENCE {

    tag "all_samples"
    label 'process_high'   // nextflow.config: 8 CPUs, 32 GB | UCL: 4 CPUs, 32 GB

    publishDir "${params.outdir}/04_oarfish_reference", mode: params.publish_mode

    input:
    path   gtf_files         // collected: all *.transcript_models.gtf
    tuple val(metas), path(sorted_bams)  // collected: all (meta, miser.sorted.bam)
    path   genome_fasta

    output:
    path "merged_expressed_transcripts.fa",  emit: merged_fa
    path "tx2gene.tsv",                      emit: tx2gene
    path "gffcmp.combined.gtf",              emit: merged_gtf
    path "gffcmp.stats",                     emit: stats
    tuple val(metas), path("*.clean.bam"),   emit: clean_bams

    script:
    """
    set -euo pipefail

    # Local temp dir scoped to the Nextflow work directory.
    # Avoids /tmp ramfs exhaustion and cross-job collisions on HPC nodes.
    TMPWORK="\${PWD}/tmp_work"
    mkdir -p "\${TMPWORK}"

    echo "============================================"
    echo " PREPARE_OARFISH_REFERENCE"
    echo " Genome   : ${genome_fasta}"
    echo " GTF files: ${gtf_files}"
    echo "============================================"

    # ── Step 1: gffcompare structural merge ───────────────────────────────
    # Only glob *.transcript_models.gtf — prevents accidentally including
    # an annotation GTF if one were ever staged into the work directory.
    echo "[1/5] Building GTF list and running gffcompare structural merge..."

    ls -1 *.transcript_models.gtf > gtf_list.txt
    echo "GTF files to merge:"
    cat gtf_list.txt

    ${params.gffcompare_exe} \\
        -i gtf_list.txt \\
        -o gffcmp

    if [ ! -f "gffcmp.combined.gtf" ]; then
        echo "ERROR: gffcompare did not produce gffcmp.combined.gtf"
        exit 1
    fi

    N_MERGED=\$(awk '\$3=="transcript"' gffcmp.combined.gtf | wc -l)
    echo "gffcompare complete — \${N_MERGED} structurally unique transcripts"
    echo "Stats:"
    cat gffcmp.stats

    # ── Step 2: gffread — build transcriptome FASTA ───────────────────────
    echo "[2/5] Building transcriptome FASTA via gffread..."

    ${params.gffread_exe} gffcmp.combined.gtf \\
        -g "${genome_fasta}" \\
        -w merged_expressed_transcripts.fa

    N_SEQS=\$(grep -c "^>" merged_expressed_transcripts.fa)
    N_TXP=\$(awk '\$3=="transcript"' gffcmp.combined.gtf | wc -l)
    echo "FASTA complete — \${N_SEQS} sequences"

    if [ "\${N_SEQS}" -ne "\${N_TXP}" ]; then
        echo "WARNING: FASTA sequences (\${N_SEQS}) != GTF transcripts (\${N_TXP})"
        echo "Check for contigs in GTF missing from genome FASTA:"
        awk '\$1!="#" {print \$1}' gffcmp.combined.gtf | sort -u > "\${TMPWORK}/gtf_chroms.txt"
        grep "^>" "${genome_fasta}" | sed 's/>//' | awk '{print \$1}' | sort -u > "\${TMPWORK}/fa_chroms.txt"
        comm -23 "\${TMPWORK}/gtf_chroms.txt" "\${TMPWORK}/fa_chroms.txt" || true
    fi

    # ── Step 3: tx2gene mapping ───────────────────────────────────────────
    # gffcompare replaces gene_id with XLOC_ internal IDs.
    # Use gene_name (biological symbol) with XLOC_ fallback for novel txps.
    echo "[3/5] Building tx2gene mapping..."

    awk '
        \$3 == "transcript" {
            transcript_id = ""
            gene_id       = ""
            gene_name     = ""
            n = split(\$9, fields, ";")
            for (i = 1; i <= n; i++) {
                gsub(/^[ \\t]+|[ \\t]+\$/, "", fields[i])
                if (fields[i] ~ /^transcript_id/) {
                    split(fields[i], a, "\\"")
                    transcript_id = a[2]
                }
                if (fields[i] ~ /^gene_id/) {
                    split(fields[i], a, "\\"")
                    gene_id = a[2]
                }
                if (fields[i] ~ /^gene_name/) {
                    split(fields[i], a, "\\"")
                    gene_name = a[2]
                }
            }
            if (transcript_id != "") {
                if (gene_name != "")
                    print transcript_id "\\t" gene_name
                else if (gene_id != "")
                    print transcript_id "\\t" gene_id
            }
        }
    ' gffcmp.combined.gtf > tx2gene.tsv

    N_TX2GENE=\$(wc -l < tx2gene.tsv)
    echo "tx2gene complete — \${N_TX2GENE} mappings"

    # ── Step 4: Validate FASTA vs tx2gene ────────────────────────────────
    echo "[4/5] Validating FASTA vs tx2gene counts..."

    N_FASTA=\$(grep -c "^>" merged_expressed_transcripts.fa)
    N_TX=\$(wc -l < tx2gene.tsv)

    if [ "\${N_FASTA}" -ne "\${N_TX}" ]; then
        echo "ERROR: FASTA sequences (\${N_FASTA}) != tx2gene entries (\${N_TX})"
        echo "Transcripts in FASTA missing from tx2gene:"
        grep "^>" merged_expressed_transcripts.fa | sed 's/>//' | awk '{print \$1}' | sort > "\${TMPWORK}/fasta_ids.txt"
        awk '{print \$1}' tx2gene.tsv | sort > "\${TMPWORK}/tx2gene_ids.txt"
        comm -23 "\${TMPWORK}/fasta_ids.txt" "\${TMPWORK}/tx2gene_ids.txt" || true
        exit 1
    fi
    echo "Validation passed — \${N_FASTA} transcripts matched in FASTA and tx2gene"

    # ── Step 5: Strip duplicate ts tags from all MisER sorted BAMs ───────
    # -F 2304 = exclude secondary (flag 256) + supplementary (flag 2048)
    # -x ts   = remove ts tag (strand indicator added twice by MisER)
    echo "[5/5] Cleaning all MisER sorted BAMs (strip ts tag, filter secondary/supplementary)..."

    for BAM in ${sorted_bams}; do
        SAMPLE=\$(basename "\${BAM}" .miser.sorted.bam)
        echo "  Cleaning \${SAMPLE}..."
        ${params.samtools_exe} view -h -F 2304 -x ts "\${BAM}" \\
            | ${params.samtools_exe} view -b -o "\${SAMPLE}.clean.bam"
        echo "  Done: \${SAMPLE}.clean.bam"
    done

    EXPECTED_CLEAN=${metas.size()}
    ACTUAL_CLEAN=\$(ls *.clean.bam | wc -l)
    if [ "\${ACTUAL_CLEAN}" -ne "\${EXPECTED_CLEAN}" ]; then
        echo "ERROR: Expected \${EXPECTED_CLEAN} clean BAMs but found \${ACTUAL_CLEAN}"
        ls -lh *.clean.bam || true
        exit 1
    fi

    # Clean up local temp dir
    rm -rf "\${TMPWORK}"

    echo "============================================"
    echo " PREPARE_OARFISH_REFERENCE complete"
    echo " Merged GTF    : gffcmp.combined.gtf"
    echo " Transcriptome : merged_expressed_transcripts.fa (\${N_FASTA} sequences)"
    echo " tx2gene       : tx2gene.tsv (\${N_TX} mappings)"
    echo " Clean BAMs    : \$(ls *.clean.bam | wc -l) samples"
    echo "============================================"
    """
}
