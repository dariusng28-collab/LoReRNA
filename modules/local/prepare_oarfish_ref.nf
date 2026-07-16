// ============================================================
// PREPARE_OARFISH_REFERENCE
//
// Runs once after all IsoQuant jobs complete (collected channel).
// Performs only the reference preparation steps — BAM cleaning has
// been moved to a per-sample CLEAN_BAM process (clean_bam.nf) so
// that cleaning can run in parallel with this GTF/FASTA prep.
//
//   Step 1  — gffcompare structural merge across all per-sample
//             transcript_models.gtf files from IsoQuant
//   Step 2  — gffread: merged GTF + genome FASTA → transcript FASTA
//             Validates FASTA count == GTF transcript count
//   Step 3  — tx2gene: extract TCONS_ID → gene_name mapping
//             Falls back to XLOC_ gene_id for novel transcripts
//   Step 4  — Validate FASTA count == tx2gene row count
//             Prevents silent tximport failures from missing mappings
//
// Inputs:
//   gtf_files     — collected channel of all transcript_models.gtf
//   genome_fasta  — GRCh38.primary_assembly.genome.fa
//
// Outputs:
//   merged_gtf    — gffcmp.combined.gtf
//   merged_fa     — merged_expressed_transcripts.fa
//   tx2gene       — tx2gene.tsv (TCONS_ID → gene_name)
//   stats         — gffcmp.stats
//
// Changelog v1.0.0:
//   Removed Step 5 (BAM cleaning) — moved to CLEAN_BAM process.
//   CLEAN_BAM runs per-sample in parallel and no longer blocked by
//   the GTF merge step, reducing wall-clock time for large cohorts.
// ============================================================

process PREPARE_OARFISH_REFERENCE {

    tag "all_samples"
    label 'process_high'   // nextflow.config: 8 CPUs, 32 GB

    publishDir "${params.outdir}/04_oarfish_reference", mode: params.publish_mode

    input:
    path   gtf_files         // collected: all *.transcript_models.gtf
    path   genome_fasta

    output:
    path "merged_expressed_transcripts.fa",  emit: merged_fa
    path "tx2gene.tsv",                      emit: tx2gene
    path "gffcmp.combined.gtf",              emit: merged_gtf
    path "gffcmp.stats",                     emit: stats

    script:
    """
    set -euo pipefail

    TMPWORK="\${PWD}/tmp_work"
    mkdir -p "\${TMPWORK}"

    echo "============================================"
    echo " PREPARE_OARFISH_REFERENCE"
    echo " Genome   : ${genome_fasta}"
    echo "============================================"

    # ── Step 1: gffcompare structural merge ───────────────────────────────
    echo "[1/4] Building GTF list and running gffcompare structural merge..."

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
    echo "[2/4] Building transcriptome FASTA via gffread..."

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
    echo "[3/4] Building tx2gene mapping..."

    # POSIX-compatible awk using match()+substr() instead of 3-arg match().
    # Priority: gene_name > cmp_ref_gene (gffcompare attribute for annotated
    # transcripts) > gene_id (XLOC locus ID for novel transcripts).
    awk '
        \$3 == "transcript" {
            tx   = ""
            gene = ""

            if (match(\$0, /transcript_id "[^"]+"/)) {
                tx = substr(\$0, RSTART+15, RLENGTH-16)
            } else next

            if (match(\$0, /gene_name "[^"]+"/)) {
                gene = substr(\$0, RSTART+11, RLENGTH-12)
            } else if (match(\$0, /cmp_ref_gene "[^"]+"/)) {
                gene = substr(\$0, RSTART+14, RLENGTH-15)
            } else if (match(\$0, /gene_id "[^"]+"/)) {
                gene = substr(\$0, RSTART+9, RLENGTH-10)
            } else next

            if (tx != "" && gene != "") print tx "\t" gene
        }
    ' gffcmp.combined.gtf > tx2gene.tsv

    sort -u tx2gene.tsv -o tx2gene.tsv

    N_TX2GENE=\$(wc -l < tx2gene.tsv)
    echo "tx2gene complete — \${N_TX2GENE} mappings"

    # ── Step 4: Fill gaps and validate FASTA vs tx2gene ──────────────────
    # Some novel transcripts from gffcompare may not parse cleanly with awk.
    # For any FASTA entry missing from tx2gene, use transcript_id as gene_name.
    # These will be quantified by oarfish and appear in DTE; in DGE they
    # aggregate under their own transcript ID (standard practice for novel tx).
    echo "[4/4] Filling tx2gene gaps and validating..."

    grep "^>" merged_expressed_transcripts.fa | sed 's/>//' | awk '{print \$1}' | sort         > "\${TMPWORK}/fasta_ids.txt"
    awk '{print \$1}' tx2gene.tsv | sort         > "\${TMPWORK}/tx2gene_ids.txt"

    N_MISSING=\$(comm -23 "\${TMPWORK}/fasta_ids.txt" "\${TMPWORK}/tx2gene_ids.txt" | wc -l)
    if [ "\${N_MISSING}" -gt 0 ]; then
        echo "WARNING: \${N_MISSING} transcripts have no gene mapping — using transcript_id as gene_name"
        comm -23 "\${TMPWORK}/fasta_ids.txt" "\${TMPWORK}/tx2gene_ids.txt" |             awk '{print \$1"\t"\$1}' >> tx2gene.tsv
    fi

    N_FASTA=\$(grep -c "^>" merged_expressed_transcripts.fa)
    N_TX=\$(wc -l < tx2gene.tsv)

    if [ "\${N_FASTA}" -ne "\${N_TX}" ]; then
        echo "ERROR: FASTA sequences (\${N_FASTA}) != tx2gene entries (\${N_TX}) after gap fill"
        exit 1
    fi
    echo "Validation passed — \${N_FASTA} transcripts in FASTA and tx2gene (\${N_MISSING} used transcript_id as gene_name)"

    rm -rf "\${TMPWORK}"

    echo "============================================"
    echo " PREPARE_OARFISH_REFERENCE complete"
    echo " Merged GTF    : gffcmp.combined.gtf (\${N_MERGED} transcripts)"
    echo " Transcriptome : merged_expressed_transcripts.fa (\${N_FASTA} sequences)"
    echo " tx2gene       : tx2gene.tsv (\${N_TX} mappings)"
    echo "============================================"
    """
}
