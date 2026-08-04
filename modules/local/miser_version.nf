// ============================================================
// MISER_VERSION
//
// Reports the tool versions inside the MisER container.
//
// Why this is a separate process: LoReRNA ships two images, and a process runs
// in exactly one of them. MISER does not report its own version because adding
// an output to that module would change its task hash and force every MisER
// job — the most expensive step in the pipeline — to re-run. A new process
// adds no hash to any existing one.
//
// Output is the pipe-delimited record format shared with SOFTWARE_VERSIONS:
//
//     <component>|<tool>|<version>
//
// which is that process's single source of truth for both the YAML and the
// MultiQC section. Pipe is safe as a delimiter: no version string or container
// URI contains one.
//
// Version lookups must never fail the pipeline. Each falls back to "unknown"
// rather than exiting non-zero, so a tool that changes its --version output
// degrades the report instead of killing the run.
// ============================================================

process MISER_VERSION {

    tag "miser_container"
    label 'process_single'

    output:
    path "versions_miser.txt", emit: versions

    script:
    """
    # Deliberately no -e: a failed version lookup must not abort the run.
    set -uo pipefail

    # MisER is a Python package: prefer pip metadata, fall back to __version__.
    MISER_V=\$(python3 -m pip show MisER 2>/dev/null | awk '/^Version:/{print \$2}')
    [ -z "\${MISER_V}" ] && MISER_V=\$(python3 -c 'import MisER; print(getattr(MisER,"__version__",""))' 2>/dev/null)
    [ -z "\${MISER_V}" ] && MISER_V="unknown"

    SAMTOOLS_V=\$(samtools --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1)
    [ -z "\${SAMTOOLS_V}" ] && SAMTOOLS_V="unknown"

    PYSAM_V=\$(python3 -c 'import pysam; print(pysam.__version__)' 2>/dev/null)
    [ -z "\${PYSAM_V}" ] && PYSAM_V="unknown"

    PYTHON_V=\$(python3 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    [ -z "\${PYTHON_V}" ] && PYTHON_V="unknown"

    {
        printf 'miser_container|image|%s\\n'    "${params.container_miser}"
        printf 'miser_container|MisER|%s\\n'    "\${MISER_V}"
        printf 'miser_container|samtools|%s\\n' "\${SAMTOOLS_V}"
        printf 'miser_container|pysam|%s\\n'    "\${PYSAM_V}"
        printf 'miser_container|python|%s\\n'   "\${PYTHON_V}"
    } > versions_miser.txt

    echo "MisER container versions:"
    cat versions_miser.txt
    """

    stub:
    """
    {
        printf 'miser_container|image|%s\\n' "${params.container_miser}"
        printf 'miser_container|MisER|stub\\n'
        printf 'miser_container|samtools|stub\\n'
        printf 'miser_container|pysam|stub\\n'
        printf 'miser_container|python|stub\\n'
    } > versions_miser.txt
    """
}
