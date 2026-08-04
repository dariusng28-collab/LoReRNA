// ============================================================
// SOFTWARE_VERSIONS
//
// Collects tool versions from the lorerna container, merges in those reported
// by MISER_VERSION for the other image, and emits:
//
//   software_versions.yml      — plain YAML, published to pipeline_info/
//   software_versions_mqc.yml  — MultiQC custom content ("Software Versions")
//
// Structure: every version is recorded once, in a pipe-delimited record
//
//     <component>|<tool>|<version>
//
// and both output formats are generated from that one list. Listing the tools
// separately per format would let the two drift, and would silently drop the
// MisER container's entries from whichever format forgot them.
//
// Design note — why one dump per image rather than nf-core's per-process
// versions.yml: the nf-core convention emits a versions.yml from every process
// because nf-core runs one container per process, so each genuinely has
// different versions to report. LoReRNA runs every process from one of two
// images, so per-process files would repeat the same versions while
// invalidating the task hash of all 13 analysis modules. One dump per image
// carries the same information and leaves every existing task hash untouched.
//
// Version lookups must never fail the pipeline: each falls back to "unknown".
// ============================================================

process SOFTWARE_VERSIONS {

    tag "all_samples"
    label 'process_single'

    input:
    path miser_versions   // versions_miser.txt from MISER_VERSION

    output:
    path "software_versions.yml",     emit: versions
    path "software_versions_mqc.yml", emit: mqc_yml

    script:
    """
    # Deliberately no -e: a failed version lookup must not abort the run.
    set -uo pipefail

    # Print the first version-like token from a command's output, or "unknown".
    # These tools disagree on format ("samtools 1.23.1", "gffcompare v0.12.6",
    # "multiqc, version 1.25.1"), so match the number rather than the layout.
    ver() {
        out=\$("\$@" 2>&1 | head -3)
        v=\$(printf '%s' "\${out}" | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1)
        [ -z "\${v}" ] && v="unknown"
        printf '%s' "\${v}"
    }

    # R package versions. fishpond performs the differential testing, so its
    # version is the one that matters most for reproducing the statistics.
    ${params.rscript_exe} -e 'for (p in c("fishpond","SummarizedExperiment","ggplot2")) {
            v <- tryCatch(as.character(packageVersion(p)), error = function(e) "unknown")
            cat(p, v, sep = "|"); cat("\\n")
        }' > r_pkgs.txt 2>/dev/null || : > r_pkgs.txt

    rpkg() {
        v=\$(awk -F'[|]' -v k="\$1" '\$1==k{print \$2}' r_pkgs.txt | head -1)
        [ -z "\${v}" ] && v="unknown"
        printf '%s' "\${v}"
    }

    # ── Single source of truth ───────────────────────────────────────────────
    {
        printf 'pipeline|LoReRNA|%s\\n'  "${workflow.manifest.version}"
        printf 'pipeline|Nextflow|%s\\n' "${workflow.nextflow.version}"

        printf 'lorerna_container|image|%s\\n'                "${params.container_lorerna}"
        printf 'lorerna_container|samtools|%s\\n'             "\$(ver ${params.samtools_exe} --version)"
        printf 'lorerna_container|IsoQuant|%s\\n'             "\$(ver ${params.isoquant_exe} --version)"
        printf 'lorerna_container|oarfish|%s\\n'              "\$(ver ${params.oarfish_exe} --version)"
        printf 'lorerna_container|gffcompare|%s\\n'           "\$(ver ${params.gffcompare_exe} --version)"
        printf 'lorerna_container|gffread|%s\\n'              "\$(ver ${params.gffread_exe} --version)"
        printf 'lorerna_container|MultiQC|%s\\n'              "\$(ver multiqc --version)"
        printf 'lorerna_container|R|%s\\n'                    "\$(ver ${params.rscript_exe} --version)"
        printf 'lorerna_container|fishpond|%s\\n'             "\$(rpkg fishpond)"
        printf 'lorerna_container|SummarizedExperiment|%s\\n' "\$(rpkg SummarizedExperiment)"
        printf 'lorerna_container|ggplot2|%s\\n'              "\$(rpkg ggplot2)"
    } > versions_all.txt

    cat "${miser_versions}" >> versions_all.txt

    # ── Plain YAML ───────────────────────────────────────────────────────────
    # Records are already grouped by component, so a header is emitted whenever
    # the component changes.
    awk -F'[|]' '
        \$1 != last { printf "%s:\\n", \$1; last = \$1 }
                    { printf "    %s: \\"%s\\"\\n", \$2, \$3 }
    ' versions_all.txt > software_versions.yml

    # ── MultiQC custom content ───────────────────────────────────────────────
    # A table rather than a definition list: samtools appears in both images,
    # so the component column is what keeps the two entries distinguishable.
    # plot_type html because MultiQC applies no styling of its own to custom
    # content, so this markup is what actually renders.
    {
        echo "id: 'lorerna_software_versions'"
        echo "section_name: 'Software Versions'"
        echo "section_href: '${workflow.manifest.homePage}'"
        echo "plot_type: 'html'"
        echo "description: 'Reported by the containers at run time. LoReRNA runs every process in one of two images.'"
        echo "data: |"
        echo "    <table class=\\"table table-condensed\\">"
        echo "    <thead><tr><th>Component</th><th>Tool</th><th>Version</th></tr></thead>"
        echo "    <tbody>"
        awk -F'[|]' '{ printf "    <tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n", \$1, \$2, \$3 }' versions_all.txt
        echo "    </tbody>"
        echo "    </table>"
    } > software_versions_mqc.yml

    echo "Collected software versions:"
    cat software_versions.yml
    """

    stub:
    """
    {
        printf 'pipeline|LoReRNA|%s\\n' "${workflow.manifest.version}"
        printf 'lorerna_container|image|%s\\n' "${params.container_lorerna}"
        printf 'lorerna_container|samtools|stub\\n'
    } > versions_all.txt

    cat "${miser_versions}" >> versions_all.txt

    awk -F'[|]' '
        \$1 != last { printf "%s:\\n", \$1; last = \$1 }
                    { printf "    %s: \\"%s\\"\\n", \$2, \$3 }
    ' versions_all.txt > software_versions.yml

    {
        echo "id: 'lorerna_software_versions'"
        echo "section_name: 'Software Versions'"
        echo "plot_type: 'html'"
        echo "data: |"
        echo "    <table><tbody>"
        awk -F'[|]' '{ printf "    <tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n", \$1, \$2, \$3 }' versions_all.txt
        echo "    </tbody></table>"
    } > software_versions_mqc.yml
    """
}
