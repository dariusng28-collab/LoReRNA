# Contributing to LoReRNA

Contributions are welcome — bug reports, documentation improvements,
new site configs, and module improvements are all appreciated.

## Adding a new site config

1. Copy the closest scheduler template:
   ```bash
   cp conf/slurm.config conf/mycluster.config
   ```
2. Edit only the scheduler-specific block (queue, account, scratch resource).
   Do **not** modify `modules/local/` or `main.nf` for site-specific reasons.
3. Test with `conf/validate.config` stacked on top before opening a PR.
4. Open a PR with a brief description of the cluster and institution.

## Reporting a bug

Please include:
- The full Nextflow command used
- The `.nextflow.log` file
- The `work/` subdirectory path for the failed process
- The scheduler and HPC system (optional but helpful)

## Code style

- One process per `.nf` file in `modules/local/`
- All shell steps inside process scripts use `set -euo pipefail`
- Resource labels (`process_low`, `process_medium`, etc.) from `nextflow.config` —
  never hardcode CPUs or memory inside a module
- New Python scripts: `argparse` CLI, docstring at the top, typed where practical
- New R scripts: `commandArgs(trailingOnly = TRUE)` CLI, `suppressPackageStartupMessages`
