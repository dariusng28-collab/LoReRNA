# Pull request

## Description

<!-- What does this change, and why? Link any related issue. -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds capability)
- [ ] Breaking change (existing behaviour or outputs change)
- [ ] Documentation / configuration only

## Checklist

- [ ] CI passes (lint + container CLI contract)
- [ ] `nextflow run main.nf -profile test -stub-run` completes, if the DAG changed
- [ ] `nextflow_schema.json` updated, if parameters were added, removed or renamed
- [ ] `docs/usage.md` / `docs/output.md` updated, if behaviour or outputs changed
- [ ] `CHANGES.md` updated
- [ ] Container `environment*.yml` pins updated, if a tool version changed

## Impact on existing runs

<!-- Does this invalidate the resume cache? Which processes re-run?
     Does it change any published output paths or file formats? -->

## Testing

<!-- How was this verified? Note the profile, executor and dataset used.
     If a tool version changed, say whether the pipeline was re-validated
     end-to-end on real data. -->
