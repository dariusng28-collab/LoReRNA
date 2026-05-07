# Configuring LoReRNA for your HPC cluster

LoReRNA uses a site config file (`-c`) for cluster-specific settings, keeping the pipeline itself portable. This is the same convention used by nf-core pipelines.

## What goes in a site config

A site config sets three things only:

1. **Executor** — SGE, SLURM, LSF, PBS
2. **Resource caps** — `max_cpus`, `max_memory`, `max_time` for your nodes
3. **Container runtime** — bind mounts and cache directory for Singularity/Apptainer

It does **not** change pipeline logic, parameters, or process definitions.

## Getting started

Copy the template and fill in your cluster details:

```bash
cp conf/cluster_template.config conf/your_cluster.config
```

Then edit:

```groovy
// 1. Set your executor
executor {
    name      = 'slurm'   // sge / slurm / lsf / pbs
    queueSize = 50
}

// 2. Set container bind mounts and cache
singularity.runOptions = "-B /your/data/filesystem,${HOME}"
singularity.cacheDir   = "/shared/path/singularity_cache"

// 3. Set resource caps to match your node limits
params {
    max_cpus   = 32
    max_memory = '128.GB'
    max_time   = '72.h'
}

// 4. Set your scheduler's memory syntax in clusterOptions
process {
    executor = 'slurm'
    queue    = 'compute'
    clusterOptions = { "--mem=${task.memory.mega}M" }
}
```

## Memory syntax by scheduler

| Scheduler | clusterOptions syntax |
|-----------|----------------------|
| SGE (standard) | `"-l h_vmem=${mem_per_core}M"` |
| SGE (UCL-style per-core) | `"-l tmem=${mem_per_core}M,h_vmem=${mem_per_core}M"` |
| SLURM | `"--mem=${task.memory.mega}M"` |
| LSF | `"-R 'rusage[mem=${task.memory.mega}MB]'"` |
| PBS | `"-l mem=${task.memory.mega}mb"` |

For SGE clusters that charge memory per core (common at UK universities), divide total memory by CPUs:
```groovy
clusterOptions = {
    def mem_per_core = (task.memory.mega / task.cpus).toInteger()
    "-S /bin/bash -l tmem=${mem_per_core}M,h_vmem=${mem_per_core}M"
}
```

## Running with your site config

```bash
nextflow run dariusng28-collab/LoReRNA \
    -profile singularity \
    -c conf/your_cluster.config \
    --samplesheet samplesheet.csv \
    --reference_fasta /path/to/genome.fa \
    --reference_gtf   /path/to/annotation.gtf \
    --reference_bed12 /path/to/annotation.bed12 \
    --outdir results
```

## Tips

- Run `nextflow run dariusng28-collab/LoReRNA -profile test,singularity -c conf/your_cluster.config` first to verify your config before submitting real data
- Set `queueSize` in the executor block to avoid overwhelming your scheduler with too many concurrent job submissions
- Point `singularity.cacheDir` at a shared filesystem so images are pulled once and shared across all users
- Use `--max_memory` and `--max_time` overrides on the command line to test with reduced resources before a full run
