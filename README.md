# Mg-Clust
## Clustering of ORF sequences in metagenomic data

Mg-Clust is a Nextflow pipeline for computing Operational Protein Units (OPUs) from metagenomic data. It takes paired-end reads from multiple samples, assembles them, predicts ORFs, and clusters them by amino acid sequence identity to produce per-cluster abundance tables.

![MG-Clust workflow](./figures/MG-Clust.png)

---

## Requirements

- [Nextflow](https://www.nextflow.io) >= 23.x
- [Docker](https://www.docker.com)

All tool dependencies (MEGAHIT, BWA, samtools, Picard, FragGeneScanRs, bedtools, bbtools, MMseqs2) are provided via Docker images.

---

## Quick start

```bash
nextflow run main.nf \
    --reads_dir data/reads \
    --output_dir results
```

By default, Nextflow looks for paired-end reads matching `*_R{1,2}*.fastq` inside `reads_dir`.

---

## Input

Paired-end FASTQ files placed in the directory specified by `--reads_dir`. Sample names are inferred automatically from the filenames using the pattern `--reads_pattern`.

Example directory structure:

```
data/reads/
    sample1_R1.fastq
    sample1_R2.fastq
    sample2_R1.fastq
    sample2_R2.fastq
```

---

## Parameters

All parameters can be set in `nextflow.config` or passed on the command line with `--param value`.

| Parameter | Default | Description |
|---|---|---|
| `reads_dir` | `data/reads` | Directory containing input paired-end reads |
| `reads_pattern` | `*_R{1,2}*.fastq` | Glob pattern used to match paired-end read files |
| `output_dir` | `results` | Directory for all pipeline outputs |
| `nslots` | `4` | Number of threads per process |
| `assem_preset` | `meta-sensitive` | MEGAHIT assembly preset |
| `min_contig_len` | `250` | Minimum contig length (bp); shorter contigs are discarded |
| `min_seq` | `5` | Minimum number of assembled contigs required to continue |
| `markdup` | `false` | Run Picard MarkDuplicates to remove PCR duplicates |
| `train_file_name` | `illumina_1` | FragGeneScanRs training model |
| `min_orf_length` | `60` | Minimum ORF length (amino acids); shorter ORFs are discarded |
| `clust_thres` | `0.7` | MMseqs2 sequence identity threshold for clustering |
| `clust_cov_len` | `0.85` | Minimum fraction of aligned residues for clustering (`-c` in MMseqs2) |

---

## Execution profiles

Two profiles are available in `nextflow.config`:

```bash
# Run locally
nextflow run main.nf -profile local

# Run on a SLURM cluster
nextflow run main.nf -profile slurm
```

The SLURM profile submits each process as a job to the `normal` queue. Edit `nextflow.config` to set the correct `--account` and queue name for your cluster.

---

## Pipeline modules

### Module 1 — De novo assembly and read mapping

**Script:** `bin/mg-clust-module-1.py`
**Container:** `ghcr.io/epereira/mg-clust/module-1:latest`
**Tools:** MEGAHIT, BWA-MEM, samtools, Picard

Assembles paired-end reads with MEGAHIT, maps them back to the assembly with BWA-MEM, filters and sorts alignments with samtools, and optionally marks and removes PCR duplicates with Picard MarkDuplicates.

**Inputs:**
- `--reads1` — R1 FASTQ file
- `--reads2` — R2 FASTQ file
- `--sample_name` — sample name used to prefix output files

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--assem_preset` | `meta-sensitive` | MEGAHIT preset (`meta-sensitive`, `meta-large`, etc.) |
| `--min_contig_length` | `250` | Discard contigs shorter than this (bp) |
| `--min_seq` | `5` | Minimum assembled sequences to continue; exits cleanly if below threshold |
| `--markdup` | `false` | Enable Picard MarkDuplicates to remove PCR duplicates |
| `--nslots` | `4` | Threads |
| `--output_dir` | — | Output directory (required) |
| `--overwrite` | `false` | Overwrite output directory if it exists |

**Outputs:**
- `assembly/<sample_name>.contigs.fa` — assembled contigs
- `<sample_name>_sorted.bam` — sorted BAM of reads mapped to contigs

---

### Module 2 — ORF prediction and coverage estimation

**Script:** `bin/mg-clust-module-2.py`
**Container:** `ghcr.io/epereira/mg-clust/module-2:latest`
**Tools:** FragGeneScanRs, bedtools

Predicts ORFs from assembled contigs using FragGeneScanRs, converts coordinates to BED format, and computes per-ORF read counts and mean depth using bedtools.

**Inputs:**
- `--assembly_file` — FASTA file of assembled contigs
- `--bam_file` — sorted BAM of reads mapped to contigs
- `--sample_name` — sample name used to prefix output files

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--train_file_name` | `illumina_1` | FragGeneScanRs training model (see FragGeneScanRs docs for options) |
| `--nslots` | `4` | Threads |
| `--output_dir` | — | Output directory (required) |
| `--overwrite` | `false` | Overwrite output directory if it exists |

**Outputs:**
- `<sample_name>_orfs.faa` — predicted ORF protein sequences
- `<sample_name>_orfs_meancov.tsv` — mean sequencing depth per ORF
- `<sample_name>_orfs_readscov.tsv` — read count per ORF

---

### Module 3 — ORF concatenation, filtering, and MMseqs2 database creation

**Script:** `bin/mg-clust-module-3.py`
**Container:** `ghcr.io/epereira/mg-clust/module-3:latest`
**Tools:** bbduk (bbtools), MMseqs2

Runs once across all samples. Concatenates per-sample ORF FASTA files (prefixing each header with the sample name), filters ORFs by minimum length using bbduk, creates an MMseqs2 sequence database, and merges per-sample coverage tables.

**Inputs:**
- `--orf_files` — list of per-sample `*_orfs.faa` files
- `--meancov_files` — list of per-sample `*_orfs_meancov.tsv` files
- `--readscov_files` — list of per-sample `*_orfs_readscov.tsv` files

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--min_orf_length` | `60` | Minimum ORF length in amino acids; shorter ORFs are discarded |
| `--nslots` | `4` | Threads |
| `--output_dir` | — | Output directory (required) |
| `--overwrite` | `false` | Overwrite output directory if it exists |

**Outputs:**
- `orfs_filt_db*` — MMseqs2 database of filtered ORFs
- `orfs_meancov.tsv` — merged mean coverage table (columns: `sample_name`, `orf_id`, `mean_coverage`)
- `orfs_readscov.tsv` — merged reads coverage table (columns: `sample_name`, `orf_id`, `read_count`)

---

### Module 4 — ORF clustering and abundance table generation

**Script:** `bin/mg-clust-module-4.py`
**Container:** `ghcr.io/epereira/mg-clust/module-4:latest`
**Tools:** MMseqs2, pandas

Runs once. Clusters the filtered ORFs using MMseqs2 `cluster` at the specified sequence identity threshold, maps per-ORF coverage to cluster IDs, and collapses the per-ORF abundance tables to per-cluster abundance summed across samples.

**Inputs:**
- `--orfs_db` — MMseqs2 database prefix (without file extensions)
- `--meancov_table` — merged mean coverage table from Module 3
- `--readscov_table` — merged reads coverage table from Module 3

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--clust_thres` | `0.7` | Sequence identity threshold for MMseqs2 clustering |
| `--clust_cov_len` | `0.85` | Minimum fraction of aligned residues (`-c` in MMseqs2, coverage mode 0) |
| `--nslots` | `4` | Threads |
| `--output_dir` | — | Output directory (required) |
| `--overwrite` | `false` | Overwrite output directory if it exists |

**Outputs (inside a subdirectory named by threshold, e.g. `clust_orfs_id70perc/`):**
- `orfs_clust_id<N>perc.tsv` — cluster membership table (columns: `cluster_id`, `orf_id`)
- `orfs_clust_id<N>perc2meancov.tsv` — per-ORF mean coverage mapped to clusters
- `orfs_clust_id<N>perc2readscov.tsv` — per-ORF read counts mapped to clusters

**Outputs (in root output directory):**
- `orfs_clust_id<N>perc_meancov_workable.tsv` — collapsed mean coverage per cluster per sample
- `orfs_clust_id<N>perc_readscov_workable.tsv` — collapsed read counts per cluster per sample

---

## Output directory structure

```
results/
    module-1/
        <sample_name>/
            assembly/<sample_name>.contigs.fa
            <sample_name>_sorted.bam
    module-2/
        <sample_name>/
            <sample_name>_orfs.faa
            <sample_name>_orfs_meancov.tsv
            <sample_name>_orfs_readscov.tsv
    module-3/
        orfs_filt_db*
        orfs_meancov.tsv
        orfs_readscov.tsv
    module-4/
        clust_orfs_id70perc/
            ...
        orfs_clust_id70perc_meancov_workable.tsv
        orfs_clust_id70perc_readscov_workable.tsv
```

---

## Building Docker images

Images are built from the Dockerfiles in `docker/` and must be run from the **repository root**:

```bash
docker build --network=host -f docker/Dockerfile.module-1 -t ghcr.io/epereira/mg-clust/module-1:latest .
docker build --network=host -f docker/Dockerfile.module-2 -t ghcr.io/epereira/mg-clust/module-2:latest .
docker build --network=host -f docker/Dockerfile.module-3 -t ghcr.io/epereira/mg-clust/module-3:latest .
docker build --network=host -f docker/Dockerfile.module-4 -t ghcr.io/epereira/mg-clust/module-4:latest .
```

To push to GitHub Container Registry:

```bash
echo $CR_PAT | docker login ghcr.io -u YOUR_GITHUB_USERNAME --password-stdin

docker push ghcr.io/epereira/mg-clust/module-1:latest
docker push ghcr.io/epereira/mg-clust/module-2:latest
docker push ghcr.io/epereira/mg-clust/module-3:latest
docker push ghcr.io/epereira/mg-clust/module-4:latest
```
