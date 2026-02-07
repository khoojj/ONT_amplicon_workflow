# Nextflow Workflow for ONT Amplicon Analysis

## Overview

This repository contains a **Nextflow workflow** for the analysis of **Oxford Nanopore Technologies (ONT)** amplicon sequencing data for **BLASTn-based identification and taxonomic classification**.

The workflow was developed and optimised mainly for identifying closely related intracellular bacteria, in particular **_Wolbachia_** and **_Rickettsia_**, from ONT amplicon data. It has also been used to analyse **insect COI** amplicon sequences. While the pipeline can be applied to other amplicon targets, its current design choices make it most suitable for **targeted taxonomic identification** rather than broad microbial community profiling.

## Upstream pipeline and credits

This workflow is based on the original **NanoCLUST** pipeline developed by the genomicsITER group:

- NanoCLUST GitHub repository: https://github.com/genomicsITER/NanoCLUST

The original NanoCLUST workflow provides the core methodology for read clustering, consensus generation, and polishing of ONT amplicon data.

For information on third-party tools, libraries, and original authorship used in this workflow and its upstream dependencies, please refer to the **`CREDITS.md`** file included in this repository.

## Key Features and Extensions

Compared to the original NanoCLUST workflow, this pipeline includes the following additions:

1. Primer trimming and read length filtering (**cutadapt**)
2. De novo chimera detection and removal (**vsearch**)
3. BLAST-based taxonomic assignment using a **best-hit strategy**
   - Only the **top BLAST hit** per consensus sequence is retained for downstream summaries
4. Excel-based result outputs (**.xlsx**)
5. OTU table generation across all samples

## Scope and Limitations

- Designed for targeted identification (e.g. closely related taxa such as **_Wolbachia_**, **_Rickettsia_**, and insect COI).
- Not intended for complex community profiling where LCA or multi-hit classification is required.

## Requirements

- Nextflow version used: **22.10.6**
- Conda
- Local or HPC execution supported
- Local NCBI BLAST database and BLAST taxonomy (**`taxdb`**) files

### Installation
Refer to INSTALL.md

### Basic run (Conda profile)
Reads file are named barcode**.fq.gz (ie. barcode06.fq.gz)

```bash
nextflow run main.nf -profile conda \
  --reads "*.fq.gz" \
  --db /path/to/blast_db/nt \
  --tax /path/to/blastdb_taxonomy_dir \
  --outdir results
```

### Change primer sequence (cutadapt)

```bash
nextflow run main.nf -profile conda \
  --reads "*.fq.gz" \
  --cutadapt_primer "YOUR_PRIMER_SEQUENCE" \
  --db /path/to/blast_db/nt \
  --tax /path/to/blastdb_taxonomy_dir
```

Example primer (replace as needed):

```text
AGAGTTTGATCMTGGCTCAG...AAGTCGTAACAAGGTAACCG
```

### Adjust clustering / polishing parameters (examples)

```bash
nextflow run main.nf -profile conda \
  --reads "*.fq.gz" \
  --db /path/to/blast_db/nt \
  --tax /path/to/blastdb_taxonomy_dir \
  --cluster_sel_epsilon 0.15 \
  --min_cluster_size 100 \
  --polishing_reads 300
```

### Running the Pipeline at Scale
For large datasets (e.g., >50 samples), the Nextflow manager requires additional memory to track parallel tasks. It is highly recommended to set the Java Virtual Machine (JVM) overhead before launching the pipeline:

```bash
# Allow the Nextflow controller to use up to 8GB of RAM
export NXF_OPTS="-Xms2g -Xmx8g"

# Run the pipeline
nextflow run main.nf -profile conda \
  --reads "*.fq.gz" \
  --db /path/to/blast_db/nt \
  --tax /path/to/blastdb_taxonomy_dir \
  --outdir results
```

## Author and Citation

Author: **Jing Jing Khoo**

If you use this workflow in your work, please cite this repository.  
If/when a publication becomes available, include the DOI here:

- Publication DOI: **TBD**
