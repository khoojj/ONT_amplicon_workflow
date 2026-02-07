# Installation & Setup Guide

This document covers the complete setup for the workflow, specifically optimized for high-performance servers (64+ cores, 256GB+ RAM).

---

## 1. Prerequisites

Ensure your system has the following installed:

* **Java 11+**: Required to run Nextflow.
* **Conda or Mamba**: To manage software dependencies automatically.
* **Nextflow**: (v21.10.3 or higher). Install using Conda.

---
## 2. Install Nextflow

```bash
conda create -n nextflow
conda activate nextflow
conda install bioconda::nextflow
```
Make sure you activate the **nextflow** environment before every run.

```bash
conda activate nextflow
```

---

## 3. Download the Pipeline

Clone the repository and enter the directory:

```bash
git clone https://github.com/khoojj/ONT_amplicon_workflow.git
cd /ONT_amplicon_workflow
```

For usage, refer to README.md



