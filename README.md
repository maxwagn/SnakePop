# SnakePop 🐍

A modular **Snakemake workflow for population genomics** from whole-genome resequencing data.

```text
                   /^\/^\
                   _|__|  O|
          \/     /~     \_/ \
           \____|__________/  \
                  \_______     \
                          `\    \
                            |    |
                           /    /
                          /    /
                        /    /
                      /    /      SnakePop
                     |    |
                     |    |
                      \    \
                       \    \______________   /
                        \                   _/
                         \_______________--'

                           ~~~ SSSSSSSSS ~~~
```

---

## Overview

SnakePop is a modular Snakemake workflow designed for population genomic analyses from whole-genome resequencing data.

The workflow currently supports:

- Reference genome preparation
- Read alignment
- BAM processing
- Variant calling
- Genotype filtering
- Site filtering
- Final VCF generation
- Principal Component Analysis (PCA)

Planned future modules include:

- FST
- π (nucleotide diversity)
- Dxy
- Sliding-window genome scans
- ADMIXTURE
- Population assignment
- Phylogenomics
- Demographic inference

---

## Installation

Clone the repository:

```bash
git clone git@github.com:maxwagn/SnakePop.git
cd SnakePop
```

Create the environment:

```bash
mamba env create -f snakepop_environment.yml
conda activate snakepop
```

---

## Input Data

### Reference Genome

Configured in:

```yaml
ref:
  original_fasta:
  assembly_report:
```

### Sample Metadata

The sample metadata file is specified in:

```yaml
sample_table: config/Metadata.tsv
```

At minimum, the metadata table must contain the following columns:

| Column | Required | Description |
|----------|----------|------------|
| id | Yes | Unique sample identifier |
| read_files | Yes | Semicolon-separated paths to paired-end FASTQ files |

Example:

```text
id      species         morphology      read_files
PM1     P. minutus      morph_1         /path/PM1_R1.fastq.gz;/path/PM1_R2.fastq.gz
PM2     P. minutus      morph_2         /path/PM2_R1.fastq.gz;/path/PM2_R2.fastq.gz
```

The `read_files` column must contain exactly two FASTQ files separated by a semicolon:

```text
/path/sample_R1.fastq.gz;/path/sample_R2.fastq.gz
```

The sample identifier in the `id` column is used throughout the workflow and becomes the sample name in BAM, BCF, VCF, and downstream population genomic analyses.

Additional metadata columns are optional but recommended. These can later be used for:

- PCA colouring (e.g. `morphology`, `species`, `country`)
- Population assignment
- Sample filtering
- Future population genomic analyses

For example:

```text
id  species         morphology  country  location
PM1 P. minutus      morph_1     Norway   Frafjord
PM2 P. minutus      morph_2     Norway   Hoegsfjord
```

---

## Configuration

Main configuration file:

```text
config/config.yaml
```

Important sections:

```yaml
ref:
sample_table:
resources:
callset:
variant_calling:
individual_filter_sets:
site_filter_sets:
popstats:
```

---

## Workflow Overview

```text
FASTQ
  |
  v
Alignment
  |
  v
Processed BAMs
  |
  v
Raw Variant Calling
  |
  v
Genotype Filtering
  |
  v
Site Filtering
  |
  v
Final Callsets
  |
  +--> PCA
  +--> Future population genomic analyses
```

---

## Usage

Display help:

```bash
./snakepop --help
```

### Alignment

```bash
./snakepop alignment --cores 32
```

Output:

```text
results/alignment/
```

---

### Raw Variant Calling

```bash
./snakepop raw_calling --cores 32
```

Output:

```text
results/variants/<callset>/bcf/raw/
```

---

### Filtering

```bash
./snakepop filtering --cores 32
```

Output:

```text
results/variants/<callset>/bcf/final/
reports/filtering_stats/
```

---

### Final Callset Generation

```bash
./snakepop final_callset --cores 32
```

Output:

```text
results/variants/<callset>/vcf/
```

Generated datasets include:

```text
all_sites
variants
biallelic_snps
```

---

### Complete Variant Pipeline

```bash
./snakepop variants --cores 32
```

Runs:

```text
raw_calling
→ filtering
→ final_callset
```

---

### PCA

```bash
./snakepop pca --cores 16
```

Output:

```text
results/popstats/<callset>/pca/
reports/popstats/<callset>/pca/
```

Generated files:

```text
plink_pca.eigenvec
plink_pca.eigenval
pca_scores.tsv
pca_PC1_PC2.png
pca_PC1_PC2.svg
```

Example PCA configuration:

```yaml
popstats:
  pca:
    n_pcs: 10
    maf: 0.05
    geno: 0.2
    mind: 0.2
    color_by: morphology
    label_samples: true
    ld_prune: false
```

---

## Cleanup

Remove intermediate variant-calling files:

```bash
./snakepop clean_intermediates --cores 1
```

Remove Snakemake metadata:

```bash
./snakepop clean_snakemake --cores 1
```

---

## Output Structure

```text
results/
├── alignment/
├── variants/
│   ├── bcf/
│   └── vcf/
└── popstats/

reports/
├── filtering_stats/
└── popstats/
```

---

## Current Status

### Implemented

- Reference preparation
- Read alignment
- BAM processing
- Variant calling
- Variant filtering
- Final VCF generation
- PCA

### Planned

- FST
- π
- Dxy
- Sliding-window scans
- ADMIXTURE
- Population assignment
- Phylogenomics
- Demographic inference

---

## Authors

**Maximilian Wagner** and ****

University of Rijeka, 
Evolutionary Genomics Group  
University of Antwerp

---

## License

