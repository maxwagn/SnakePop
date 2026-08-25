# SnakePop 🐍

> **Status:** 🚧 Work in progress... Use at your own risk!!

> A modular **Snakemake workflow for population genomics and phylogenomics** from whole-genome resequencing data.

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
                      /    /      SnakePop v1
                     |    |
                     |    |
                      \    \
                       \    \______________   /
                        \                   _/
                         \_______________--'
                           ~~~ SSSSSSSSS ~~~
```

## Overview

SnakePop is a modular, reproducible **Snakemake** workflow for analysing chromosome- and scaffold-level whole-genome resequencing data. It takes users from raw FASTQ files to high-quality variant datasets and provides an integrated collection of downstream population genomic and phylogenomic analyses.

### Features

- End-to-end workflow from FASTQ to filtered VCFs
- Modular execution of individual analyses
- Configuration through a single `config/config.yaml`
- Scalable from local workstations to HPC clusters
- Reproducible software environments via Conda/Mamba
- Easily extensible with additional downstream analyses

---

# Workflow

```text
FASTQ
  │
  ▼
Reference preparation
  │
  ▼
Read alignment (BWA)
  │
  ▼
Variant calling (bcftools)
  │
  ▼
Filtering
  │
  ▼
Final VCF
  │
  ├── PCA
  ├── Heterozygosity
  ├── ROH
  ├── PopGenWindows
  ├── Watterson's theta (θW)
  ├── Manhattan plots
  ├── WinPCA
  ├── Site-frequency spectra
  ├── SNP trees (IQ-TREE)
  ├── ASTRAL
  ├── Dsuite
  └── Twisst
```

---

# Implemented modules

## Core workflow

- Reference preparation
- Read alignment (BWA)
- BAM processing
- Variant calling (bcftools)
- Genotype and site filtering
- Final VCF generation

## Population genomics

- Principal component analysis (PLINK)
- Genome-wide heterozygosity and inbreeding coefficients
- Runs of homozygosity (ROH)
- PopGenWindows
  - FST
  - dXY
  - dA
  - π
- Watterson's theta (θW)
  - Window-based estimation
  - Population-specific estimates
  - Missing-data-aware chromosome counts
  - Callable-site normalization
- Manhattan plots
- Candidate-region detection
- Windowed PCA
- One- and two-dimensional site-frequency spectra
- HTML summary reports

## Phylogenomics

- Window-based SNP trees (IQ-TREE)
- Species-tree inference with ASTRAL
- Dsuite introgression analyses
- Twisst topology weighting
- Rooted topology-only trees for downstream analyses

## Planned

- ADMIXTURE
- LD decay
- IBS / IBD
- Tajima's D
- PBS
- Selection scans
- Demographic inference

---

# Installation

```bash
git clone git@github.com:maxwagn/SnakePop.git

cd SnakePop

mamba env create -f snakepop_environment.yml

conda activate snakepop

chmod +x snakepop
```

Update an existing environment:

```bash
mamba env update -n snakepop -f snakepop_environment.yml
```

---

# Configuration

SnakePop is configured through `config/config.yaml`.

The configuration file specifies:

- reference genome
- sample metadata
- sequencing libraries
- population assignments
- filtering thresholds
- analysis-specific parameters
- output locations

---

# Wrapper targets

## Core workflow

- `alignment`
- `raw_calling`
- `filtering`
- `final_callset`
- `variants`
- `all`

## Population genomics

- `pca`
- `heterozygosity`
- `roh`
- `popgenwindows`
- `watterson_theta`
- `manhattan`
- `winpca`
- `sfs`
- `popstats`

## Phylogenomics

- `snptrees_iqtree`
- `astral`
- `dsuite`
- `twisst`

## Utilities

- `clean_intermediates`
- `clean_snakemake`

---

# Example usage

```bash
./snakepop alignment --cores 16

./snakepop variants --cores 32

./snakepop pca --cores 16

./snakepop heterozygosity --cores 4

./snakepop roh --cores 4

./snakepop popgenwindows --cores 16

./snakepop watterson_theta --cores 8

./snakepop manhattan --cores 8

./snakepop winpca --cores 16

./snakepop sfs --cores 4

./snakepop snptrees_iqtree --cores 8

./snakepop astral --cores 1

./snakepop dsuite --cores 4

./snakepop twisst --cores 4

./snakepop popstats --cores 32
```

---

# Watterson's theta

SnakePop calculates windowed **Watterson's theta (θW)** directly from the filtered all-sites VCFs.

Watterson's theta is an estimator of the population mutation parameter based on the number of segregating sites:

\[
\theta_W = \frac{S}{a_n}
\]

where:

- \(S\) is the number of segregating sites
- \(n\) is the number of sampled chromosomes
- \(a_n\) is the harmonic sum

\[
a_n = \sum_{i=1}^{n-1}\frac{1}{i}
\]

For genomic windows, SnakePop reports θW per callable nucleotide site.

The implementation accounts for variable missing data by using the number of called chromosomes at each site. Sites below the configured minimum called fraction are excluded.

Example configuration:

```yaml
popstats:

  watterson_theta:
    enabled: true

    population_column: morphology

    populations:
      include: all

    window_size: 60000
    window_step: 30000

    ploidy: 2

    min_called_fraction: 0.8
    min_callable_sites: 1000

    snps_only: true

    script: bin/watterson_theta_snakepop.py
```

Run the analysis with:

```bash
./snakepop watterson_theta --cores 8
```

The final merged output is written to:

```text
results/popstats/<CALLSET>_<REFERENCE>/watterson_theta/watterson_theta.tsv
```

Chromosome-specific results are written to:

```text
results/popstats/<CALLSET>_<REFERENCE>/watterson_theta/per_chrom/
```

The output includes:

- chromosome
- window start and end
- population
- number of samples
- maximum number of chromosomes
- minimum required called chromosomes
- number of callable sites
- number of segregating sites
- mean number of called chromosomes
- harmonic denominator
- Watterson's theta (θW)

---

# Authors

**Maximilian Wagner**  
**Alex Hooft van Huysduynen**  
**Hannes Svardal**

Evolutionary Genomics Group  
University of Antwerp
