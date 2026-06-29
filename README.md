# SnakePop 🐍

A modular **Snakemake workflow for population genomics** from
whole-genome resequencing data.

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

SnakePop is a modular Snakemake workflow designed for chromosome-level
and scaffold-level whole-genome resequencing projects. It provides an
end-to-end workflow from raw FASTQ files to filtered VCFs and a growing
suite of downstream population genomic and phylogenomic analyses.

## Implemented modules

### Core workflow

-   Reference preparation
-   Read alignment (BWA)
-   BAM processing
-   Variant calling (bcftools)
-   Individual and site filtering
-   Final VCF generation

### Population genomics

-   PCA (PLINK)
-   Genome-wide heterozygosity and inbreeding coefficient (vcftools)
-   Runs of Homozygosity (PLINK)
-   PopGenWindows
    -   FST
    -   dXY
    -   dA
    -   π
-   Manhattan plots
-   Candidate region detection
-   WinPCA
-   HTML summary reports

### Phylogenomics

-   Window-based SNP trees (IQ-TREE)
-   ASTRAL species tree inference
-   Rooted topology-only tree for downstream analyses

### Planned

-   Dsuite
-   ADMIXTURE
-   LD decay
-   IBS/IBD
-   Tajima's D
-   PBS
-   Selection scans
-   Demographic inference

------------------------------------------------------------------------

# Installation

``` bash
git clone git@github.com:maxwagn/SnakePop.git
cd SnakePop

mamba env create -f snakepop_environment.yml
conda activate snakepop
chmod +x snakepop
```

Update an existing environment:

``` bash
mamba env update -n snakepop -f snakepop_environment.yml
```

------------------------------------------------------------------------

# General workflow

    FASTQ
      │
      ▼
    Alignment
      │
      ▼
    Variant calling
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
      ├── Manhattan plots
      ├── WinPCA
      ├── SNP Trees
      └── ASTRAL

------------------------------------------------------------------------

# Wrapper targets

## Core workflow

    alignment
    raw_calling
    filtering
    final_callset
    variants
    all

## Population genomics

    pca
    heterozygosity
    roh
    popgenwindows
    manhattan
    winpca
    popstats

## Phylogenomics

    snptrees_iqtree
    astral
    dsuite (planned)

## Utilities

    clean_intermediates
    clean_snakemake

------------------------------------------------------------------------

# Example commands

``` bash
./snakepop alignment --cores 16

./snakepop variants --cores 32

./snakepop pca --cores 16

./snakepop heterozygosity --cores 4

./snakepop roh --cores 4

./snakepop popgenwindows --cores 16

./snakepop manhattan --cores 8

./snakepop winpca --cores 16

./snakepop snptrees_iqtree --cores 8

./snakepop astral --cores 1
```

------------------------------------------------------------------------

# Population genomics

## PCA

Principal component analysis using PLINK.

Outputs

    results/popstats/<callset>/pca/

------------------------------------------------------------------------

## Genome-wide heterozygosity

Calculated using `vcftools --het`.

Outputs

    heterozygosity.per_sample.tsv
    heterozygosity.per_sample.pdf
    heterozygosity.report.html

------------------------------------------------------------------------

## Runs of Homozygosity

Calculated with PLINK.

Outputs

    plink_roh.hom
    plink_roh.hom.indiv
    roh.per_sample.tsv
    roh.per_sample.pdf
    roh.report.html

------------------------------------------------------------------------

## PopGenWindows

Calculates

-   FST
-   dXY
-   dA
-   π

Outputs genome-wide TSV files.

------------------------------------------------------------------------

## Manhattan module

Produces

-   Manhattan plots
-   Outlier windows
-   Candidate regions
-   WinPCA comparison
-   HTML report

------------------------------------------------------------------------

## WinPCA

Window-based PCA.

Outputs HTML plots and per-window statistics.

------------------------------------------------------------------------

# Phylogenomics

## SNP window trees

Window-based maximum likelihood trees inferred with IQ-TREE.

Temporary files are written to a configurable scratch directory to avoid
generating thousands of permanent files.

Output

    chr*.window_trees.tsv.gz
    all.window_trees.tsv.gz
    snptrees_iqtree.summary.tsv

Each tree stores

-   window coordinates
-   SNP count
-   alignment length
-   Newick tree

------------------------------------------------------------------------

## ASTRAL

ASTRAL summarizes window trees into a species tree.

Outputs

    astral.tree
    astral.topology.tree
    astral_mapping.tsv
    window_trees.newick

`astral.tree`

-   original ASTRAL output

`astral.topology.tree`

-   rooted topology only
-   no branch lengths
-   no support values
-   intended for downstream software such as Dsuite

------------------------------------------------------------------------

# Planned modules

-   Dsuite
-   ADMIXTURE
-   LD decay
-   IBS / IBD
-   PBS
-   Tajima's D
-   Selection scans
-   Demographic inference

------------------------------------------------------------------------

# Authors

Maximilian Wagner

Alex Hooft van Huysduynen

Hannes Svardal

University of Rijeka

Evolutionary Genomics Group

University of Antwerp

