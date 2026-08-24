#!/usr/bin/env python3

import argparse
import math
import os

import numpy as np
import pandas as pd
import pysam


def clean_pop_name(value):
    """
    Make population names compatible with SnakePop paths/wildcards.
    """
    return (
        str(value)
        .replace(" ", "_")
        .replace("/", "_")
        .replace(";", "_")
    )


def harmonic_a1(n):
    """
    a_n = sum(1/i), i = 1 ... n-1
    """
    if n < 2:
        return 0.0

    return float(
        np.sum(
            1.0 / np.arange(1, n, dtype=float)
        )
    )


def is_nucleotide_site(record):
    """
    Return True for invariant nucleotide positions and SNP positions.

    Indels, symbolic alleles, etc. are excluded when --snps-only is used.
    """

    if len(record.ref) != 1:
        return False

    # Invariant site.
    if not record.alts:
        return True

    for alt in record.alts:
        if alt is None:
            return False

        if len(alt) != 1:
            return False

        if alt.startswith("<"):
            return False

    return True


def called_alleles(record, sample_names):
    """
    Extract all called chromosome-level alleles for a population.
    """

    alleles = []

    for sample in sample_names:

        gt = record.samples[sample].get("GT")

        if gt is None:
            continue

        for allele in gt:
            if allele is not None and allele >= 0:
                alleles.append(allele)

    return alleles


def window_indices(
    pos,
    n_windows,
    window_size,
    window_step,
):
    """
    Return indices of all overlapping windows containing a 1-based VCF
    position.
    """

    k_lo = max(
        0,
        math.ceil(
            (pos - window_size) / window_step
        ),
    )

    k_hi = min(
        n_windows - 1,
        (pos - 1) // window_step,
    )

    if k_lo > k_hi:
        return range(0)

    return range(
        k_lo,
        k_hi + 1,
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate windowed Watterson's theta from "
            "a SnakePop all-sites VCF."
        )
    )

    parser.add_argument(
        "--vcf",
        required=True,
    )

    parser.add_argument(
        "--metadata",
        required=True,
    )

    parser.add_argument(
        "--chrom",
        required=True,
    )

    parser.add_argument(
        "--sample-col",
        required=True,
    )

    parser.add_argument(
        "--population-col",
        required=True,
    )

    parser.add_argument(
        "--populations",
        default="all",
    )

    parser.add_argument(
        "--window-size",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--window-step",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--min-called-fraction",
        type=float,
        default=0.8,
    )

    parser.add_argument(
        "--min-callable-sites",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--ploidy",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--snps-only",
        action="store_true",
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    ###########################################################################
    # Validate arguments
    ###########################################################################

    if args.window_size <= 0:
        raise ValueError(
            "--window-size must be positive"
        )

    if args.window_step <= 0:
        raise ValueError(
            "--window-step must be positive"
        )

    if not 0 < args.min_called_fraction <= 1:
        raise ValueError(
            "--min-called-fraction must be in (0, 1]"
        )

    if args.min_callable_sites < 1:
        raise ValueError(
            "--min-callable-sites must be >= 1"
        )

    if args.ploidy < 1:
        raise ValueError(
            "--ploidy must be >= 1"
        )

    ###########################################################################
    # Metadata
    ###########################################################################

    meta = pd.read_csv(
        args.metadata,
        sep="\t",
        dtype=str,
    )

    if args.sample_col not in meta.columns:
        raise ValueError(
            f"Missing sample column in metadata: "
            f"{args.sample_col}"
        )

    if args.population_col not in meta.columns:
        raise ValueError(
            f"Missing population column in metadata: "
            f"{args.population_col}"
        )

    meta = (
        meta[
            [
                args.sample_col,
                args.population_col,
            ]
        ]
        .dropna()
        .copy()
    )

    meta[args.population_col] = (
        meta[args.population_col]
        .map(clean_pop_name)
    )

    ###########################################################################
    # VCF
    ###########################################################################

    vcf = pysam.VariantFile(
        args.vcf
    )

    vcf_samples = set(
        vcf.header.samples
    )

    if not vcf_samples:
        raise ValueError(
            "VCF contains no samples"
        )

    ###########################################################################
    # Populations
    ###########################################################################

    if args.populations == "all":

        populations = sorted(
            meta[
                args.population_col
            ].unique()
        )

    else:

        populations = [
            clean_pop_name(x.strip())
            for x
            in args.populations.split(",")
            if x.strip()
        ]

    pop_samples = {}

    for population in populations:

        samples = meta.loc[
            meta[args.population_col]
            == population,
            args.sample_col,
        ].tolist()

        samples = [
            sample
            for sample in samples
            if sample in vcf_samples
        ]

        if not samples:
            raise ValueError(
                f"No VCF samples found for population: "
                f"{population}"
            )

        pop_samples[
            population
        ] = samples

    ###########################################################################
    # Chromosome
    ###########################################################################

    chrom = args.chrom

    if chrom not in vcf.header.contigs:
        raise ValueError(
            f"Chromosome {chrom} not found "
            f"in VCF header"
        )

    chrom_length = (
        vcf.header.contigs[
            chrom
        ].length
    )

    if chrom_length is None:
        raise ValueError(
            f"No contig length found in "
            f"VCF header for {chrom}"
        )

    ###########################################################################
    # Windows
    ###########################################################################

    starts = np.arange(
        1,
        chrom_length + 1,
        args.window_step,
        dtype=np.int64,
    )

    ends = np.minimum(
        starts
        + args.window_size
        - 1,
        chrom_length,
    )

    n_windows = len(
        starts
    )

    ###########################################################################
    # Per-population arrays
    ###########################################################################

    stats = {}
    thresholds = {}

    for population, samples in pop_samples.items():

        max_chromosomes = (
            args.ploidy
            * len(samples)
        )

        thresholds[
            population
        ] = max(
            2,
            math.ceil(
                max_chromosomes
                * args.min_called_fraction
            ),
        )

        stats[
            population
        ] = {

            "callable_sites":
                np.zeros(
                    n_windows,
                    dtype=np.int64,
                ),

            "segregating_sites":
                np.zeros(
                    n_windows,
                    dtype=np.int64,
                ),

            "harmonic_denominator":
                np.zeros(
                    n_windows,
                    dtype=float,
                ),

            "called_chromosome_sum":
                np.zeros(
                    n_windows,
                    dtype=np.int64,
                ),
        }

    ###########################################################################
    # Cache harmonic numbers
    ###########################################################################

    max_n = (
        args.ploidy
        * len(vcf_samples)
    )

    a1_cache = {
        n: harmonic_a1(n)
        for n in range(
            2,
            max_n + 1,
        )
    }

    ###########################################################################
    # Process sites
    ###########################################################################

    for record in vcf.fetch(
        chrom
    ):

        pos = record.pos

        if (
            args.snps_only
            and not is_nucleotide_site(record)
        ):
            continue

        indices = list(
            window_indices(
                pos=pos,
                n_windows=n_windows,
                window_size=args.window_size,
                window_step=args.window_step,
            )
        )

        if not indices:
            continue

        for (
            population,
            samples,
        ) in pop_samples.items():

            alleles = called_alleles(
                record,
                samples,
            )

            n_called = len(
                alleles
            )

            if (
                n_called
                < thresholds[
                    population
                ]
            ):
                continue

            a1 = a1_cache.get(
                n_called
            )

            if a1 is None:

                a1 = harmonic_a1(
                    n_called
                )

                a1_cache[
                    n_called
                ] = a1

            if a1 <= 0:
                continue

            # A site is segregating when more than one
            # allele is observed in the population.
            segregating = (
                len(set(alleles))
                > 1
            )

            s = stats[
                population
            ]

            for i in indices:

                s[
                    "callable_sites"
                ][i] += 1

                s[
                    "harmonic_denominator"
                ][i] += a1

                s[
                    "called_chromosome_sum"
                ][i] += n_called

                if segregating:

                    s[
                        "segregating_sites"
                    ][i] += 1

    ###########################################################################
    # Output
    ###########################################################################

    rows = []

    for (
        population,
        samples,
    ) in pop_samples.items():

        max_chromosomes = (
            args.ploidy
            * len(samples)
        )

        s = stats[
            population
        ]

        for i, (
            start,
            end,
        ) in enumerate(
            zip(
                starts,
                ends,
            )
        ):

            callable_sites = int(
                s[
                    "callable_sites"
                ][i]
            )

            segregating_sites = int(
                s[
                    "segregating_sites"
                ][i]
            )

            denominator = float(
                s[
                    "harmonic_denominator"
                ][i]
            )

            if (
                callable_sites
                >= args.min_callable_sites
                and denominator > 0
            ):

                theta_w = (
                    segregating_sites
                    / denominator
                )

                mean_called = (
                    s[
                        "called_chromosome_sum"
                    ][i]
                    / callable_sites
                )

            else:

                theta_w = np.nan
                mean_called = np.nan

            rows.append(
                {
                    "chrom":
                        chrom,

                    "start":
                        int(start),

                    "end":
                        int(end),

                    "population":
                        population,

                    "n_samples":
                        len(samples),

                    "max_chromosomes":
                        max_chromosomes,

                    "min_called_chromosomes":
                        thresholds[
                            population
                        ],

                    "callable_sites":
                        callable_sites,

                    "segregating_sites":
                        segregating_sites,

                    "mean_called_chromosomes":
                        mean_called,

                    "harmonic_denominator":
                        denominator,

                    "theta_w":
                        theta_w,
                }
            )

    out = pd.DataFrame(
        rows
    )

    os.makedirs(
        os.path.dirname(
            args.output
        )
        or ".",
        exist_ok=True,
    )

    out.to_csv(
        args.output,
        sep="\t",
        index=False,
        na_rep="NA",
    )


if __name__ == "__main__":
    main()
