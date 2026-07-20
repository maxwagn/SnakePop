#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import matplotlib.pyplot as plt


def clean_pop_name(x):
    return str(x).replace(" ", "_").replace("/", "_").replace(";", "_")


def read_regions(path):
    if path is None or path == "none":
        return None

    df = pd.read_csv(path, sep="\t")
    required = {"scaffold", "start", "end"}

    if not required.issubset(df.columns):
        raise ValueError(f"Region file missing columns {required}: {path}")

    regions = {}
    for _, row in df.iterrows():
        chrom = str(row["scaffold"])
        regions.setdefault(chrom, []).append((int(row["start"]), int(row["end"])))

    for chrom in regions:
        regions[chrom].sort()

    return regions


def in_regions(chrom, pos, regions):
    if regions is None:
        return True

    if chrom not in regions:
        return False

    for start, end in regions[chrom]:
        if start <= pos <= end:
            return True
        if pos < start:
            return False

    return False


def folded_count(alt_count, n_called):
    return min(alt_count, n_called - alt_count)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--vcf", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--sample-col", required=True)
    ap.add_argument("--population-col", required=True)
    ap.add_argument("--pop1", required=True)
    ap.add_argument("--pop2", default=None)
    ap.add_argument("--mode", choices=["1d", "2d"], required=True)
    ap.add_argument("--regions", default="none")
    ap.add_argument("--polarisation", choices=["folded"], default="folded")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--formats", default="pdf,svg,png")

    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    if args.sample_col not in meta.columns:
        raise ValueError(f"Missing sample column: {args.sample_col}")

    if args.population_col not in meta.columns:
        raise ValueError(f"Missing population column: {args.population_col}")

    meta[args.population_col] = meta[args.population_col].map(clean_pop_name)

    vcf = VCF(args.vcf)
    samples = list(vcf.samples)

    pop1_samples = set(meta.loc[meta[args.population_col] == args.pop1, args.sample_col])
    pop1_idx = [i for i, s in enumerate(samples) if s in pop1_samples]

    if len(pop1_idx) == 0:
        raise ValueError(f"No VCF samples found for population: {args.pop1}")

    if args.mode == "2d":
        if args.pop2 is None:
            raise ValueError("--pop2 is required for 2D SFS")

        pop2_samples = set(meta.loc[meta[args.population_col] == args.pop2, args.sample_col])
        pop2_idx = [i for i, s in enumerate(samples) if s in pop2_samples]

        if len(pop2_idx) == 0:
            raise ValueError(f"No VCF samples found for population: {args.pop2}")
    else:
        pop2_idx = []

    regions = read_regions(args.regions)

    n1 = 2 * len(pop1_idx)

    if args.mode == "1d":
        sfs = np.zeros(n1 + 1, dtype=np.int64)
    else:
        n2 = 2 * len(pop2_idx)
        sfs = np.zeros((n1 + 1, n2 + 1), dtype=np.int64)

    n_sites_seen = 0
    n_sites_used = 0

    for var in vcf:
        if len(var.ALT) != 1:
            continue

        if not in_regions(var.CHROM, var.POS, regions):
            continue

        n_sites_seen += 1
        gt = np.array(var.genotypes)[:, :2]

        def count_alt(indices):
            sub = gt[indices, :]
            called = sub[sub >= 0]

            if len(called) == 0:
                return None, None

            alt = int(np.sum(called == 1))
            n_called = int(len(called))

            return alt, n_called

        alt1, called1 = count_alt(pop1_idx)

        if called1 is None or called1 != n1:
            continue

        c1 = folded_count(alt1, called1)

        if args.mode == "1d":
            sfs[c1] += 1
            n_sites_used += 1

        else:
            alt2, called2 = count_alt(pop2_idx)

            if called2 is None or called2 != n2:
                continue

            c2 = folded_count(alt2, called2)
            sfs[c1, c2] += 1
            n_sites_used += 1

    meta_out = pd.DataFrame([{
        "mode": args.mode,
        "polarisation": args.polarisation,
        "pop1": args.pop1,
        "pop2": args.pop2 if args.pop2 else "",
        "n_pop1_samples": len(pop1_idx),
        "n_pop1_chromosomes": n1,
        "n_pop2_samples": len(pop2_idx) if args.mode == "2d" else "",
        "n_pop2_chromosomes": n2 if args.mode == "2d" else "",
        "n_sites_seen_in_region": n_sites_seen,
        "n_sites_used_complete_data": n_sites_used,
    }])

    meta_out.to_csv(args.out_prefix + ".summary.tsv", sep="\t", index=False)

    if args.mode == "1d":
        out = pd.DataFrame({
            "allele_count": np.arange(len(sfs)),
            "n_sites": sfs,
        })

        out.to_csv(args.out_prefix + ".sfs.tsv", sep="\t", index=False)

        for fmt in args.formats.split(","):
            fmt = fmt.strip()
            plt.figure(figsize=(7, 4))
            plt.bar(out["allele_count"], out["n_sites"])
            plt.xlabel("Folded allele count")
            plt.ylabel("Number of SNPs")
            plt.title(f"1D folded SFS: {args.pop1}")
            plt.tight_layout()
            plt.savefig(args.out_prefix + f".{fmt}", dpi=300)
            plt.close()

    else:
        out = pd.DataFrame(sfs)
        out.index.name = args.pop1
        out.columns.name = args.pop2
        out.to_csv(args.out_prefix + ".sfs.tsv", sep="\t")

        for fmt in args.formats.split(","):
            fmt = fmt.strip()
            plt.figure(figsize=(6, 5))
            plt.imshow(np.log10(sfs + 1), origin="lower", aspect="auto")
            plt.colorbar(label="log10(SNPs + 1)")
            plt.xlabel(args.pop2)
            plt.ylabel(args.pop1)
            plt.title(f"2D folded SFS: {args.pop1} vs {args.pop2}")
            plt.tight_layout()
            plt.savefig(args.out_prefix + f".{fmt}", dpi=300)
            plt.close()


if __name__ == "__main__":
    main()
