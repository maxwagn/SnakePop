###############################################################################
# SnakePop: population statistics
###############################################################################

import pandas as pd

###############################################################################
# Config
###############################################################################

REF_NAME = config["ref"]["name"]
CALLSET_ID = config["callset"]["id"]

IND_FILTER_ID = config["variant_calling"]["ind_filter_id"]
SITE_FILTER_ID = config["variant_calling"]["site_filter_id"]

VC_ROOT = f"results/variants/{CALLSET_ID}_{REF_NAME}"
POP_ROOT = f"results/popstats/{CALLSET_ID}_{REF_NAME}"
PCA_ROOT = f"{POP_ROOT}/pca"
REPORT_ROOT = f"reports/popstats/{CALLSET_ID}_{REF_NAME}/pca"

CHROM_LIST = config["ref"]["fasta"] + ".chromosomes.txt"
CHROMOSOMES = [line.strip() for line in open(CHROM_LIST) if line.strip()]

SAMPLE_TABLE = config["sample_table"]
SAMPLE_COL = config.get("sample_id_column", "id")

RES = config.get("resources", {})
PCA_RES = RES.get("pca", {})
PCA_CFG = config.get("popstats", {}).get("pca", {})

N_PCS = PCA_CFG.get("n_pcs", 10)
MAF = PCA_CFG.get("maf", 0.05)
MIND = PCA_CFG.get("mind", 0.2)
GENO = PCA_CFG.get("geno", 0.2)
COLOR_BY = PCA_CFG.get("color_by", "species")
LABEL_SAMPLES = PCA_CFG.get("label_samples", False)
LD_PRUNE = PCA_CFG.get("ld_prune", False)

LD = PCA_CFG.get("indep_pairwise", {})
LD_WINDOW = LD.get("window_kb", 50)
LD_STEP = LD.get("step", 5)
LD_R2 = LD.get("r2", 0.2)

PLINK_CHR_SET = PCA_CFG.get("plink_chr_set", "23 no-xy no-mt")

BIALLELIC_SNPS = expand(
    f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
    chrom=CHROMOSOMES,
)


###############################################################################
# Public target
###############################################################################

rule pca:
    input:
        eigenvec = f"{PCA_ROOT}/plink_pca.eigenvec",
        eigenval = f"{PCA_ROOT}/plink_pca.eigenval",
        table = f"{REPORT_ROOT}/pca_scores.tsv",
        png = f"{REPORT_ROOT}/pca_PC1_PC2.png",
        svg = f"{REPORT_ROOT}/pca_PC1_PC2.svg"


###############################################################################
# PCA workflow
###############################################################################

rule concat_biallelic_snps:
    input:
        vcfs = BIALLELIC_SNPS
    output:
        vcf = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz",
        tbi = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz.tbi"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    shell:
        r"""
        bcftools concat \
          --threads {threads} \
          -Oz \
          -o {output.vcf} \
          {input.vcfs}

        bcftools index --tbi -f {output.vcf}
        """


rule make_plink_pgen:
    input:
        vcf = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz",
        tbi = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz.tbi"
    output:
        pgen = f"{PCA_ROOT}/plink_raw.pgen",
        pvar = f"{PCA_ROOT}/plink_raw.pvar",
        psam = f"{PCA_ROOT}/plink_raw.psam"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    params:
        chr_set = PLINK_CHR_SET
    shell:
        r"""
        plink2 \
          --threads {threads} \
          --memory {resources.mem_mb} \
          --vcf {input.vcf} \
          --double-id \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --set-all-var-ids @:#:\$r:\$a \
          --make-pgen \
          --out {PCA_ROOT}/plink_raw
        """


rule qc_plink_pgen:
    input:
        pgen = f"{PCA_ROOT}/plink_raw.pgen",
        pvar = f"{PCA_ROOT}/plink_raw.pvar",
        psam = f"{PCA_ROOT}/plink_raw.psam"
    output:
        pgen = f"{PCA_ROOT}/plink_qc.pgen",
        pvar = f"{PCA_ROOT}/plink_qc.pvar",
        psam = f"{PCA_ROOT}/plink_qc.psam"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    params:
        chr_set = PLINK_CHR_SET
    shell:
        r"""
        plink2 \
          --threads {threads} \
          --memory {resources.mem_mb} \
          --pfile {PCA_ROOT}/plink_raw \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --maf {MAF} \
          --geno {GENO} \
          --mind {MIND} \
          --make-pgen \
          --out {PCA_ROOT}/plink_qc
        """


rule ld_prune:
    input:
        pgen = f"{PCA_ROOT}/plink_qc.pgen",
        pvar = f"{PCA_ROOT}/plink_qc.pvar",
        psam = f"{PCA_ROOT}/plink_qc.psam"
    output:
        prune_in = f"{PCA_ROOT}/plink_qc.prune.in",
        prune_out = f"{PCA_ROOT}/plink_qc.prune.out"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    params:
        chr_set = PLINK_CHR_SET
    shell:
        r"""
        plink2 \
          --threads {threads} \
          --memory {resources.mem_mb} \
          --pfile {PCA_ROOT}/plink_qc \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --indep-pairwise {LD_WINDOW} {LD_STEP} {LD_R2} \
          --out {PCA_ROOT}/plink_qc
        """


rule calc_allele_freq:
    input:
        pgen = f"{PCA_ROOT}/plink_qc.pgen",
        pvar = f"{PCA_ROOT}/plink_qc.pvar",
        psam = f"{PCA_ROOT}/plink_qc.psam"
    output:
        afreq = f"{PCA_ROOT}/plink_qc.afreq"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    params:
        chr_set = PLINK_CHR_SET
    shell:
        r"""
        plink2 \
          --threads {threads} \
          --memory {resources.mem_mb} \
          --pfile {PCA_ROOT}/plink_qc \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --freq \
          --out {PCA_ROOT}/plink_qc
        """


def pca_inputs(wildcards):
    inputs = {
        "pgen": f"{PCA_ROOT}/plink_qc.pgen",
        "pvar": f"{PCA_ROOT}/plink_qc.pvar",
        "psam": f"{PCA_ROOT}/plink_qc.psam",
        "afreq": f"{PCA_ROOT}/plink_qc.afreq",
    }

    if LD_PRUNE:
        inputs["prune_in"] = f"{PCA_ROOT}/plink_qc.prune.in"

    return inputs

def pca_extract_param(wildcards, input):
    if LD_PRUNE:
        return f"--extract {input.prune_in}"
    return ""


rule run_pca:
    input:
        unpack(pca_inputs)
    output:
        eigenvec = f"{PCA_ROOT}/plink_pca.eigenvec",
        eigenval = f"{PCA_ROOT}/plink_pca.eigenval"
    threads: PCA_RES.get("threads", 8)
    resources:
        mem_mb = PCA_RES.get("mem_mb", 32000),
        walltime = PCA_RES.get("walltime", 8)
    params:
        chr_set = PLINK_CHR_SET,
        extract = pca_extract_param
    shell:
        r"""
        plink2 \
          --threads {threads} \
          --memory {resources.mem_mb} \
          --pfile {PCA_ROOT}/plink_qc \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --read-freq {input.afreq} \
          {params.extract} \
          --pca {N_PCS} \
          --out {PCA_ROOT}/plink_pca
        """

rule plot_pca:
    input:
        eigenvec = f"{PCA_ROOT}/plink_pca.eigenvec",
        eigenval = f"{PCA_ROOT}/plink_pca.eigenval",
        metadata = SAMPLE_TABLE
    output:
        table = f"{REPORT_ROOT}/pca_scores.tsv",
        png = f"{REPORT_ROOT}/pca_PC1_PC2.png",
        svg = f"{REPORT_ROOT}/pca_PC1_PC2.svg"
    threads: 1
    resources:
        mem_mb = 4000,
        walltime = 1
    params:
        color_by = COLOR_BY,
        label_samples = LABEL_SAMPLES
    run:
        import matplotlib.pyplot as plt

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        pcs = pd.read_csv(input.eigenvec, sep=r"\s+", dtype=str)
        pcs.columns = [c.replace("#", "") for c in pcs.columns]

        eigenvals = pd.read_csv(input.eigenval, header=None)[0].astype(float)

        pcs = pcs.rename(columns={"IID": SAMPLE_COL})

        pc_cols = [c for c in pcs.columns if c.startswith("PC")]
        for col in pc_cols:
            pcs[col] = pcs[col].astype(float)

        df = pcs.merge(meta, on=SAMPLE_COL, how="left")
        df.to_csv(output.table, sep="\t", index=False)

        pc1_var = 100 * eigenvals.iloc[0] / eigenvals.sum()
        pc2_var = 100 * eigenvals.iloc[1] / eigenvals.sum()

        plt.figure(figsize=(6, 5))

        if params.color_by in df.columns:
            for group, sub in df.groupby(params.color_by, dropna=False):
                plt.scatter(
                    sub["PC1"],
                    sub["PC2"],
                    label=str(group),
                    s=45,
                    alpha=0.85,
                    edgecolors="black",
                    linewidths=0.4,
                )

            plt.legend(
                title=params.color_by,
                frameon=False,
                fontsize=8,
                title_fontsize=9,
            )
        else:
            plt.scatter(
                df["PC1"],
                df["PC2"],
                s=45,
                alpha=0.85,
                edgecolors="black",
                linewidths=0.4,
            )

        if params.label_samples:
            for _, row in df.iterrows():
                plt.text(
                    row["PC1"],
                    row["PC2"],
                    str(row[SAMPLE_COL]),
                    fontsize=7,
                    ha="left",
                    va="bottom",
                )

        plt.xlabel(f"PC1 ({pc1_var:.2f}%)")
        plt.ylabel(f"PC2 ({pc2_var:.2f}%)")
        plt.title(f"PCA colored by {params.color_by}")
        plt.tight_layout()

        plt.savefig(output.png, dpi=300)
        plt.savefig(output.svg)
        plt.close()


###############################################################################
# popgenWindows
###############################################################################

PGW_CFG = config["popstats"].get("popgenwindows", {})

PGW_ROOT = f"{POP_ROOT}/popgenwindows"
PGW_GENO_ROOT = f"{PGW_ROOT}/geno"
PGW_POP_ROOT = f"{POP_ROOT}/populations"

PARSE_VCF = "bin/parseVCF.py"
POPGENWINDOWS = "bin/popgenWindows.py"

PGW_WINDOW_SIZE = PGW_CFG.get("window_size", 60000)
PGW_WINDOW_STEP = PGW_CFG.get("window_step", 30000)
PGW_WINDOW_MIN_SITES = PGW_CFG.get("window_min_sites", 100)
PGW_FORMAT = PGW_CFG.get("input_format", "phased")

POP_COL = config["popstats"].get("population_column", "morphology")

PGW_LABEL = f"w{PGW_WINDOW_SIZE}.m{PGW_WINDOW_MIN_SITES}.s{PGW_WINDOW_STEP}"

PGW_RES = RES.get("popgenwindows", {})
PARSE_RES = RES.get("parse_vcf", {})


def clean_pop_name(x):
    return str(x).replace(" ", "_").replace("/", "_").replace(";", "_")


meta = pd.read_csv(SAMPLE_TABLE, sep="\t", dtype=str)

if POP_COL not in meta.columns:
    raise ValueError(f"Population column not found in metadata: {POP_COL}")

POPS = sorted(clean_pop_name(x) for x in meta[POP_COL].dropna().unique())


def population_string():
    return " ".join([f"-p {pop}" for pop in POPS])


PGW_CHROM_TARGETS = expand(
    f"{PGW_ROOT}/{{chrom}}.{PGW_LABEL}.Fst.Dxy.pi.csv.gz",
    chrom=CHROMOSOMES,
)

PGW_MERGED_TARGET = f"{PGW_ROOT}/popgenwindows.{PGW_LABEL}.Fst.Dxy.pi.dA.csv.gz"


rule make_popgenwindows_population_file:
    input:
        metadata = SAMPLE_TABLE
    output:
        popfile = f"{PGW_POP_ROOT}/popgenwindows_populations.tsv"
    run:
        import os
        import pandas as pd

        df = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in df.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if POP_COL not in df.columns:
            raise ValueError(f"Missing population column in metadata: {POP_COL}")

        df = df[[SAMPLE_COL, POP_COL]].dropna()
        df[POP_COL] = df[POP_COL].map(clean_pop_name)

        os.makedirs(PGW_POP_ROOT, exist_ok=True)

        df.to_csv(
            output.popfile,
            sep="\t",
            header=False,
            index=False,
        )


rule prep_geno:
    input:
        vcf = f"{VC_ROOT}/vcf/all_sites.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz"
    output:
        geno = f"{PGW_GENO_ROOT}/all_sites.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.geno.gz"
    threads: PARSE_RES.get("threads", 1)
    resources:
        mem_mb = PARSE_RES.get("mem_mb", 4000),
        walltime = PARSE_RES.get("walltime", 2)
    shell:
        r"""
        mkdir -p {PGW_GENO_ROOT}

        python {PARSE_VCF} \
          -i {input.vcf} \
        | bgzip -c > {output.geno}
        """


rule popgenwindows_chrom:
    input:
        geno = f"{PGW_GENO_ROOT}/all_sites.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.geno.gz",
        popfile = f"{PGW_POP_ROOT}/popgenwindows_populations.tsv"
    output:
        stats = f"{PGW_ROOT}/{{chrom}}.{PGW_LABEL}.Fst.Dxy.pi.csv.gz"
    threads: PGW_RES.get("threads", 8)
    resources:
        mem_mb = PGW_RES.get("mem_mb", 16000),
        walltime = PGW_RES.get("walltime", 2)
    params:
        populations = population_string()
    shell:
        r"""
        mkdir -p {PGW_ROOT}

        python {POPGENWINDOWS} \
          -g {input.geno} \
          -o {output.stats} \
          -f {PGW_FORMAT} \
          -w {PGW_WINDOW_SIZE} \
          -m {PGW_WINDOW_MIN_SITES} \
          -s {PGW_WINDOW_STEP} \
          {params.populations} \
          --popsFile {input.popfile} \
          -T {threads} \
          --writeFailedWindows
        """


rule merge_popgenwindows:
    input:
        stats = PGW_CHROM_TARGETS
    output:
        merged = PGW_MERGED_TARGET
    run:
        import os
        import pandas as pd

        dfs = []
        reference_header = None

        for fn in input.stats:
            df = pd.read_csv(fn, compression="gzip")

            if reference_header is None:
                reference_header = list(df.columns)

            elif set(reference_header) != set(df.columns):
                replace_dic = {}

                for col in set(df.columns) - set(reference_header):
                    parts = col.split("_")
                    if len(parts) == 3:
                        entry, pop1, pop2 = parts
                        swapped = f"{entry}_{pop2}_{pop1}"

                        if swapped in reference_header:
                            replace_dic[col] = swapped

                df.rename(columns=replace_dic, inplace=True)

                if set(reference_header) != set(df.columns):
                    raise ValueError(f"Header mismatch after repair for file: {fn}")

            df = df[reference_header]
            dfs.append(df)

        merged = pd.concat(dfs, axis=0)

        dxy_cols = [col for col in merged.columns if col.startswith("dxy_")]


        for col in dxy_cols:
            pair = col.replace("dxy_", "")

            matched = None

            for pop1 in POPS:
                for pop2 in POPS:
                    if pop1 == pop2:
                        continue
                    if pair == f"{pop1}_{pop2}":
                        matched = (pop1, pop2)
                        break
                if matched:
                    break

            if matched is None:
                print(f"WARNING: Could not parse dxy column for dA: {col}")
                continue

            pop1, pop2 = matched
            pi1 = f"pi_{pop1}"
            pi2 = f"pi_{pop2}"

            if pi1 in merged.columns and pi2 in merged.columns:
                merged[f"dA_{pop1}_{pop2}"] = (
                        merged[col] - (merged[pi1] + merged[pi2]) / 2)
        
        os.makedirs(PGW_ROOT, exist_ok=True)
        merged.to_csv(output.merged, sep="\t", compression="gzip",index=False,)


rule summarize_popgenwindows:
    input:
        merged = PGW_MERGED_TARGET
    output:
        summary = PGW_MERGED_TARGET.replace(".csv.gz", ".summary.tsv")
    run:
        import pandas as pd

        df = pd.read_csv(input.merged, sep="\t")

        skip_cols = {
            "scaffold",
            "start",
            "end",
            "mid",
            "sites",
        }

        stat_cols = [c for c in df.columns if c not in skip_cols]

        rows = []

        for col in stat_cols:
            x = pd.to_numeric(df[col], errors="coerce").dropna()

            if len(x) == 0:
                continue

            rows.append({
                "statistic": col,
                "n_windows": len(x),
                "mean": x.mean(),
                "median": x.median(),
                "sd": x.std(),
                "min": x.min(),
                "p05": x.quantile(0.05),
                "p95": x.quantile(0.95),
                "max": x.max(),
            })

        summary = pd.DataFrame(rows)

        summary.sort_values(
            "statistic",
            inplace=True,
        )

        summary.to_csv(
            output.summary,
            sep="\t",
            index=False,
        )





###############################################################################
# Manhattan-style popgenWindows reports
###############################################################################

MAN_CFG = config["popstats"].get("manhattan", {})

MAN_ROOT = f"{POP_ROOT}/popgenwindows/manhattan"

MAN_PAIRS = MAN_CFG.get("pairs", [])
MAN_METRICS = MAN_CFG.get("metrics", ["Fst", "dxy", "dA"])
MAN_INCLUDE_PI = MAN_CFG.get("include_pi", True)
MAN_OUTLIER_Q = MAN_CFG.get("outlier_quantile", 0.99)
MAN_COLORS = MAN_CFG.get("colors", {})
MAN_YLIM = MAN_CFG.get("ylim", {})


# Optional: compare merged Manhattan outlier regions against existing WinPCA HTMLs.
# If false, Manhattan runs without requiring WinPCA outputs.
MAN_COMPARE_WINPCA = MAN_CFG.get("compare_winpca", False)
#MAN_WINPCA_ROOT = MAN_CFG.get("winpca_root", f"{POP_ROOT}/winpca")

CAND_CFG = MAN_CFG.get("candidate_regions", {})
CAND_PRIMARY_METRIC = CAND_CFG.get("primary_metric", "Fst")
CAND_MERGE_GAP = CAND_CFG.get("merge_gap_bp", PGW_WINDOW_STEP)
CAND_MIN_WINDOWS = CAND_CFG.get("min_windows", 2)


def manhattan_pair_names():
    return [f"{p['pop1']}_{p['pop2']}" for p in MAN_PAIRS]


def manhattan_targets():
    targets = []

    for p in MAN_PAIRS:
        pair = f"{p['pop1']}_{p['pop2']}"

        targets.extend([
            f"{MAN_ROOT}/{pair}/{pair}.manhattan.pdf",
            f"{MAN_ROOT}/{pair}/{pair}.manhattan.png",
            f"{MAN_ROOT}/{pair}/{pair}.outlier_windows.tsv",
            f"{MAN_ROOT}/{pair}/{pair}.summary.html",
        ])

        if MAN_COMPARE_WINPCA:
            targets.extend([
                f"{MAN_ROOT}/{pair}/{pair}.candidate_regions.tsv",
                f"{MAN_ROOT}/{pair}/{pair}.candidate_regions.html",
            ])

    targets.append(f"{MAN_ROOT}/manhattan_index.html")
    return targets

rule plot_pair_manhattan:
    input:
        merged = PGW_MERGED_TARGET
    output:
        pdf = f"{MAN_ROOT}/{{pair}}/{{pair}}.manhattan.pdf",
        png = f"{MAN_ROOT}/{{pair}}/{{pair}}.manhattan.png",
        outliers = f"{MAN_ROOT}/{{pair}}/{{pair}}.outlier_windows.tsv",
        html = f"{MAN_ROOT}/{{pair}}/{{pair}}.summary.html"
    params:
        outlier_q = MAN_OUTLIER_Q,
        metrics = MAN_METRICS,
        include_pi = MAN_INCLUDE_PI,
        colors = MAN_COLORS,
        ylim = MAN_YLIM
    run:
        import os
        import numpy as np
        import pandas as pd
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        os.makedirs(os.path.dirname(output.pdf), exist_ok=True)

        pair = wildcards.pair

        found = None
        for p in MAN_PAIRS:
            if f"{p['pop1']}_{p['pop2']}" == pair:
                found = p
                break

        if found is None:
            raise ValueError(f"Could not find pair in config: {pair}")

        pop1 = found["pop1"]
        pop2 = found["pop2"]

        df = pd.read_csv(input.merged, sep="\t")

        required = {"scaffold", "start", "end", "mid", "sites"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        chrom_order = list(CHROMOSOMES)
        df = df[df["scaffold"].isin(chrom_order)].copy()
        df["scaffold"] = pd.Categorical(
            df["scaffold"],
            categories=chrom_order,
            ordered=True,
        )
        df = df.sort_values(["scaffold", "start"]).reset_index(drop=True)

        chrom_lengths = (
            df.groupby("scaffold", observed=True)["end"]
            .max()
            .reindex(chrom_order)
            .dropna()
        )

        offsets = {}
        chrom_centers = {}
        current = 0

        for chrom, length in chrom_lengths.items():
            offsets[str(chrom)] = current
            chrom_centers[str(chrom)] = current + length / 2
            current += length

        total_len = current

        df["genome_pos"] = df.apply(
            lambda r: offsets[str(r["scaffold"])] + float(r["mid"]),
            axis=1,
        )

        metric_cols = []

        for metric in params.metrics:
            col1 = f"{metric}_{pop1}_{pop2}"
            col2 = f"{metric}_{pop2}_{pop1}"

            if col1 in df.columns:
                metric_cols.append((metric, col1))
            elif col2 in df.columns:
                metric_cols.append((metric, col2))
            else:
                raise ValueError(
                    f"Missing column for {metric}: expected {col1} or {col2}"
                )

        pi_cols = []

        if params.include_pi:
            for pop in [pop1, pop2]:
                col = f"pi_{pop}"
                if col in df.columns:
                    pi_cols.append((f"pi_{pop}", col, pop))
                else:
                    raise ValueError(f"Missing pi column: {col}")

        # Separate pi panels: one panel per population.
        plot_cols = metric_cols + [(label, col) for label, col, pop in pi_cols]

        outlier_rows = []

        for label, col in plot_cols:
            y = pd.to_numeric(df[col], errors="coerce")
            threshold = y.quantile(params.outlier_q)

            sub = df.loc[
                y >= threshold,
                ["scaffold", "start", "end", "mid", "sites"]
            ].copy()

            sub["pair"] = pair
            sub["metric"] = label
            sub["value"] = y.loc[sub.index]
            sub["threshold_quantile"] = params.outlier_q
            sub["threshold_value"] = threshold

            outlier_rows.append(sub)

        outliers = pd.concat(outlier_rows, axis=0) if outlier_rows else pd.DataFrame()
        outliers.to_csv(output.outliers, sep="\t", index=False)

        # Plot settings: A4 landscape, editable text, rasterized points.
        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42
        mpl.rcParams["font.family"] = "DejaVu Sans"

        n_panels = len(plot_cols)

        fig, axes = plt.subplots(
            n_panels,
            1,
            figsize=(11.69, 8.27),
            sharex=True,
            gridspec_kw={"height_ratios": [1] * n_panels},
        )

        if n_panels == 1:
            axes = [axes]

        chroms = list(chrom_lengths.index)

        band_colors = [
            "#cfe2f3", "#d9ead3", "#fff2cc", "#fce5cd",
            "#f4cccc", "#d9e2f3", "#d9ead3", "#fff2cc",
            "#fce5cd", "#eadcf8",
        ]

        def apply_ylim(ax, label):
            key = "pi" if label.startswith("pi_") else label

            if key in params.ylim:
                low, high = params.ylim[key]
                cur_low, cur_high = ax.get_ylim()

                if low is None:
                    low = cur_low
                if high is None:
                    high = cur_high

                ax.set_ylim(low, high)

        def decorate_chromosomes(ax):
            ymin, ymax = ax.get_ylim()
            yrange = ymax - ymin
            band_h = yrange * 0.13
            band_y0 = ymin

            for idx, chrom in enumerate(chroms):
                start = offsets[str(chrom)]
                end = offsets[str(chrom)] + chrom_lengths[chrom]
                center = chrom_centers[str(chrom)]

                ax.axvspan(
                    start,
                    end,
                    ymin=0,
                    ymax=1,
                    color="0.985" if idx % 2 == 0 else "white",
                    zorder=0,
                )

                ax.axvline(
                    start,
                    color="0.75",
                    linestyle=":",
                    linewidth=0.5,
                    zorder=1,
                )

                ax.axvspan(
                    start,
                    end,
                    ymin=0,
                    ymax=0.08,
                    color=band_colors[idx % len(band_colors)],
                    alpha=0.85,
                    zorder=1,
                )

                ax.text(
                    center,
                    band_y0 + band_h * 0.45,
                    str(chrom),
                    ha="center",
                    va="center",
                    fontsize=5.5,
                    fontweight="bold",
                    zorder=5,
                )

            ax.axvline(
                total_len,
                color="0.75",
                linestyle=":",
                linewidth=0.5,
                zorder=1,
            )

        for ax_i, (label, col) in enumerate(plot_cols):
            ax = axes[ax_i]

            y = pd.to_numeric(df[col], errors="coerce")
            threshold = y.quantile(params.outlier_q)
            is_outlier = y >= threshold

            ax.scatter(
                df.loc[~is_outlier, "genome_pos"],
                y.loc[~is_outlier],
                s=4,
                c="0.45",
                alpha=0.45,
                linewidths=0,
                rasterized=True,
                zorder=2,
            )

            ax.scatter(
                df.loc[is_outlier, "genome_pos"],
                y.loc[is_outlier],
                s=7,
                c="red",
                alpha=0.9,
                linewidths=0,
                rasterized=True,
                zorder=4,
            )

            ax.axhline(
                threshold,
                color="red",
                linestyle="--",
                linewidth=0.8,
                zorder=3,
            )

            apply_ylim(ax, label)
            decorate_chromosomes(ax)

            ax.set_xlim(0, total_len)

            pretty_label = label
            if label == "dxy":
                pretty_label = "dXY"
            elif label == "Fst":
                pretty_label = "FST"
            elif label.startswith("pi_"):
                pretty_label = f"π {label.replace('pi_', '')}"

            ax.set_ylabel(pretty_label, fontsize=8)

            ax.text(
                0.995,
                0.92,
                f"{pretty_label} outlier threshold ({params.outlier_q * 100:.1f}%): {threshold:.4g}",
                color="red",
                fontsize=6,
                ha="right",
                va="top",
                transform=ax.transAxes,
            )

            ax.tick_params(axis="both", labelsize=6)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        axes[-1].set_xticks([])
        axes[-1].set_xlabel("Genomic position", fontsize=8)

        fig.suptitle(
            f"{pop1} vs {pop2} – Genome scan",
            fontsize=12,
            fontweight="bold",
        )

        fig.tight_layout(rect=[0, 0, 1, 0.965])
        fig.savefig(output.pdf, dpi=300)
        fig.savefig(output.png, dpi=200)
        plt.close(fig)

        outlier_html = outliers.head(500).to_html(index=False, classes="outliers")

        html = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{pair} Manhattan summary</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    img {{
      max-width: 100%;
      border: 1px solid #ccc;
      background: white;
    }}
    table {{
      border-collapse: collapse;
      font-size: 12px;
      background: white;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 4px 8px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>{pair}</h1>

  <p>
    <b>PDF:</b> <a href="{os.path.basename(output.pdf)}">{os.path.basename(output.pdf)}</a><br>
    <b>PNG:</b> <a href="{os.path.basename(output.png)}">{os.path.basename(output.png)}</a><br>
    <b>Outlier table:</b> <a href="{os.path.basename(output.outliers)}">{os.path.basename(output.outliers)}</a>
  </p>

  <h2>Genome scan</h2>
  <img src="{os.path.basename(output.png)}">

  <h2>Outlier windows</h2>
  <p>Showing first 500 outlier windows. Full table is in the TSV file.</p>
  {outlier_html}
</body>
</html>
"""

        with open(output.html, "w") as out:
            out.write(html)


        import pandas as pd


        with open(output.html, "w") as out_html:
            out_html.write(html)


rule compare_manhattan_winpca:
    input:
        merged = PGW_MERGED_TARGET,
        pc1 = expand(
            f"{POP_ROOT}/winpca/{{chrom}}/{{chrom}}.pc_1.tsv.gz",
            chrom=CHROMOSOMES,
        ),
        stat = expand(
            f"{POP_ROOT}/winpca/{{chrom}}/{{chrom}}.stat.tsv.gz",
            chrom=CHROMOSOMES,
        ),
        metadata = SAMPLE_TABLE
    output:
        tsv = f"{MAN_ROOT}/{{pair}}/{{pair}}.candidate_regions.tsv",
        html = f"{MAN_ROOT}/{{pair}}/{{pair}}.candidate_regions.html"
    params:
        primary_metric = CAND_PRIMARY_METRIC,
        outlier_q = MAN_OUTLIER_Q,
        merge_gap = CAND_MERGE_GAP,
        min_windows = CAND_MIN_WINDOWS
    run:
        import os
        import numpy as np
        import pandas as pd

        os.makedirs(os.path.dirname(output.tsv), exist_ok=True)

        pair = wildcards.pair

        found = None
        for p in MAN_PAIRS:
            if f"{p['pop1']}_{p['pop2']}" == pair:
                found = p
                break

        if found is None:
            raise ValueError(f"Pair not found in config: {pair}")

        pop1 = found["pop1"]
        pop2 = found["pop2"]

        df = pd.read_csv(input.merged, sep="\t")

        def find_pair_col(metric):
            c1 = f"{metric}_{pop1}_{pop2}"
            c2 = f"{metric}_{pop2}_{pop1}"
            if c1 in df.columns:
                return c1
            if c2 in df.columns:
                return c2
            return None

        primary_col = find_pair_col(params.primary_metric)
        fst_col = find_pair_col("Fst")
        dxy_col = find_pair_col("dxy")
        da_col = find_pair_col("dA")

        if primary_col is None:
            raise ValueError(f"Could not find primary metric column: {params.primary_metric}")

        pi1_col = f"pi_{pop1}"
        pi2_col = f"pi_{pop2}"

        y = pd.to_numeric(df[primary_col], errors="coerce")
        threshold = y.quantile(params.outlier_q)

        hits = df.loc[
            y >= threshold,
            ["scaffold", "start", "end", "mid", "sites"]
        ].copy()

        hits["value"] = y.loc[hits.index]
        hits = hits.sort_values(["scaffold", "start"])

        regions = []

        for chrom, sub in hits.groupby("scaffold", sort=False):
            sub = sub.sort_values("start")

            cur_start = None
            cur_end = None
            cur_values = []
            cur_windows = 0

            for _, row in sub.iterrows():
                if cur_start is None:
                    cur_start = int(row["start"])
                    cur_end = int(row["end"])
                    cur_values = [row["value"]]
                    cur_windows = 1
                    continue

                if int(row["start"]) <= cur_end + params.merge_gap:
                    cur_end = max(cur_end, int(row["end"]))
                    cur_values.append(row["value"])
                    cur_windows += 1
                else:
                    if cur_windows >= params.min_windows:
                        regions.append({
                            "scaffold": chrom,
                            "start": cur_start,
                            "end": cur_end,
                            "n_outlier_windows": cur_windows,
                            "max_primary": np.nanmax(cur_values),
                            "mean_primary": np.nanmean(cur_values),
                        })

                    cur_start = int(row["start"])
                    cur_end = int(row["end"])
                    cur_values = [row["value"]]
                    cur_windows = 1

            if cur_start is not None and cur_windows >= params.min_windows:
                regions.append({
                    "scaffold": chrom,
                    "start": cur_start,
                    "end": cur_end,
                    "n_outlier_windows": cur_windows,
                    "max_primary": np.nanmax(cur_values),
                    "mean_primary": np.nanmean(cur_values),
                })

        regions = pd.DataFrame(regions)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if POP_COL not in meta.columns:
            raise ValueError(f"Missing population column in metadata: {POP_COL}")

        meta = meta[[SAMPLE_COL, POP_COL]].dropna().copy()
        meta[POP_COL] = meta[POP_COL].map(clean_pop_name)

        winpca_rows = []

        for pc1_file, stat_file in zip(input.pc1, input.stat):
            chrom = os.path.basename(os.path.dirname(pc1_file))

            pc = pd.read_csv(pc1_file, sep="\t", compression="gzip")
            st = pd.read_csv(stat_file, sep="\t", compression="gzip")

            sample_cols = [c for c in pc.columns if c != "pos"]

            long = pc.melt(
                id_vars="pos",
                value_vars=sample_cols,
                var_name=SAMPLE_COL,
                value_name="pc1",
            )

            long = long.merge(meta, on=SAMPLE_COL, how="left")
            long = long[long[POP_COL].isin([pop1, pop2])].copy()
            long["pc1"] = pd.to_numeric(long["pc1"], errors="coerce")

            means = (
                long.groupby(["pos", POP_COL])["pc1"]
                .mean()
                .reset_index()
                .pivot(index="pos", columns=POP_COL, values="pc1")
                .reset_index()
            )

            if pop1 not in means.columns or pop2 not in means.columns:
                continue

            means["abs_delta_pc1"] = (means[pop1] - means[pop2]).abs()
            means["scaffold"] = chrom

            keep_stat = st[["pos", "pc_1_ve", "n_var"]].copy()
            keep_stat["pc_1_ve"] = pd.to_numeric(keep_stat["pc_1_ve"], errors="coerce")
            keep_stat["n_var"] = pd.to_numeric(keep_stat["n_var"], errors="coerce")

            means = means.merge(keep_stat, on="pos", how="left")
            winpca_rows.append(means)

        if winpca_rows:
            winpca = pd.concat(winpca_rows, axis=0, ignore_index=True)
        else:
            winpca = pd.DataFrame(
                columns=["scaffold", "pos", "abs_delta_pc1", "pc_1_ve", "n_var"]
            )

        rows = []

        for _, reg in regions.iterrows():
            chrom = str(reg["scaffold"])
            start = int(reg["start"])
            end = int(reg["end"])

            sub = df[
                (df["scaffold"].astype(str) == chrom) &
                (df["mid"] >= start) &
                (df["mid"] <= end)
            ].copy()

            wsub = winpca[
                (winpca["scaffold"].astype(str) == chrom) &
                (winpca["pos"] >= start) &
                (winpca["pos"] <= end)
            ].copy()

            row = reg.to_dict()
            row["pair"] = pair
            row["length_bp"] = end - start + 1

            for label, col in [
                ("Fst", fst_col),
                ("dxy", dxy_col),
                ("dA", da_col),
                (f"pi_{pop1}", pi1_col if pi1_col in df.columns else None),
                (f"pi_{pop2}", pi2_col if pi2_col in df.columns else None),
            ]:
                if col is not None and col in sub.columns:
                    vals = pd.to_numeric(sub[col], errors="coerce")
                    row[f"mean_{label}"] = vals.mean()
                    row[f"max_{label}"] = vals.max()
                else:
                    row[f"mean_{label}"] = np.nan
                    row[f"max_{label}"] = np.nan

            if len(wsub) > 0:
                row["winpca_windows"] = len(wsub)
                row["mean_abs_delta_pc1"] = wsub["abs_delta_pc1"].mean()
                row["max_abs_delta_pc1"] = wsub["abs_delta_pc1"].max()
                row["mean_pc1_var_explained"] = wsub["pc_1_ve"].mean()
                row["max_pc1_var_explained"] = wsub["pc_1_ve"].max()
                row["mean_n_variants_winpca"] = wsub["n_var"].mean()

                row["winpca_score"] = (
                    row["mean_abs_delta_pc1"]
                    * (row["mean_pc1_var_explained"] / 100)
                    * np.log1p(row["winpca_windows"])
                )

                if row["mean_abs_delta_pc1"] >= 30 and row["winpca_windows"] >= 5:
                    row["winpca_class"] = "strong"
                elif row["mean_abs_delta_pc1"] >= 15 and row["winpca_windows"] >= 3:
                    row["winpca_class"] = "moderate"
                elif row["mean_abs_delta_pc1"] >= 8:
                    row["winpca_class"] = "weak"
                else:
                    row["winpca_class"] = "little_or_none"
            else:
                row["winpca_windows"] = 0
                row["mean_abs_delta_pc1"] = np.nan
                row["max_abs_delta_pc1"] = np.nan
                row["mean_pc1_var_explained"] = np.nan
                row["max_pc1_var_explained"] = np.nan
                row["mean_n_variants_winpca"] = np.nan
                row["winpca_score"] = np.nan
                row["winpca_class"] = "no_overlap"

            rows.append(row)

        out = pd.DataFrame(rows)

        if not out.empty:
            out = out.sort_values(
                ["winpca_score", "max_primary"],
                ascending=[False, False],
                na_position="last",
            )
            out.insert(0, "rank", range(1, len(out) + 1))

        out.to_csv(output.tsv, sep="\t", index=False)

        html_table = (
            out.to_html(index=False, classes="candidates")
            if not out.empty
            else "<p>No candidate regions found.</p>"
        )

        html = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{pair} Manhattan × WinPCA candidate regions</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    table {{
      border-collapse: collapse;
      font-size: 12px;
      background: white;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 4px 8px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>{pair}: Manhattan × WinPCA candidate regions</h1>

  <p>
    Candidate regions are merged {params.primary_metric} outlier windows
    above the {params.outlier_q} quantile. WinPCA support is estimated from
    local PC1 separation between {pop1} and {pop2}.
  </p>

  <p>
    <b>Download:</b>
    <a href="{os.path.basename(output.tsv)}">{os.path.basename(output.tsv)}</a>
  </p>

  {html_table}
</body>
</html>
"""

        with open(output.html, "w") as out_html:
            out_html.write(html)

rule summarize_manhattan_reports:
    input:
        htmls = lambda wc: [
            f"{MAN_ROOT}/{pair}/{pair}.summary.html"
            for pair in manhattan_pair_names()
        ]
    output:
        html = f"{MAN_ROOT}/manhattan_index.html"
    run:
        import os

        os.makedirs(MAN_ROOT, exist_ok=True)

        candidate_headers = ""
        if MAN_COMPARE_WINPCA:
            candidate_headers = """
              <th>Candidate regions</th>
              <th>Candidate TSV</th>
            """

        links = []
        for html in input.htmls:
            pair = os.path.basename(os.path.dirname(html))
            rel = os.path.relpath(html, MAN_ROOT)
            pdf = f"{pair}/{pair}.manhattan.pdf"
            png = f"{pair}/{pair}.manhattan.png"
            outliers = f"{pair}/{pair}.outlier_windows.tsv"

            candidate_cells = ""
            if MAN_COMPARE_WINPCA:
                candidates = f"{pair}/{pair}.candidate_regions.html"
                candidate_tsv = f"{pair}/{pair}.candidate_regions.tsv"
                candidate_cells = f"""
                  <td><a href="{candidates}">Candidate regions</a></td>
                  <td><a href="{candidate_tsv}">Candidate TSV</a></td>
                """

            links.append(f"""
            <tr>
              <td>{pair}</td>
              <td><a href="{rel}">HTML summary</a></td>
              <td><a href="{pdf}">PDF plot</a></td>
              <td><a href="{png}">PNG plot</a></td>
              <td><a href="{outliers}">Outlier TSV</a></td>
              {candidate_cells}
            </tr>
            """)

        page = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>SnakePop Manhattan reports</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    table {{
      border-collapse: collapse;
      background: white;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 8px 12px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>SnakePop Manhattan reports</h1>
  <table>
    <tr>
      <th>Pair</th>
      <th>HTML</th>
      <th>PDF</th>
      <th>PNG</th>
      <th>Outliers</th>
      {candidate_headers}
    </tr>
    {''.join(links)}
  </table>
</body>
</html>
"""

        with open(output.html, "w") as out:
            out.write(page)


rule manhattan:
    input:
        manhattan_targets()


rule popgenwindows:
    input:
        merged = PGW_MERGED_TARGET,
        summary = PGW_MERGED_TARGET.replace(".csv.gz", ".summary.tsv"),
        manhattan = f"{MAN_ROOT}/manhattan_index.html"


###############################################################################
# Heterozygosity and inbreeding coefficient
###############################################################################

HET_CFG = config["popstats"].get("heterozygosity", {})
HET_ROOT = f"{POP_ROOT}/heterozygosity"
HET_COLOR_BY = HET_CFG.get("color_by", config["popstats"].get("population_column", "morphology"))

HET_VCF = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz"
HET_RAW = f"{HET_ROOT}/vcftools.het"
HET_TABLE = f"{HET_ROOT}/heterozygosity.per_sample.tsv"
HET_PDF = f"{HET_ROOT}/heterozygosity.per_sample.pdf"
HET_HTML = f"{HET_ROOT}/heterozygosity.report.html"


rule vcftools_heterozygosity:
    input:
        vcf = HET_VCF
    output:
        het = HET_RAW
    threads: 1
    params:
        out = f"{HET_ROOT}/vcftools"
    shell:
        r"""
        mkdir -p {HET_ROOT}

        vcftools \
          --gzvcf {input.vcf} \
          --het \
          --out {params.out}
        """


rule summarize_heterozygosity:
    input:
        het = HET_RAW,
        metadata = SAMPLE_TABLE
    output:
        table = HET_TABLE,
        pdf = HET_PDF,
        html = HET_HTML
    params:
        color_by = HET_COLOR_BY
    run:
        import os
        import pandas as pd
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        os.makedirs(HET_ROOT, exist_ok=True)

        het = pd.read_csv(input.het, sep=r"\s+")
        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        het = het.rename(columns={"INDV": SAMPLE_COL})

        het["observed_heterozygosity"] = (
            het["N_SITES"] - het["O(HOM)"]
        ) / het["N_SITES"]

        het["expected_heterozygosity"] = (
            het["N_SITES"] - het["E(HOM)"]
        ) / het["N_SITES"]

        het = het.rename(columns={
            "O(HOM)": "observed_homozygotes",
            "E(HOM)": "expected_homozygotes",
            "N_SITES": "n_sites",
            "F": "inbreeding_coefficient_F",
        })

        df = het.merge(meta, on=SAMPLE_COL, how="left")
        df.to_csv(output.table, sep="\t", index=False)

        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42
        mpl.rcParams["font.family"] = "DejaVu Sans"

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(8.27, 11.69),
            constrained_layout=True,
        )

        color_by = params.color_by

        def plot_by_group(ax, ycol, ylabel):
            if color_by in df.columns:
                groups = list(df[color_by].fillna("NA").unique())
                positions = range(1, len(groups) + 1)

                for pos, group in zip(positions, groups):
                    sub = df[df[color_by].fillna("NA") == group]
                    x = [pos] * len(sub)

                    ax.scatter(
                        x,
                        sub[ycol],
                        s=35,
                        alpha=0.75,
                        edgecolors="black",
                        linewidths=0.4,
                        rasterized=True,
                    )

                ax.set_xticks(list(positions))
                ax.set_xticklabels(groups, rotation=45, ha="right")
                ax.set_xlabel(color_by)
            else:
                ax.scatter(
                    range(len(df)),
                    df[ycol],
                    s=35,
                    alpha=0.75,
                    edgecolors="black",
                    linewidths=0.4,
                    rasterized=True,
                )
                ax.set_xlabel("Sample")

            ax.set_ylabel(ylabel)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        plot_by_group(
            axes[0],
            "observed_heterozygosity",
            "Observed heterozygosity",
        )

        plot_by_group(
            axes[1],
            "expected_heterozygosity",
            "Expected heterozygosity",
        )

        plot_by_group(
            axes[2],
            "inbreeding_coefficient_F",
            "Inbreeding coefficient F",
        )

        fig.suptitle("Genome-wide heterozygosity and inbreeding", fontsize=12)
        fig.savefig(output.pdf, dpi=300)
        plt.close(fig)

        html_table = df.to_html(index=False, classes="het")

        html = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>SnakePop heterozygosity report</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    table {{
      border-collapse: collapse;
      font-size: 12px;
      background: white;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 4px 8px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>SnakePop heterozygosity report</h1>

  <p><b>PDF:</b> <a href="{os.path.basename(output.pdf)}">{os.path.basename(output.pdf)}</a></p>
  <p><b>Table:</b> <a href="{os.path.basename(output.table)}">{os.path.basename(output.table)}</a></p>

  <h2>Per-sample heterozygosity and F</h2>
  {html_table}
</body>
</html>
"""

        with open(output.html, "w") as out:
            out.write(html)


rule heterozygosity:
    input:
        table = HET_TABLE,
        pdf = HET_PDF,
        html = HET_HTML




###############################################################################
# WinPCA
###############################################################################

WINPCA_CFG = config["popstats"].get("winpca", {})

WINPCA_ROOT = f"{POP_ROOT}/winpca"
WINPCA_REPORT_ROOT = f"reports/popstats/{CALLSET_ID}_{REF_NAME}/winpca"

WINPCA_BIN = WINPCA_CFG.get("executable", "bin/winpca/winpca")
WINPCA_WINDOW_SIZE = WINPCA_CFG.get("window_size", 1000000)
WINPCA_INCREMENT = WINPCA_CFG.get("increment", 50000)
WINPCA_MIN_MAF = WINPCA_CFG.get("min_maf", 0.01)
WINPCA_COLOR_BY = WINPCA_CFG.get(
    "color_by",
    config["popstats"].get("population_column", "morphology"),
)
WINPCA_PLOT_VAR = WINPCA_CFG.get("plot_var", 1)
WINPCA_PLOT_INTERVAL = WINPCA_CFG.get("plot_interval", 5)
WINPCA_PLOT_FORMAT = WINPCA_CFG.get("plot_format", "HTML,SVG")

WINPCA_FILTER = WINPCA_CFG.get("sample_filter", {})
WINPCA_FILTER_COLUMN = WINPCA_FILTER.get("column", None)
WINPCA_FILTER_VALUES = WINPCA_FILTER.get("values", [])

WINPCA_RES = RES.get("winpca", {})

WINPCA_SAMPLE_LIST = f"{POP_ROOT}/populations/winpca.samples.txt"


def chrom_region(wildcards):
    fai = config["ref"]["fasta"] + ".fai"

    with open(fai) as f:
        for line in f:
            chrom, length = line.split("\t")[:2]
            if chrom == wildcards.chrom:
                return f"{chrom}:1-{length}"

    raise ValueError(f"Chromosome not found in FASTA index: {wildcards.chrom}")


rule make_winpca_sample_list:
    input:
        metadata = SAMPLE_TABLE
    output:
        samples = WINPCA_SAMPLE_LIST
    run:
        import os
        import pandas as pd

        df = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in df.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if WINPCA_FILTER_COLUMN is not None:
            if WINPCA_FILTER_COLUMN not in df.columns:
                raise ValueError(
                    f"Missing WinPCA filter column in metadata: {WINPCA_FILTER_COLUMN}"
                )

            if not WINPCA_FILTER_VALUES:
                raise ValueError(
                    "WinPCA sample_filter.values is empty, but sample_filter.column is set."
                )

            df = df[df[WINPCA_FILTER_COLUMN].isin(WINPCA_FILTER_VALUES)]

        samples = df[SAMPLE_COL].dropna().drop_duplicates()

        if len(samples) < 3:
            raise ValueError(
                f"WinPCA sample list has only {len(samples)} samples. "
                "Check sample_filter settings."
            )

        os.makedirs(os.path.dirname(output.samples), exist_ok=True)

        samples.to_csv(
            output.samples,
            index=False,
            header=False,
        )


rule run_winpca_chrom:
    input:
        vcf = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi",
        samples = WINPCA_SAMPLE_LIST
    output:
        done = f"{WINPCA_ROOT}/{{chrom}}/{{chrom}}.winpca.done"
    threads: WINPCA_RES.get("threads", 4)
    resources:
        mem_mb = WINPCA_RES.get("mem_mb", 16000),
        walltime = WINPCA_RES.get("walltime", 4)
    params:
        prefix = lambda wc: f"{WINPCA_ROOT}/{wc.chrom}/{wc.chrom}",
        region = chrom_region,
        window_size = WINPCA_WINDOW_SIZE,
        increment = WINPCA_INCREMENT,
        min_maf = WINPCA_MIN_MAF
    shell:
        r"""
        mkdir -p {WINPCA_ROOT}/{wildcards.chrom}

        {WINPCA_BIN} pca \
          {params.prefix} \
          {input.vcf} \
          {params.region} \
          --threads {threads} \
          --samples {input.samples} \
          --window_size {params.window_size} \
          --increment {params.increment} \
          --min_maf {params.min_maf}

        touch {output.done}
        """


rule plot_winpca_chrom:
    input:
        done = f"{WINPCA_ROOT}/{{chrom}}/{{chrom}}.winpca.done",
        metadata = SAMPLE_TABLE
    output:
        html = f"{WINPCA_ROOT}/{{chrom}}/{{chrom}}.pc_1.html",
        done = f"{WINPCA_ROOT}/{{chrom}}/{{chrom}}.chromplot.done"
    threads: WINPCA_RES.get("threads", 4)
    resources:
        mem_mb = WINPCA_RES.get("mem_mb", 16000),
        walltime = WINPCA_RES.get("walltime", 4)
    params:
        prefix = lambda wc: f"{WINPCA_ROOT}/{wc.chrom}/{wc.chrom}",
        region = chrom_region,
        color_by = WINPCA_COLOR_BY,
        plot_var = WINPCA_PLOT_VAR,
        interval = WINPCA_PLOT_INTERVAL,
        fmt = "HTML"
    shell:
        r"""
        {WINPCA_BIN} chromplot \
          {params.prefix} \
          {params.region} \
          --threads {threads} \
          --metadata {input.metadata} \
          --groups {params.color_by} \
          --plot_var {params.plot_var} \
          --interval {params.interval} \
          --format {params.fmt}

        touch {output.done}
        """

rule merge_winpca_html:
    input:
        htmls = expand(
            f"{WINPCA_ROOT}/{{chrom}}/{{chrom}}.pc_1.html",
            chrom=CHROMOSOMES
        )
    output:
        html = f"{WINPCA_ROOT}/winpca_merged.html"
    run:
        import os

        rows = []

        for html_file in input.htmls:
            chrom = os.path.basename(os.path.dirname(html_file))
            rel = os.path.relpath(html_file, WINPCA_ROOT)

            rows.append(f"""
<section class="chrom-section">
  <h2>{chrom}</h2>
  <iframe
    src="{rel}"
    width="100%"
    height="720"
    loading="lazy"
    style="border:1px solid #ccc; background:white;"
  ></iframe>
</section>
""")

        page = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>SnakePop WinPCA report</title>
</head>
<body>
  <h1>SnakePop WinPCA report</h1>
  {''.join(rows)}
</body>
</html>
"""

        with open(output.html, "w", encoding="utf-8") as out:
            out.write(page)

rule winpca:
    input:
        f"{WINPCA_ROOT}/winpca_merged.html"


###############################################################################
# Runs of homozygosity
###############################################################################

ROH_CFG = config["popstats"].get("roh", {})
ROH_ROOT = f"{POP_ROOT}/roh"

ROH_COLOR_BY = ROH_CFG.get(
    "color_by",
    config["popstats"].get("population_column", "morphology"),
)

ROH_MIN_KB = ROH_CFG.get("min_kb", 250)
ROH_MIN_SNP = ROH_CFG.get("min_snp", 100)
ROH_DENSITY = ROH_CFG.get("density_kb_per_snp", 50)
ROH_MAX_GAP = ROH_CFG.get("max_gap_kb", 500)

ROH_WINDOW_SNP = ROH_CFG.get("window_snp", 50)
ROH_WINDOW_HET = ROH_CFG.get("window_het", 1)
ROH_WINDOW_MISSING = ROH_CFG.get("window_missing", 5)
ROH_WINDOW_THRESHOLD = ROH_CFG.get("window_threshold", 0.05)

ROH_LENGTH_CLASSES = ROH_CFG.get("length_classes", {
    "short": [250, 500],
    "medium": [500, 1000],
    "long": [1000, 5000],
    "very_long": [5000, None],
})

ROH_RES = RES.get("roh", {})

ROH_RAW = f"{ROH_ROOT}/plink_roh.hom"
ROH_INDIV = f"{ROH_ROOT}/plink_roh.hom.indiv"
ROH_TABLE = f"{ROH_ROOT}/roh.per_sample.tsv"
ROH_SEGMENTS = f"{ROH_ROOT}/roh.segments.tsv"
ROH_PDF = f"{ROH_ROOT}/roh.per_sample.pdf"
ROH_HTML = f"{ROH_ROOT}/roh.report.html"


rule plink_roh:
    input:
        vcf = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz"
    output:
        hom = ROH_RAW,
        indiv = ROH_INDIV
    threads: ROH_RES.get("threads", 1)
    resources:
        mem_mb = ROH_RES.get("mem_mb", 8000),
        walltime = ROH_RES.get("walltime", 4)
    params:
        out = f"{ROH_ROOT}/plink_roh",
        chr_set = PLINK_CHR_SET,
        min_kb = ROH_MIN_KB,
        min_snp = ROH_MIN_SNP,
        density = ROH_DENSITY,
        max_gap = ROH_MAX_GAP,
        window_snp = ROH_WINDOW_SNP,
        window_het = ROH_WINDOW_HET,
        window_missing = ROH_WINDOW_MISSING,
        window_threshold = ROH_WINDOW_THRESHOLD
    shell:
        r"""
        mkdir -p {ROH_ROOT}

        plink \
          --vcf {input.vcf} \
          --double-id \
          --allow-extra-chr \
          --chr-set {params.chr_set} \
          --homozyg \
          --homozyg-kb {params.min_kb} \
          --homozyg-snp {params.min_snp} \
          --homozyg-density {params.density} \
          --homozyg-gap {params.max_gap} \
          --homozyg-window-snp {params.window_snp} \
          --homozyg-window-het {params.window_het} \
          --homozyg-window-missing {params.window_missing} \
          --homozyg-window-threshold {params.window_threshold} \
          --out {params.out}
        """


rule summarize_roh:
    input:
        hom = ROH_RAW,
        indiv = ROH_INDIV,
        metadata = SAMPLE_TABLE
    output:
        table = ROH_TABLE,
        segments = ROH_SEGMENTS,
        pdf = ROH_PDF,
        html = ROH_HTML
    params:
        color_by = ROH_COLOR_BY
    run:
        import os
        import pandas as pd
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        os.makedirs(ROH_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        indiv = pd.read_csv(input.indiv, sep=r"\s+")
        indiv = indiv.rename(columns={
            "IID": SAMPLE_COL,
            "NSEG": "n_roh",
            "KB": "total_roh_kb",
            "KBAVG": "mean_roh_kb",
        })

        keep_cols = [
            SAMPLE_COL,
            "n_roh",
            "total_roh_kb",
            "mean_roh_kb",
        ]

        indiv = indiv[keep_cols].copy()

        if os.path.getsize(input.hom) > 0:
            seg = pd.read_csv(input.hom, sep=r"\s+")
            seg = seg.rename(columns={
                "IID": SAMPLE_COL,
                "CHR": "scaffold",
                "POS1": "start",
                "POS2": "end",
                "KB": "length_kb",
                "NSNP": "n_snp",
            })
        else:
            seg = pd.DataFrame(columns=[
                SAMPLE_COL,
                "scaffold",
                "start",
                "end",
                "length_kb",
                "n_snp",
            ])

        df = indiv.merge(meta, on=SAMPLE_COL, how="left")

        if not seg.empty:
            max_roh = (
                seg.groupby(SAMPLE_COL)["length_kb"]
                .max()
                .reset_index()
                .rename(columns={"length_kb": "max_roh_kb"})
            )
            df = df.merge(max_roh, on=SAMPLE_COL, how="left")
        else:
            df["max_roh_kb"] = 0

        df["max_roh_kb"] = df["max_roh_kb"].fillna(0)

        df["total_roh_mb"] = df["total_roh_kb"] / 1000
        df["mean_roh_mb"] = df["mean_roh_kb"] / 1000
        df["max_roh_mb"] = df["max_roh_kb"] / 1000

        genome_size_bp = 0
        with open(CHROM_LIST) as f:
            for line in f:
                chrom = line.strip()
                if not chrom:
                    continue

        fai = config["ref"]["fasta"] + ".fai"
        chrom_set = set(CHROMOSOMES)

        with open(fai) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts[0] in chrom_set:
                    genome_size_bp += int(parts[1])

        if genome_size_bp > 0:
            df["FROH"] = (df["total_roh_kb"] * 1000) / genome_size_bp
        else:
            df["FROH"] = pd.NA

        for class_name, bounds in ROH_LENGTH_CLASSES.items():
            low, high = bounds
            low = float(low)

            if seg.empty:
                df[f"n_roh_{class_name}"] = 0
                df[f"total_roh_mb_{class_name}"] = 0.0
                continue

            if high is None:
                class_seg = seg[seg["length_kb"] >= low].copy()
            else:
                high = float(high)
                class_seg = seg[
                    (seg["length_kb"] >= low) &
                    (seg["length_kb"] < high)
                ].copy()

            if class_seg.empty:
                df[f"n_roh_{class_name}"] = 0
                df[f"total_roh_mb_{class_name}"] = 0.0
                continue

            tmp = (
                class_seg.groupby(SAMPLE_COL)
                .agg(
                    **{
                        f"n_roh_{class_name}": ("length_kb", "count"),
                        f"total_roh_mb_{class_name}": (
                            "length_kb",
                            lambda x: x.sum() / 1000,
                        ),
                    }
                )
                .reset_index()
            )

            df = df.merge(tmp, on=SAMPLE_COL, how="left")
            df[f"n_roh_{class_name}"] = (
                df[f"n_roh_{class_name}"]
                .fillna(0)
                .astype(int)
            )
            df[f"total_roh_mb_{class_name}"] = (
                df[f"total_roh_mb_{class_name}"]
                .fillna(0.0)
            )

        df.to_csv(output.table, sep="\t", index=False)
        seg.to_csv(output.segments, sep="\t", index=False)

        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42
        mpl.rcParams["font.family"] = "DejaVu Sans"

        fig, axes = plt.subplots(
            4,
            1,
            figsize=(8.27, 11.69),
            constrained_layout=True,
        )

        color_by = params.color_by

        def plot_grouped(ax, ycol, ylabel):
            if color_by in df.columns:
                groups = list(df[color_by].fillna("NA").unique())

                for pos, group in enumerate(groups, start=1):
                    sub = df[df[color_by].fillna("NA") == group]

                    ax.scatter(
                        [pos] * len(sub),
                        sub[ycol],
                        s=35,
                        alpha=0.75,
                        edgecolors="black",
                        linewidths=0.4,
                        rasterized=True,
                    )

                ax.set_xticks(range(1, len(groups) + 1))
                ax.set_xticklabels(groups, rotation=45, ha="right")
                ax.set_xlabel(color_by)
            else:
                ax.scatter(
                    range(len(df)),
                    df[ycol],
                    s=35,
                    alpha=0.75,
                    edgecolors="black",
                    linewidths=0.4,
                    rasterized=True,
                )
                ax.set_xlabel("Sample")

            ax.set_ylabel(ylabel)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        plot_grouped(axes[0], "n_roh", "Number of ROH")
        plot_grouped(axes[1], "total_roh_mb", "Total ROH [Mb]")
        plot_grouped(axes[2], "max_roh_mb", "Longest ROH [Mb]")
        plot_grouped(axes[3], "FROH", "FROH")

        fig.suptitle("Runs of homozygosity", fontsize=12)
        fig.savefig(output.pdf, dpi=300)
        plt.close(fig)

        html_table = df.to_html(index=False, classes="roh")

        html = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>SnakePop ROH report</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    table {{
      border-collapse: collapse;
      font-size: 12px;
      background: white;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 4px 8px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>SnakePop ROH report</h1>

  <p><b>PDF:</b> <a href="{os.path.basename(output.pdf)}">{os.path.basename(output.pdf)}</a></p>
  <p><b>Per-sample table:</b> <a href="{os.path.basename(output.table)}">{os.path.basename(output.table)}</a></p>
  <p><b>ROH segments:</b> <a href="{os.path.basename(output.segments)}">{os.path.basename(output.segments)}</a></p>

  <h2>Per-sample ROH summary</h2>
  {html_table}
</body>
</html>
"""

        with open(output.html, "w") as out:
            out.write(html)


rule roh:
    input:
        table = ROH_TABLE,
        segments = ROH_SEGMENTS,
        pdf = ROH_PDF,
        html = ROH_HTML



###############################################################################
# SNP window trees with IQ-TREE
###############################################################################

SNPTREE_CFG = config["popstats"].get("snptrees_iqtree", {})

SNPTREE_ROOT = f"{POP_ROOT}/snptrees_iqtree"

SNPTREE_SAMPLE_LIST = f"{SNPTREE_ROOT}/samples.txt"
SNPTREE_WINDOWS = f"{SNPTREE_ROOT}/windows.tsv"

SNPTREE_WINDOW_SIZE = SNPTREE_CFG.get("window_size", 1000000)
SNPTREE_INCREMENT = SNPTREE_CFG.get("increment", 1000000)
SNPTREE_MIN_SNPS = SNPTREE_CFG.get("min_snps", 500)

SNPTREE_VCF2PHYLIP = SNPTREE_CFG.get(
    "vcf2phylip",
    "bin/vcf2phylip/vcf2phylip.py",
)
SNPTREE_IQTREE = SNPTREE_CFG.get("iqtree", "iqtree")
SNPTREE_MODEL = SNPTREE_CFG.get("model", "GTR")
SNPTREE_UFBOOTS = SNPTREE_CFG.get("UFBoots", 0)
SNPTREE_FAST = SNPTREE_CFG.get("fast", True)
SNPTREE_CLEANUP = SNPTREE_CFG.get("cleanup", True)
SNPTREE_TMPDIR = SNPTREE_CFG.get("tmpdir", None)

SNPTREE_FILTER = SNPTREE_CFG.get("sample_filter", {})
SNPTREE_FILTER_IDS = SNPTREE_FILTER.get("ids", [])
SNPTREE_FILTER_COLUMN = SNPTREE_FILTER.get("column", None)
SNPTREE_FILTER_VALUES = SNPTREE_FILTER.get("values", [])

SNPTREE_RES = RES.get("snptrees_iqtree", {})

SNPTREE_CHROM_TREES = expand(
    f"{SNPTREE_ROOT}/{{chrom}}.window_trees.tsv.gz",
    chrom=CHROMOSOMES,
)

SNPTREE_ALL_TREES = f"{SNPTREE_ROOT}/all.window_trees.tsv.gz"
SNPTREE_SUMMARY = f"{SNPTREE_ROOT}/snptrees_iqtree.summary.tsv"


rule make_snptree_sample_list:
    input:
        metadata = SAMPLE_TABLE
    output:
        samples = SNPTREE_SAMPLE_LIST
    run:
        import os
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if SNPTREE_FILTER_IDS:
            samples = pd.Series(SNPTREE_FILTER_IDS).dropna().drop_duplicates()

        elif SNPTREE_FILTER_COLUMN is not None:
            if SNPTREE_FILTER_COLUMN not in meta.columns:
                raise ValueError(
                    f"Missing snptrees_iqtree sample_filter column: {SNPTREE_FILTER_COLUMN}"
                )

            if not SNPTREE_FILTER_VALUES:
                raise ValueError(
                    "snptrees_iqtree sample_filter.column is set, but values is empty."
                )

            samples = (
                meta.loc[
                    meta[SNPTREE_FILTER_COLUMN].isin(SNPTREE_FILTER_VALUES),
                    SAMPLE_COL,
                ]
                .dropna()
                .drop_duplicates()
            )

        else:
            samples = meta[SAMPLE_COL].dropna().drop_duplicates()

        if len(samples) < 4:
            raise ValueError(
                f"SNP tree sample list has only {len(samples)} samples. "
                "Use at least 4 samples for tree inference."
            )

        samples.to_csv(output.samples, index=False, header=False)


rule make_snptree_windows:
    input:
        fai = config["ref"]["fasta"] + ".fai"
    output:
        windows = SNPTREE_WINDOWS
    run:
        import os
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        chrom_set = set(CHROMOSOMES)
        rows = []

        with open(input.fai) as f:
            for line in f:
                chrom, length = line.split("\t")[:2]
                length = int(length)

                if chrom not in chrom_set:
                    continue

                start = 1
                while start <= length:
                    end = min(start + SNPTREE_WINDOW_SIZE - 1, length)

                    rows.append({
                        "tree_id": f"{chrom}:{start}-{end}",
                        "scaffold": chrom,
                        "start": start,
                        "end": end,
                    })

                    start += SNPTREE_INCREMENT

        pd.DataFrame(rows).to_csv(output.windows, sep="\t", index=False)


rule snptrees_iqtree_chrom:
    input:
        windows = SNPTREE_WINDOWS,
        samples = SNPTREE_SAMPLE_LIST,
        vcf = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi"
    output:
        trees = f"{SNPTREE_ROOT}/{{chrom}}.window_trees.tsv.gz"
    threads: SNPTREE_RES.get("threads", 4)
    resources:
        mem_mb = SNPTREE_RES.get("mem_mb", 16000),
        walltime = SNPTREE_RES.get("walltime", 24)
    params:
        model = SNPTREE_MODEL,
        min_snps = SNPTREE_MIN_SNPS,
        ufboots = SNPTREE_UFBOOTS,
        fast = SNPTREE_FAST,
        cleanup = SNPTREE_CLEANUP,
        tmpdir = SNPTREE_TMPDIR
    run:
        import os
        import gzip
        import glob
        import shutil
        import tempfile
        import subprocess
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        chrom = wildcards.chrom

        tmp_base = params.tmpdir
        if tmp_base is not None:
            os.makedirs(tmp_base, exist_ok=True)

        workdir = tempfile.mkdtemp(
            prefix=f"snakepop_snptrees_{chrom}_",
            dir=tmp_base,
        )

        windows = pd.read_csv(input.windows, sep="\t")
        windows = windows[windows["scaffold"].astype(str) == str(chrom)].copy()

        tmp_vcf = os.path.join(workdir, "window.vcf.gz")
        phy_prefix = "window"

        tmp_out = str(output.trees) + ".tmp"

        def clean_window_files():
            for fn in glob.glob(os.path.join(workdir, "window*")):
                try:
                    os.remove(fn)
                except FileNotFoundError:
                    pass

        def count_vcf_snps(vcf_file):
            cmd = ["bcftools", "view", "-H", vcf_file]
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
            p2 = subprocess.Popen(
                ["wc", "-l"],
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                text=True,
            )
            p1.stdout.close()
            out, _ = p2.communicate()
            return int(out.strip())

        def phylip_n_sites(phy_file):
            with open(phy_file) as f:
                first = f.readline().strip().split()

            if len(first) < 2:
                return 0

            try:
                return int(first[1])
            except ValueError:
                return 0

        try:
            with gzip.open(tmp_out, "wt") as out:
                out.write(
                    "tree_id\tscaffold\tstart\tend\tn_vcf_snps\tn_phy_sites\ttree\n"
                )

                for _, row in windows.iterrows():
                    tree_id = row["tree_id"]
                    start = int(row["start"])
                    end = int(row["end"])
                    region = f"{chrom}:{start}-{end}"

                    clean_window_files()

                    subprocess.run(
                        [
                            "bcftools", "view",
                            "-S", input.samples,
                            "-r", region,
                            "-m2", "-M2",
                            "-v", "snps",
                            "-Oz",
                            "-o", tmp_vcf,
                            input.vcf,
                        ],
                        check=True,
                    )

                    subprocess.run(
                        ["bcftools", "index", "-t", "-f", tmp_vcf],
                        check=True,
                    )

                    n_vcf_snps = count_vcf_snps(tmp_vcf)

                    if n_vcf_snps == 0:
                        clean_window_files()
                        continue

                    subprocess.run(
                        [
                            "python",
                            SNPTREE_VCF2PHYLIP,
                            "-i",
                            tmp_vcf,
                            "--output-folder",
                            workdir,
                            "--output-prefix",
                            phy_prefix,
                        ],
                        check=True,
                    )

                    phy_files = glob.glob(
                        os.path.join(workdir, f"{phy_prefix}*.phy")
                    )

                    if not phy_files:
                        clean_window_files()
                        continue

                    phy = phy_files[0]
                    n_phy_sites = phylip_n_sites(phy)

                    if n_phy_sites < params.min_snps:
                        clean_window_files()
                        continue

                    iq_prefix = os.path.join(workdir, "window.iqtree")

                    iq_cmd = [
                        SNPTREE_IQTREE,
                        "-s", phy,
                        "--seqtype", "DNA",
                        "-m", params.model,
                        "-T", "1",
                        "--prefix", iq_prefix,
                        "-quiet",
                        "-redo",
                    ]

                    if params.fast:
                        iq_cmd.append("--fast")

                    if int(params.ufboots) > 0:
                        iq_cmd += ["-B", str(params.ufboots)]

                    subprocess.run(iq_cmd, check=True)

                    treefile = f"{iq_prefix}.treefile"

                    if not os.path.exists(treefile):
                        raise ValueError(
                            f"IQ-TREE produced no treefile for {tree_id}"
                        )

                    with open(treefile) as tf:
                        tree = tf.read().strip()

                    out.write(
                        f"{tree_id}\t{chrom}\t{start}\t{end}\t"
                        f"{n_vcf_snps}\t{n_phy_sites}\t{tree}\n"
                    )

                    out.flush()
                    clean_window_files()

            os.replace(tmp_out, output.trees)

        finally:
            if params.cleanup:
                shutil.rmtree(workdir, ignore_errors=True)
                try:
                    os.remove(tmp_out)
                except FileNotFoundError:
                    pass
            else:
                print(f"Temporary SNP-tree files kept in: {workdir}")
                print(f"Temporary chromosome output kept as: {tmp_out}")


rule merge_snptrees_iqtree:
    input:
        trees = SNPTREE_CHROM_TREES
    output:
        all_trees = SNPTREE_ALL_TREES,
        summary = SNPTREE_SUMMARY
    run:
        import os
        import gzip
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        summaries = []

        with gzip.open(output.all_trees, "wt") as out:
            out.write(
                "tree_id\tscaffold\tstart\tend\tn_vcf_snps\tn_phy_sites\ttree\n"
            )

            for fn in input.trees:
                df = pd.read_csv(fn, sep="\t", compression="gzip")

                summaries.append({
                    "file": fn,
                    "n_trees": len(df),
                    "mean_phy_sites": df["n_phy_sites"].mean() if len(df) else 0,
                    "min_phy_sites": df["n_phy_sites"].min() if len(df) else 0,
                    "max_phy_sites": df["n_phy_sites"].max() if len(df) else 0,
                    "mean_vcf_snps": df["n_vcf_snps"].mean() if len(df) else 0,
                    "min_vcf_snps": df["n_vcf_snps"].min() if len(df) else 0,
                    "max_vcf_snps": df["n_vcf_snps"].max() if len(df) else 0,
                })

                for _, row in df.iterrows():
                    out.write(
                        f"{row['tree_id']}\t{row['scaffold']}\t{row['start']}\t"
                        f"{row['end']}\t{row['n_vcf_snps']}\t"
                        f"{row['n_phy_sites']}\t{row['tree']}\n"
                    )

        pd.DataFrame(summaries).to_csv(output.summary, sep="\t", index=False)


rule snptrees_iqtree:
    input:
        samples = SNPTREE_SAMPLE_LIST,
        windows = SNPTREE_WINDOWS,
        all_trees = SNPTREE_ALL_TREES,
        summary = SNPTREE_SUMMARY


###############################################################################
# ASTRAL summary tree from SNP window trees
###############################################################################

ASTRAL_CFG = config["popstats"].get("astral", {})
ASTRAL_ROOT = f"{POP_ROOT}/astral"

ASTRAL_EXE = ASTRAL_CFG.get("executable", "astral")
ASTRAL_MAPPING_COLUMN = ASTRAL_CFG.get(
    "mapping_column",
    config["popstats"].get("population_column", "morphology"),
)
ASTRAL_BRANCH_ANNOTATE = ASTRAL_CFG.get("branch_annotate", 3)

ASTRAL_ROOT_CFG = ASTRAL_CFG.get("root", {})
ASTRAL_DO_ROOT = ASTRAL_ROOT_CFG.get("enabled", True)
ASTRAL_OUTGROUPS = ASTRAL_ROOT_CFG.get("outgroups", [])

ASTRAL_RES = RES.get("astral", {})

ASTRAL_GENE_TREES = f"{ASTRAL_ROOT}/window_trees.newick"
ASTRAL_MAPPING = f"{ASTRAL_ROOT}/astral_mapping.tsv"
ASTRAL_TREE = f"{ASTRAL_ROOT}/astral.tree"
ASTRAL_TOPOLOGY = f"{ASTRAL_ROOT}/astral.topology.tree"
ASTRAL_LOG = f"{ASTRAL_ROOT}/astral.log"


rule prepare_astral_input:
    input:
        trees = SNPTREE_ALL_TREES,
        metadata = SAMPLE_TABLE
    output:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING
    run:
        import os
        import pandas as pd

        os.makedirs(ASTRAL_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if ASTRAL_MAPPING_COLUMN not in meta.columns:
            raise ValueError(
                f"Missing ASTRAL mapping column in metadata: {ASTRAL_MAPPING_COLUMN}"
            )

        meta = meta[[SAMPLE_COL, ASTRAL_MAPPING_COLUMN]].dropna().copy()
        meta[ASTRAL_MAPPING_COLUMN] = meta[ASTRAL_MAPPING_COLUMN].map(clean_pop_name)

        # ASTRAL mapping format:
        # species: sample1,sample2,sample3
        with open(output.mapping, "w") as out:
            for group, sub in meta.groupby(ASTRAL_MAPPING_COLUMN):
                group_samples = sorted(sub[SAMPLE_COL].dropna().unique())
                if group_samples:
                    out.write(f"{group}: {','.join(group_samples)}\n")

        df = pd.read_csv(input.trees, sep="\t", compression="gzip")

        if "tree" not in df.columns:
            raise ValueError("Missing tree column in SNP-tree table.")

        with open(output.gene_trees, "w") as out:
            for tree in df["tree"].dropna():
                tree = str(tree).strip()
                if tree:
                    out.write(tree + "\n")


rule astral_tree:
    input:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING
    output:
        tree = ASTRAL_TREE,
        log = ASTRAL_LOG
    threads: ASTRAL_RES.get("threads", 1)
    resources:
        mem_mb = ASTRAL_RES.get("mem_mb", 16000),
        walltime = ASTRAL_RES.get("walltime", 4)
    params:
        exe = ASTRAL_EXE,
        branch_annotate = ASTRAL_BRANCH_ANNOTATE,
        do_root = ASTRAL_DO_ROOT,
        outgroups = ASTRAL_OUTGROUPS
    run:
        import os
        import subprocess

        os.makedirs(ASTRAL_ROOT, exist_ok=True)

        cmd = [
            params.exe,
            "-i", input.gene_trees,
            "-a", input.mapping,
            "-o", output.tree,
            "-t", str(params.branch_annotate),
        ]

        if params.do_root and params.outgroups:
            cmd += ["--outgroup", ",".join(params.outgroups)]

        with open(output.log, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)



rule astral_topology_tree:
    input:
        tree = ASTRAL_TREE
    output:
        topology = ASTRAL_TOPOLOGY
    params:
        do_root = ASTRAL_DO_ROOT,
        outgroups = ASTRAL_OUTGROUPS
    run:
        import re
        from ete3 import Tree

        with open(input.tree) as f:
            newick = f.read().strip()

        # Remove ASTRAL bracket annotations, e.g. [q1=...;q2=...]
        newick = re.sub(r"\[[^\[\]]*\]", "", newick)

        # Parse cleaned Newick
        t = Tree(newick, format=1)

        if params.do_root and params.outgroups:
            outgroup = params.outgroups[0]

            try:
                t.set_outgroup(t & outgroup)
            except Exception:
                tips = sorted([leaf.name for leaf in t.iter_leaves()])
                raise ValueError(
                    f"Outgroup '{outgroup}' not found in ASTRAL tree. "
                    f"Available tips: {tips}"
                )

        # Write plain topology only: no branch lengths, no support values
        topology = t.write(format=9).strip()

        # Optional: normalize order so outgroup is visibly first
        if params.do_root and params.outgroups:
            outgroup = params.outgroups[0]
            if not topology.startswith(f"({outgroup},"):
                # topology is still correctly rooted, this only affects display order
                pass

        with open(output.topology, "w") as out:
            out.write(topology + "\n")

rule astral:
    input:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING,
        tree = ASTRAL_TREE,
        topology = ASTRAL_TOPOLOGY,
        log = ASTRAL_LOG



###############################################################################
# Dsuite Dtrios: per-chromosome runs + DtriosCombine
###############################################################################

DSUITE_CFG = config["popstats"].get("dsuite", {})
DSUITE_ROOT = f"{POP_ROOT}/dsuite"

DSUITE_EXE = DSUITE_CFG.get("executable", "Dsuite")
DSUITE_MAPPING_COLUMN = DSUITE_CFG.get(
    "mapping_column",
    config["popstats"].get("population_column", "morphology"),
)
DSUITE_OUTGROUPS = DSUITE_CFG.get("outgroups", [])
DSUITE_TREE = DSUITE_CFG.get("tree", ASTRAL_TOPOLOGY)

DSUITE_JKNUM = DSUITE_CFG.get("jknum", 20)
DSUITE_JKWINDOW = DSUITE_CFG.get("jkwindow", None)
DSUITE_NO_F4_RATIO = DSUITE_CFG.get("no_f4_ratio", False)
DSUITE_ABBA_CLUSTERING = DSUITE_CFG.get("abba_clustering", False)
DSUITE_USE_GENOTYPE_PROBABILITIES = DSUITE_CFG.get("use_genotype_probabilities", False)

DSUITE_RES = RES.get("dsuite", {})

DSUITE_SETS = f"{DSUITE_ROOT}/sets.txt"
DSUITE_COMBINED_PREFIX = f"{DSUITE_ROOT}/combined/dsuite"
DSUITE_COMBINED_LOG = f"{DSUITE_ROOT}/combined/dsuite_combine.log"

DSUITE_CHROM_DMIN = expand(
    f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_Dmin.txt",
    chrom=CHROMOSOMES,
)

DSUITE_COMBINED_BBAA = f"{DSUITE_COMBINED_PREFIX}_combined_BBAA.txt"
DSUITE_COMBINED_DMIN = f"{DSUITE_COMBINED_PREFIX}_combined_Dmin.txt"
DSUITE_COMBINED_TREE = f"{DSUITE_COMBINED_PREFIX}_combined_tree.txt"


rule prepare_dsuite_sets:
    input:
        metadata = SAMPLE_TABLE
    output:
        sets = DSUITE_SETS
    run:
        import os
        import pandas as pd

        os.makedirs(DSUITE_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if DSUITE_MAPPING_COLUMN not in meta.columns:
            raise ValueError(
                f"Missing Dsuite mapping column in metadata: {DSUITE_MAPPING_COLUMN}"
            )

        if not DSUITE_OUTGROUPS:
            raise ValueError(
                "No Dsuite outgroup specified. Add popstats: dsuite: outgroups:"
            )

        df = meta[[SAMPLE_COL, DSUITE_MAPPING_COLUMN]].dropna().copy()
        df[DSUITE_MAPPING_COLUMN] = df[DSUITE_MAPPING_COLUMN].map(clean_pop_name)

        outgroups = set(clean_pop_name(x) for x in DSUITE_OUTGROUPS)

        with open(output.sets, "w") as out:
            for _, row in df.iterrows():
                sample = row[SAMPLE_COL]
                group = row[DSUITE_MAPPING_COLUMN]

                if group in outgroups:
                    group = "Outgroup"

                out.write(f"{sample}\t{group}\n")


rule dsuite_dtrios_chrom:
    input:
        vcf = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi",
        sets = DSUITE_SETS,
        tree = DSUITE_TREE
    output:
        bbaa = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_BBAA.txt",
        dmin = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_Dmin.txt",
        tree = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_tree.txt",
        combine = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_combine.txt",
        combine_stderr = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}_combine_stderr.txt",
        log = f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}.log"
    threads: DSUITE_RES.get("threads", 1)
    resources:
        mem_mb = DSUITE_RES.get("mem_mb", 16000),
        walltime = DSUITE_RES.get("walltime", 24)
    params:
        exe = DSUITE_EXE,
        prefix = lambda wc: f"{DSUITE_ROOT}/{wc.chrom}/dsuite_{wc.chrom}",
        jknum = DSUITE_JKNUM,
        jkwindow = DSUITE_JKWINDOW,
        no_f4_ratio = DSUITE_NO_F4_RATIO,
        abba_clustering = DSUITE_ABBA_CLUSTERING,
        use_gp = DSUITE_USE_GENOTYPE_PROBABILITIES
    run:
        import os
        import subprocess

        os.makedirs(os.path.dirname(output.log), exist_ok=True)

        cmd = [
            params.exe,
            "Dtrios",
            "-t", input.tree,
            "-o", params.prefix,
        ]

        if params.jkwindow is not None:
            cmd += ["-j", str(params.jkwindow)]
        else:
            cmd += ["-k", str(params.jknum)]

        if params.no_f4_ratio:
            cmd.append("--no-f4-ratio")

        if params.abba_clustering:
            cmd.append("--ABBAclustering")

        if params.use_gp:
            cmd.append("-g")

        cmd += [
            input.vcf,
            input.sets,
        ]

        with open(output.log, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


rule combine_dsuite_dtrios:
    input:
        dmins = DSUITE_CHROM_DMIN,
        tree = DSUITE_TREE
    output:
        check = f"{DSUITE_ROOT}/combined/dsuite_combine.checkpoint",
        bbaa = DSUITE_COMBINED_BBAA,
        dmin = DSUITE_COMBINED_DMIN,
        tree = DSUITE_COMBINED_TREE,
        log = DSUITE_COMBINED_LOG
    params:
        exe = DSUITE_EXE,
        out_prefix = os.getcwd() + "/" + DSUITE_COMBINED_PREFIX,
        in_prefixes = expand(
            os.getcwd() + "/" + f"{DSUITE_ROOT}/{{chrom}}/dsuite_{{chrom}}",
            chrom = CHROMOSOMES,
        )
    run:
        import os
        import subprocess

        os.makedirs(os.path.dirname(output.log), exist_ok=True)

        cmd = [
            params.exe,
            "DtriosCombine",
            "-o", params.out_prefix,
            "-t", input.tree,
        ] + list(params.in_prefixes)

        with open(output.log, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)

        if not os.path.exists(output.tree) or os.path.getsize(output.tree) == 0:
            raise ValueError(f"Dsuite combine did not create non-empty {output.tree}")

        with open(output.check, "w") as f:
            f.write("OK\n")



rule dsuite:
    input:
        sets = DSUITE_SETS,
        chrom_dmins = DSUITE_CHROM_DMIN,
        combined_check = f"{DSUITE_ROOT}/combined/dsuite_combine.checkpoint",
        combined_bbaa = DSUITE_COMBINED_BBAA,
        combined_dmin = DSUITE_COMBINED_DMIN,
        combined_tree = DSUITE_COMBINED_TREE,
        log = DSUITE_COMBINED_LOG


###############################################################################
# Twisst topology weighting from SNP-window trees
###############################################################################

TWISST_CFG = config["popstats"].get("twisst", {})
TWISST_ROOT = f"{POP_ROOT}/twisst"

TWISST_EXE = TWISST_CFG.get("executable", "bin/twisst/twisst.py")
TWISST_MAPPING_COLUMN = TWISST_CFG.get(
    "mapping_column",
    config["popstats"].get("population_column", "morphology"),
)
TWISST_OUTGROUP = TWISST_CFG.get("outgroup", None)
TWISST_METHOD = TWISST_CFG.get("method", "complete")
TWISST_ITERATIONS = TWISST_CFG.get("iterations", 100)
TWISST_ABORT_CUTOFF = TWISST_CFG.get("abort_cutoff", 100)

TWISST_PLOT_CFG = TWISST_CFG.get("plotting", {})
TWISST_PLOT_ENABLED = TWISST_PLOT_CFG.get("enabled", True)
TWISST_PLOT_TOPOLOGIES = TWISST_PLOT_CFG.get("topologies", 3)  # integer or "all"
TWISST_PLOT_COMBINE_OTHER = TWISST_PLOT_CFG.get("combine_other", True)
TWISST_PLOT_SMOOTH = TWISST_PLOT_CFG.get("smooth", 1)
TWISST_PLOT_STACKED_AREA = TWISST_PLOT_CFG.get("stacked_area", True)
TWISST_PLOT_LINE = TWISST_PLOT_CFG.get("line_plot", True)
TWISST_PLOT_WIDTH = TWISST_PLOT_CFG.get("width", 14)
TWISST_PLOT_HEIGHT = TWISST_PLOT_CFG.get("height", 8)
TWISST_PLOT_DPI = TWISST_PLOT_CFG.get("dpi", 300)
TWISST_PLOT_FORMATS = TWISST_PLOT_CFG.get("formats", ["pdf", "svg"])
TWISST_PLOT_CHR_LABELS = TWISST_PLOT_CFG.get("chromosome_labels", True)
TWISST_PLOT_CHR_LINES = TWISST_PLOT_CFG.get("chromosome_lines", True)
TWISST_PLOT_SHOW_PERCENTAGES = TWISST_PLOT_CFG.get("show_percentages", True)
TWISST_PLOT_TOPOLOGY_PANEL = TWISST_PLOT_CFG.get("topology_panel", True)

TWISST_RES = RES.get("twisst", {})

TWISST_TREEFILE = f"{TWISST_ROOT}/window_trees.newick"
TWISST_GROUPS_FILE = f"{TWISST_ROOT}/groups.tsv"
TWISST_TREE_METADATA = f"{TWISST_ROOT}/window_tree_metadata.tsv"
TWISST_WEIGHTS = f"{TWISST_ROOT}/twisst.weights.tsv"
TWISST_DISTANCES = f"{TWISST_ROOT}/twisst.distances.tsv"
TWISST_TOPOLOGIES = f"{TWISST_ROOT}/twisst.topologies.tsv"
TWISST_LOG = f"{TWISST_ROOT}/twisst.log"

TWISST_PLOTS_ROOT = f"{TWISST_ROOT}/plots"
TWISST_PLOT_PREFIX = f"{TWISST_PLOTS_ROOT}/twisst_topology_weights"

TWISST_PLOT_OUTPUTS = [
    f"{TWISST_PLOT_PREFIX}.{fmt.lower()}"
    for fmt in TWISST_PLOT_FORMATS
] if TWISST_PLOT_ENABLED else []


rule prepare_twisst_input:
    input:
        trees = SNPTREE_ALL_TREES,
        metadata = SAMPLE_TABLE
    output:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE,
        tree_metadata = TWISST_TREE_METADATA
    run:
        import os
        import pandas as pd

        os.makedirs(TWISST_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if TWISST_MAPPING_COLUMN not in meta.columns:
            raise ValueError(
                f"Missing Twisst mapping column in metadata: {TWISST_MAPPING_COLUMN}"
            )

        mapping = meta[[SAMPLE_COL, TWISST_MAPPING_COLUMN]].dropna().copy()
        mapping[TWISST_MAPPING_COLUMN] = mapping[TWISST_MAPPING_COLUMN].map(clean_pop_name)

        if mapping.empty:
            raise ValueError(
                f"No samples found for Twisst mapping column '{TWISST_MAPPING_COLUMN}'."
            )

        mapping.to_csv(output.groups, sep="\t", index=False, header=False)

        df = pd.read_csv(input.trees, sep="\t", compression="gzip")

        if "tree" not in df.columns:
            raise ValueError("Missing tree column in SNP-tree table.")

        required = ["tree_id", "scaffold", "start", "end", "n_vcf_snps", "n_phy_sites"]
        existing = [x for x in required if x in df.columns]

        if not {"scaffold", "start", "end"}.issubset(set(existing)):
            raise ValueError(
                "Twisst plotting requires scaffold/start/end columns in SNP-tree metadata."
            )

        df[existing].to_csv(output.tree_metadata, sep="\t", index=False)

        with open(output.treefile, "w") as out:
            for tree in df["tree"].dropna():
                tree = str(tree).strip()
                if tree:
                    out.write(tree + "\n")


rule run_twisst:
    input:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE
    output:
        weights = TWISST_WEIGHTS,
        distances = TWISST_DISTANCES,
        topologies = TWISST_TOPOLOGIES,
        log = TWISST_LOG
    threads: TWISST_RES.get("threads", 1)
    resources:
        mem_mb = TWISST_RES.get("mem_mb", 8000),
        walltime = TWISST_RES.get("walltime", 8)
    params:
        exe = TWISST_EXE,
        outgroup = TWISST_OUTGROUP,
        method = TWISST_METHOD,
        iterations = TWISST_ITERATIONS,
        abort_cutoff = TWISST_ABORT_CUTOFF
    run:
        import os
        import subprocess
        import pandas as pd

        os.makedirs(TWISST_ROOT, exist_ok=True)

        groups = pd.read_csv(
            input.groups,
            sep="\t",
            header=None,
            names=["sample", "group"],
            dtype=str,
        )

        group_args = []
        for group, sub in groups.groupby("group"):
            samples = sorted(sub["sample"].dropna().unique())
            if samples:
                group_args += ["-g", group] + list(samples)

        cmd = [
            "python",
            params.exe,
            "-t", input.treefile,
            "-w", output.weights,
            "-D", output.distances,
            "--outputTopos", output.topologies,
            "--method", params.method,
            "--abortCutoff", str(params.abort_cutoff),
            "--silent",
        ] + group_args

        if params.method == "fixed":
            cmd += ["--iterations", str(params.iterations)]

        if params.outgroup:
            cmd += ["--outgroup", clean_pop_name(params.outgroup)]

        with open(output.log, "w") as log:
            log.write("COMMAND:\n")
            log.write(" ".join(cmd) + "\n\n")
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


rule plot_twisst_topologies:
    input:
        weights = TWISST_WEIGHTS,
        metadata = TWISST_TREE_METADATA,
        topologies = TWISST_TOPOLOGIES
    output:
        plots = TWISST_PLOT_OUTPUTS
    params:
        script = "bin/plot_twisst_snakepop.R",
        out_prefix = TWISST_PLOT_PREFIX,
        topologies = TWISST_PLOT_TOPOLOGIES,
        combine_other = TWISST_PLOT_COMBINE_OTHER,
        smooth = TWISST_PLOT_SMOOTH,
        width = TWISST_PLOT_WIDTH,
        height = TWISST_PLOT_HEIGHT,
        formats = ",".join(TWISST_PLOT_FORMATS)
    shell:
        r"""
        Rscript {params.script} \
          --weights {input.weights} \
          --metadata {input.metadata} \
          --topologies {input.topologies} \
          --out-prefix {params.out_prefix} \
          --top-n {params.topologies} \
          --combine-other {params.combine_other} \
          --smooth {params.smooth} \
          --width {params.width} \
          --height {params.height} \
          --formats {params.formats}
        """

rule twisst:
    input:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE,
        metadata = TWISST_TREE_METADATA,
        weights = TWISST_WEIGHTS,
        distances = TWISST_DISTANCES,
        topologies = TWISST_TOPOLOGIES,
        log = TWISST_LOG,
        plots = TWISST_PLOT_OUTPUTS

############################################################

rule popstats:
    input:
        rules.pca.input,
        rules.popgenwindows.input,
        rules.winpca.input,
        rules.manhattan.input,
        rules.heterozygosity.input,
        rules.roh.input,
        rules.snptrees_iqtree.input,
        rules.astral.input,
        rules.dsuite.input,
        rules.twisst.input




