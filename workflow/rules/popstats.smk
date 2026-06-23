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
            _, pop1, pop2 = col.split("_", 2)

            pi1 = f"pi_{pop1}"
            pi2 = f"pi_{pop2}"

            if pi1 in merged.columns and pi2 in merged.columns:
                merged[f"dA_{pop1}_{pop2}"] = (
                    merged[col] - (merged[pi1] + merged[pi2]) / 2
                )

        os.makedirs(PGW_ROOT, exist_ok=True)

        merged.to_csv(
            output.merged,
            sep="\t",
            compression="gzip",
            index=False,
        )


rule popgenwindows:
    input:
        PGW_MERGED_TARGET

rule popstats:
    input:
        rules.pca.input,
        rules.popgenwindows.input
