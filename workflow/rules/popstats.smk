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

rule popgenwindows:
    input:
        merged = PGW_MERGED_TARGET,
        summary = PGW_MERGED_TARGET.replace(".csv.gz", ".summary.tsv")


rule popstats:
    input:
        rules.pca.input,
        rules.popgenwindows.input


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
        done = f"{WINPCA_REPORT_ROOT}/{{chrom}}/{{chrom}}.chromplot.done"
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
        fmt = WINPCA_PLOT_FORMAT
    shell:
        r"""
        mkdir -p {WINPCA_REPORT_ROOT}/{wildcards.chrom}

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
        import html as html_lib

        os.makedirs(WINPCA_ROOT, exist_ok=True)

        sections = []

        for html_file in input.htmls:
            chrom = os.path.basename(os.path.dirname(html_file))

            with open(html_file, "r", encoding="utf-8") as f:
                content = f.read()

            escaped = html_lib.escape(content, quote=True)

            sections.append(f"""
<section class="chrom-section">
  <h2>{chrom}</h2>
  <iframe
    srcdoc="{escaped}"
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
  <title>SnakePop WinPCA merged report</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    h1 {{
      margin-bottom: 10px;
    }}
    h2 {{
      margin-top: 40px;
      padding-bottom: 8px;
      border-bottom: 1px solid #ccc;
    }}
    .chrom-section {{
      background: white;
      padding: 20px;
      margin-bottom: 40px;
      border: 1px solid #ddd;
      border-radius: 8px;
    }}
    iframe {{
      display: block;
    }}
  </style>
</head>
<body>
  <h1>SnakePop WinPCA merged report</h1>
  <p>Windowed PCA plots for all chromosomes.</p>

  {''.join(sections)}

</body>
</html>
"""

        with open(output.html, "w", encoding="utf-8") as out:
            out.write(page)

rule winpca:
    input:
        f"{WINPCA_ROOT}/winpca_merged.html"
