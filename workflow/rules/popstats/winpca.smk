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
