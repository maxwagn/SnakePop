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
