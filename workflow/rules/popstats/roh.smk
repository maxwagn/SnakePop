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
