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
