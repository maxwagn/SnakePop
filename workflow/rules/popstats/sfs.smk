###############################################################################
# SnakePop: Site Frequency Spectrum
###############################################################################

import itertools
import os
import pandas as pd


###############################################################################
# Config
###############################################################################

SFS_CFG = config["popstats"].get("sfs", {})
SFS_ENABLED = SFS_CFG.get("enabled", False)

SFS_ROOT = f"{POP_ROOT}/sfs"
SFS_SCRIPT = SFS_CFG.get("script", "bin/sfs_snakepop.py")

SFS_RES = RES.get("sfs", {})

SFS_POP_COL = SFS_CFG.get(
    "population_column",
    config["popstats"].get("population_column", "morphology"),
)

SFS_POLARISATION = SFS_CFG.get("polarisation", {}).get("mode", "folded")

if SFS_POLARISATION != "folded":
    raise ValueError(
        "SnakePop SFS currently supports only polarisation.mode: folded. "
        "Unfolded SFS will be added later."
    )

# Reuse the concatenated filtered biallelic SNP VCF produced by the PCA module.
SFS_VCF = f"{PCA_ROOT}/biallelic_snps.concat.vcf.gz"

SFS_PLOT_CFG = SFS_CFG.get("plotting", {})
SFS_FORMATS = SFS_PLOT_CFG.get("formats", ["pdf", "svg", "png"])
SFS_FORMATS = [fmt.lower().strip() for fmt in SFS_FORMATS]

SFS_REGION_CFG = SFS_CFG.get("regions", {})
SFS_USE_GENOME = SFS_REGION_CFG.get("genome", True)
SFS_USE_OUTLIERS = SFS_REGION_CFG.get("outlier_windows", False)
SFS_USE_CANDIDATES = SFS_REGION_CFG.get("candidate_regions", False)


###############################################################################
# Population and comparison helpers
###############################################################################

def sfs_all_metadata_pops():
    meta = pd.read_csv(SAMPLE_TABLE, sep="\t", dtype=str)

    if SFS_POP_COL not in meta.columns:
        raise ValueError(f"Missing SFS population column in metadata: {SFS_POP_COL}")

    return sorted(clean_pop_name(x) for x in meta[SFS_POP_COL].dropna().unique())


def sfs_included_pops():
    all_pops = sfs_all_metadata_pops()
    include_pops = SFS_CFG.get("populations", {}).get("include", "all")

    if include_pops == "all":
        return all_pops

    include_pops = [clean_pop_name(x) for x in include_pops]
    missing = sorted(set(include_pops) - set(all_pops))

    if missing:
        raise ValueError(f"SFS populations not found in metadata: {missing}")

    return [p for p in all_pops if p in include_pops]


SFS_POPS = sfs_included_pops()


def sfs_2d_pairs():
    two_d_cfg = SFS_CFG.get("spectra", {}).get("two_d", {})
    if not two_d_cfg.get("enabled", True):
        return []

    comparisons = two_d_cfg.get("comparisons", "all")

    if comparisons == "all":
        return list(itertools.combinations(SFS_POPS, 2))

    pairs = []
    for pair in comparisons:
        if len(pair) != 2:
            raise ValueError(f"Invalid SFS 2D comparison: {pair}")

        p1 = clean_pop_name(pair[0])
        p2 = clean_pop_name(pair[1])

        if p1 not in SFS_POPS:
            raise ValueError(f"SFS comparison population not included: {p1}")
        if p2 not in SFS_POPS:
            raise ValueError(f"SFS comparison population not included: {p2}")

        pairs.append((p1, p2))

    return pairs


SFS_2D_PAIRS = sfs_2d_pairs()
SFS_PAIR_MAP = {f"{p1}_{p2}": (p1, p2) for p1, p2 in SFS_2D_PAIRS}


###############################################################################
# Region helpers
###############################################################################

def sfs_region_file_for(region_class, pair):
    if region_class == "genome":
        return "none"

    if region_class == "outlier_windows":
        return f"{MAN_ROOT}/{pair}/{pair}.outlier_windows.tsv"

    if region_class == "candidate_regions":
        return f"{MAN_ROOT}/{pair}/{pair}.candidate_regions.tsv"

    raise ValueError(f"Unknown SFS region class: {region_class}")


###############################################################################
# Target construction
###############################################################################

def sfs_targets():
    if not SFS_ENABLED:
        return []

    targets = []

    one_d_enabled = SFS_CFG.get("spectra", {}).get("one_d", {}).get("enabled", True)
    two_d_enabled = SFS_CFG.get("spectra", {}).get("two_d", {}).get("enabled", True)

    if SFS_USE_GENOME:
        if one_d_enabled:
            for pop_name in SFS_POPS:
                prefix = f"{SFS_ROOT}/genome/1d/{pop_name}"
                targets.append(prefix + ".sfs.tsv")
                targets.append(prefix + ".summary.tsv")
                targets.extend([prefix + f".{fmt}" for fmt in SFS_FORMATS])

        if two_d_enabled:
            for pair_name in SFS_PAIR_MAP:
                prefix = f"{SFS_ROOT}/genome/2d/{pair_name}"
                targets.append(prefix + ".sfs.tsv")
                targets.append(prefix + ".summary.tsv")
                targets.extend([prefix + f".{fmt}" for fmt in SFS_FORMATS])

    for region_class in ["outlier_windows", "candidate_regions"]:
        if region_class == "outlier_windows" and not SFS_USE_OUTLIERS:
            continue

        if region_class == "candidate_regions" and not SFS_USE_CANDIDATES:
            continue

        for pair_name, (p1, p2) in SFS_PAIR_MAP.items():
            if one_d_enabled:
                for pop_name in [p1, p2]:
                    prefix = f"{SFS_ROOT}/{region_class}/{pair_name}/1d/{pop_name}"
                    targets.append(prefix + ".sfs.tsv")
                    targets.append(prefix + ".summary.tsv")
                    targets.extend([prefix + f".{fmt}" for fmt in SFS_FORMATS])

            if two_d_enabled:
                prefix = f"{SFS_ROOT}/{region_class}/{pair_name}/2d/{pair_name}"
                targets.append(prefix + ".sfs.tsv")
                targets.append(prefix + ".summary.tsv")
                targets.extend([prefix + f".{fmt}" for fmt in SFS_FORMATS])

    targets.append(f"{SFS_ROOT}/sfs_index.html")

    return targets


SFS_TARGETS = sfs_targets()


###############################################################################
# Rules: genome-wide SFS
###############################################################################

rule sfs_genome_1d:
    input:
        vcf = SFS_VCF,
        tbi = SFS_VCF + ".tbi",
        metadata = SAMPLE_TABLE
    output:
        table = f"{SFS_ROOT}/genome/1d/{{group}}.sfs.tsv",
        summary = f"{SFS_ROOT}/genome/1d/{{group}}.summary.tsv",
        plots = expand(
            f"{SFS_ROOT}/genome/1d/{{group}}.{{fmt}}",
            fmt = SFS_FORMATS,
            allow_missing = True,
        )
    threads: SFS_RES.get("threads", 1)
    resources:
        mem_mb = SFS_RES.get("mem_mb", 8000),
        walltime = SFS_RES.get("walltime", 6)
    params:
        script = SFS_SCRIPT,
        sample_col = SAMPLE_COL,
        pop_col = SFS_POP_COL,
        polarisation = SFS_POLARISATION,
        formats = ",".join(SFS_FORMATS),
        out_prefix = lambda wc: f"{SFS_ROOT}/genome/1d/{wc.group}"
    shell:
        r"""
        python {params.script} \
          --vcf {input.vcf} \
          --metadata {input.metadata} \
          --sample-col {params.sample_col} \
          --population-col {params.pop_col} \
          --pop1 {wildcards.group} \
          --mode 1d \
          --polarisation {params.polarisation} \
          --regions none \
          --out-prefix {params.out_prefix} \
          --formats {params.formats}
        """


rule sfs_genome_2d:
    input:
        vcf = SFS_VCF,
        tbi = SFS_VCF + ".tbi",
        metadata = SAMPLE_TABLE
    output:
        table = f"{SFS_ROOT}/genome/2d/{{pair}}.sfs.tsv",
        summary = f"{SFS_ROOT}/genome/2d/{{pair}}.summary.tsv",
        plots = expand(
            f"{SFS_ROOT}/genome/2d/{{pair}}.{{fmt}}",
            fmt = SFS_FORMATS,
            allow_missing = True,
        )
    threads: SFS_RES.get("threads", 1)
    resources:
        mem_mb = SFS_RES.get("mem_mb", 8000),
        walltime = SFS_RES.get("walltime", 6)
    params:
        script = SFS_SCRIPT,
        sample_col = SAMPLE_COL,
        pop_col = SFS_POP_COL,
        polarisation = SFS_POLARISATION,
        formats = ",".join(SFS_FORMATS),
        out_prefix = lambda wc: f"{SFS_ROOT}/genome/2d/{wc.pair}",
        pop1 = lambda wc: SFS_PAIR_MAP[wc.pair][0],
        pop2 = lambda wc: SFS_PAIR_MAP[wc.pair][1]
    shell:
        r"""
        python {params.script} \
          --vcf {input.vcf} \
          --metadata {input.metadata} \
          --sample-col {params.sample_col} \
          --population-col {params.pop_col} \
          --pop1 {params.pop1} \
          --pop2 {params.pop2} \
          --mode 2d \
          --polarisation {params.polarisation} \
          --regions none \
          --out-prefix {params.out_prefix} \
          --formats {params.formats}
        """


###############################################################################
# Rules: region-restricted SFS
###############################################################################

rule sfs_region_1d:
    input:
        vcf = SFS_VCF,
        tbi = SFS_VCF + ".tbi",
        metadata = SAMPLE_TABLE,
        regions = lambda wc: sfs_region_file_for(wc.region_class, wc.pair)
    output:
        table = f"{SFS_ROOT}/{{region_class}}/{{pair}}/1d/{{group}}.sfs.tsv",
        summary = f"{SFS_ROOT}/{{region_class}}/{{pair}}/1d/{{group}}.summary.tsv",
        plots = expand(
            f"{SFS_ROOT}/{{region_class}}/{{pair}}/1d/{{group}}.{{fmt}}",
            fmt = SFS_FORMATS,
            allow_missing = True,
        )
    threads: SFS_RES.get("threads", 1)
    resources:
        mem_mb = SFS_RES.get("mem_mb", 8000),
        walltime = SFS_RES.get("walltime", 6)
    params:
        script = SFS_SCRIPT,
        sample_col = SAMPLE_COL,
        pop_col = SFS_POP_COL,
        polarisation = SFS_POLARISATION,
        formats = ",".join(SFS_FORMATS),
        out_prefix = lambda wc: f"{SFS_ROOT}/{wc.region_class}/{wc.pair}/1d/{wc.group}"
    shell:
        r"""
        python {params.script} \
          --vcf {input.vcf} \
          --metadata {input.metadata} \
          --sample-col {params.sample_col} \
          --population-col {params.pop_col} \
          --pop1 {wildcards.group} \
          --mode 1d \
          --polarisation {params.polarisation} \
          --regions {input.regions} \
          --out-prefix {params.out_prefix} \
          --formats {params.formats}
        """


rule sfs_region_2d:
    input:
        vcf = SFS_VCF,
        tbi = SFS_VCF + ".tbi",
        metadata = SAMPLE_TABLE,
        regions = lambda wc: sfs_region_file_for(wc.region_class, wc.pair)
    output:
        table = f"{SFS_ROOT}/{{region_class}}/{{pair}}/2d/{{pair}}.sfs.tsv",
        summary = f"{SFS_ROOT}/{{region_class}}/{{pair}}/2d/{{pair}}.summary.tsv",
        plots = expand(
            f"{SFS_ROOT}/{{region_class}}/{{pair}}/2d/{{pair}}.{{fmt}}",
            fmt = SFS_FORMATS,
            allow_missing = True,
        )
    threads: SFS_RES.get("threads", 1)
    resources:
        mem_mb = SFS_RES.get("mem_mb", 8000),
        walltime = SFS_RES.get("walltime", 6)
    params:
        script = SFS_SCRIPT,
        sample_col = SAMPLE_COL,
        pop_col = SFS_POP_COL,
        polarisation = SFS_POLARISATION,
        formats = ",".join(SFS_FORMATS),
        out_prefix = lambda wc: f"{SFS_ROOT}/{wc.region_class}/{wc.pair}/2d/{wc.pair}",
        pop1 = lambda wc: SFS_PAIR_MAP[wc.pair][0],
        pop2 = lambda wc: SFS_PAIR_MAP[wc.pair][1]
    shell:
        r"""
        python {params.script} \
          --vcf {input.vcf} \
          --metadata {input.metadata} \
          --sample-col {params.sample_col} \
          --population-col {params.pop_col} \
          --pop1 {params.pop1} \
          --pop2 {params.pop2} \
          --mode 2d \
          --polarisation {params.polarisation} \
          --regions {input.regions} \
          --out-prefix {params.out_prefix} \
          --formats {params.formats}
        """


###############################################################################
# HTML index
###############################################################################

rule summarize_sfs:
    input:
        targets = SFS_TARGETS[:-1] if SFS_TARGETS else []
    output:
        html = f"{SFS_ROOT}/sfs_index.html"
    run:
        import os

        os.makedirs(SFS_ROOT, exist_ok=True)

        links = []

        for path in input.targets:
            if not str(path).endswith(".sfs.tsv"):
                continue

            rel = os.path.relpath(path, SFS_ROOT)
            summary = rel.replace(".sfs.tsv", ".summary.tsv")

            links.append(f"""
            <tr>
              <td>{rel}</td>
              <td><a href="{rel}">SFS table</a></td>
              <td><a href="{summary}">Summary</a></td>
            </tr>
            """)

        html = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>SnakePop SFS report</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 30px;
      background: #fafafa;
    }}
    table {{
      border-collapse: collapse;
      background: white;
      font-size: 12px;
    }}
    th, td {{
      border: 1px solid #ddd;
      padding: 6px 10px;
    }}
    th {{
      background: #eee;
    }}
  </style>
</head>
<body>
  <h1>SnakePop SFS report</h1>

  <p>
    Site frequency spectra were calculated from the filtered biallelic SNP VCF.
    Region-restricted spectra use Manhattan outlier/candidate-region files where enabled.
  </p>

  <table>
    <tr>
      <th>Output</th>
      <th>SFS</th>
      <th>Summary</th>
    </tr>
    {''.join(links)}
  </table>
</body>
</html>
"""

        with open(output.html, "w") as out:
            out.write(html)


###############################################################################
# Public target
###############################################################################

rule sfs:
    input:
        SFS_TARGETS

