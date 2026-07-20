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
