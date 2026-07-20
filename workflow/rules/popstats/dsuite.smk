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
