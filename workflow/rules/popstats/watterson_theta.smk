###############################################################################
# SnakePop: Watterson's theta
###############################################################################

import pandas as pd


###############################################################################
# Configuration
###############################################################################

WAT_CFG = (
    config["popstats"]
    .get(
        "watterson_theta",
        {},
    )
)

WAT_ENABLED = WAT_CFG.get(
    "enabled",
    False,
)

WAT_ROOT = (
    f"{POP_ROOT}/watterson_theta"
)

WAT_SCRIPT = WAT_CFG.get(
    "script",
    "bin/watterson_theta_snakepop.py",
)

WAT_RES = RES.get(
    "watterson_theta",
    {},
)

WAT_POP_COL = WAT_CFG.get(
    "population_column",
    config[
        "popstats"
    ].get(
        "population_column",
        "morphology",
    ),
)

WAT_WINDOW_SIZE = int(
    WAT_CFG.get(
        "window_size",
        60000,
    )
)

WAT_WINDOW_STEP = int(
    WAT_CFG.get(
        "window_step",
        30000,
    )
)

WAT_PLOIDY = int(
    WAT_CFG.get(
        "ploidy",
        2,
    )
)

WAT_MIN_CALLED = float(
    WAT_CFG.get(
        "min_called_fraction",
        0.8,
    )
)

WAT_MIN_CALLABLE = int(
    WAT_CFG.get(
        "min_callable_sites",
        1000,
    )
)

WAT_SNPS_ONLY = WAT_CFG.get(
    "snps_only",
    True,
)


###############################################################################
# Population configuration
###############################################################################

WAT_INCLUDE = (
    WAT_CFG
    .get(
        "populations",
        {},
    )
    .get(
        "include",
        "all",
    )
)


def watterson_population_argument():

    if WAT_INCLUDE == "all":
        return "all"

    return ",".join(
        clean_pop_name(x)
        for x in WAT_INCLUDE
    )


WAT_POP_ARG = (
    watterson_population_argument()
)


###############################################################################
# Input VCF
###############################################################################

def watterson_vcf(wc):

    return (
        f"{VC_ROOT}/vcf/"
        f"all_sites."
        f"{IND_FILTER_ID}."
        f"{SITE_FILTER_ID}."
        f"{wc.chrom}.vcf.gz"
    )


###############################################################################
# Targets
###############################################################################

WAT_CHROM_TARGETS = (
    expand(
        f"{WAT_ROOT}/per_chrom/"
        "{{chrom}}.tsv",
        chrom=CHROMOSOMES,
    )
    if WAT_ENABLED
    else []
)


WAT_TARGETS = (
    [
        f"{WAT_ROOT}/watterson_theta.tsv"
    ]
    if WAT_ENABLED
    else []
)


###############################################################################
# Per-chromosome calculation
###############################################################################

rule watterson_theta_chrom:

    input:
        vcf=watterson_vcf,

        tbi=lambda wc:
            watterson_vcf(wc)
            + ".tbi",

        metadata=SAMPLE_TABLE

    output:
        table=(
            f"{WAT_ROOT}/per_chrom/"
            "{{chrom}}.tsv"
        )

    threads:
        WAT_RES.get(
            "threads",
            1,
        )

    resources:
        mem_mb=WAT_RES.get(
            "mem_mb",
            8000,
        ),

        walltime=WAT_RES.get(
            "walltime",
            6,
        )

    params:
        script=WAT_SCRIPT,

        sample_col=SAMPLE_COL,

        population_col=WAT_POP_COL,

        populations=WAT_POP_ARG,

        window_size=WAT_WINDOW_SIZE,

        window_step=WAT_WINDOW_STEP,

        ploidy=WAT_PLOIDY,

        min_called=WAT_MIN_CALLED,

        min_callable=WAT_MIN_CALLABLE,

        snps_only=(
            "--snps-only"
            if WAT_SNPS_ONLY
            else ""
        )

    shell:
        r"""
        python {params.script} \
          --vcf {input.vcf} \
          --metadata {input.metadata} \
          --chrom {wildcards.chrom} \
          --sample-col {params.sample_col} \
          --population-col {params.population_col} \
          --populations {params.populations} \
          --window-size {params.window_size} \
          --window-step {params.window_step} \
          --ploidy {params.ploidy} \
          --min-called-fraction {params.min_called} \
          --min-callable-sites {params.min_callable} \
          {params.snps_only} \
          --output {output.table}
        """


###############################################################################
# Merge chromosomes
###############################################################################

rule merge_watterson_theta:

    input:
        WAT_CHROM_TARGETS

    output:
        table=(
            f"{WAT_ROOT}/"
            "watterson_theta.tsv"
        )

    run:

        import os
        import pandas as pd

        os.makedirs(
            WAT_ROOT,
            exist_ok=True,
        )

        if not input:
            raise ValueError(
                "No chromosome-level "
                "Watterson theta files found."
            )

        frames = [
            pd.read_csv(
                fn,
                sep="\t",
            )
            for fn in input
        ]

        result = pd.concat(
            frames,
            ignore_index=True,
        )

        result.to_csv(
            output.table,
            sep="\t",
            index=False,
            na_rep="NA",
        )


###############################################################################
# Public target
###############################################################################

rule watterson_theta:

    input:
        WAT_TARGETS
