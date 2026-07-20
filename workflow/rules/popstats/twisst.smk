###############################################################################
# Twisst topology weighting from SNP-window trees
###############################################################################

TWISST_CFG = config["popstats"].get("twisst", {})
TWISST_ROOT = f"{POP_ROOT}/twisst"

TWISST_EXE = TWISST_CFG.get("executable", "bin/twisst/twisst.py")
TWISST_MAPPING_COLUMN = TWISST_CFG.get(
    "mapping_column",
    config["popstats"].get("population_column", "morphology"),
)
TWISST_OUTGROUP = TWISST_CFG.get("outgroup", None)
TWISST_METHOD = TWISST_CFG.get("method", "complete")
TWISST_ITERATIONS = TWISST_CFG.get("iterations", 100)
TWISST_ABORT_CUTOFF = TWISST_CFG.get("abort_cutoff", 100)

TWISST_RES = RES.get("twisst", {})

TWISST_TREEFILE = f"{TWISST_ROOT}/window_trees.newick"
TWISST_GROUPS_FILE = f"{TWISST_ROOT}/groups.tsv"
TWISST_TREE_METADATA = f"{TWISST_ROOT}/window_tree_metadata.tsv"
TWISST_WEIGHTS = f"{TWISST_ROOT}/twisst.weights.tsv"
TWISST_DISTANCES = f"{TWISST_ROOT}/twisst.distances.tsv"
TWISST_TOPOLOGIES = f"{TWISST_ROOT}/twisst.topologies.tsv"
TWISST_LOG = f"{TWISST_ROOT}/twisst.log"

###############################################################################
# Twisst plotting
###############################################################################

TWISST_PLOT_CFG = TWISST_CFG.get("plotting", {})
TWISST_PLOT_ENABLED = TWISST_PLOT_CFG.get("enabled", True)

TWISST_PLOTS_ROOT = f"{TWISST_ROOT}/plots"
TWISST_PLOT_PREFIX = f"{TWISST_PLOTS_ROOT}/twisst"

TWISST_PLOT_SCRIPT = TWISST_PLOT_CFG.get(
    "script",
    "bin/plot_twisst_snakepop.R",
)

TWISST_PLOT_WIDTH = TWISST_PLOT_CFG.get("width", 14)
TWISST_PLOT_HEIGHT = TWISST_PLOT_CFG.get("height", 8)
TWISST_PLOT_TOP_N = TWISST_PLOT_CFG.get("top_n", 6)
TWISST_PLOT_SMOOTH_K = TWISST_PLOT_CFG.get("smooth_k", 15)

TWISST_PLOT_FORMATS = TWISST_PLOT_CFG.get("formats", ["pdf", "svg", "png"])
TWISST_PLOT_FORMATS = [fmt.lower().strip() for fmt in TWISST_PLOT_FORMATS]

TWISST_PLOT_NAMES = [
    "summary_barplot",
    "summary_boxplot",
    "summary_topN_barplot",
    "summary_topN_boxplot",
    "topologies_only",
    "genomewide_raw_all_overlay",
    "genomewide_raw_all_stacked",
    "genomewide_smooth_all_overlay",
    "genomewide_smooth_all_stacked",
    "genomewide_raw_topN_overlay",
    "genomewide_raw_topN_stacked",
    "genomewide_smooth_topN_overlay",
    "genomewide_smooth_topN_stacked",
]

TWISST_PLOT_OUTPUTS = [
    f"{TWISST_PLOT_PREFIX}.{plot}.{fmt}"
    for plot in TWISST_PLOT_NAMES
    for fmt in TWISST_PLOT_FORMATS
] if TWISST_PLOT_ENABLED else []


rule prepare_twisst_input:
    input:
        trees = SNPTREE_ALL_TREES,
        metadata = SAMPLE_TABLE
    output:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE,
        tree_metadata = TWISST_TREE_METADATA
    run:
        import os
        import pandas as pd

        os.makedirs(TWISST_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if TWISST_MAPPING_COLUMN not in meta.columns:
            raise ValueError(
                f"Missing Twisst mapping column in metadata: {TWISST_MAPPING_COLUMN}"
            )

        mapping = meta[[SAMPLE_COL, TWISST_MAPPING_COLUMN]].dropna().copy()
        mapping[TWISST_MAPPING_COLUMN] = mapping[TWISST_MAPPING_COLUMN].map(clean_pop_name)

        if mapping.empty:
            raise ValueError(
                f"No samples found for Twisst mapping column '{TWISST_MAPPING_COLUMN}'."
            )

        mapping.to_csv(output.groups, sep="\t", index=False, header=False)

        df = pd.read_csv(input.trees, sep="\t", compression="gzip")

        if "tree" not in df.columns:
            raise ValueError("Missing tree column in SNP-tree table.")

        required = ["tree_id", "scaffold", "start", "end", "n_vcf_snps", "n_phy_sites"]
        existing = [x for x in required if x in df.columns]

        if not {"scaffold", "start", "end"}.issubset(set(existing)):
            raise ValueError(
                "Twisst plotting requires scaffold/start/end columns in SNP-tree metadata."
            )

        df[existing].to_csv(output.tree_metadata, sep="\t", index=False)

        with open(output.treefile, "w") as out:
            for tree in df["tree"].dropna():
                tree = str(tree).strip()
                if tree:
                    out.write(tree + "\n")


rule run_twisst:
    input:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE
    output:
        weights = TWISST_WEIGHTS,
        distances = TWISST_DISTANCES,
        topologies = TWISST_TOPOLOGIES,
        log = TWISST_LOG
    threads: TWISST_RES.get("threads", 1)
    resources:
        mem_mb = TWISST_RES.get("mem_mb", 8000),
        walltime = TWISST_RES.get("walltime", 8)
    params:
        exe = TWISST_EXE,
        outgroup = TWISST_OUTGROUP,
        method = TWISST_METHOD,
        iterations = TWISST_ITERATIONS,
        abort_cutoff = TWISST_ABORT_CUTOFF
    run:
        import os
        import subprocess
        import pandas as pd

        os.makedirs(TWISST_ROOT, exist_ok=True)

        groups = pd.read_csv(
            input.groups,
            sep="\t",
            header=None,
            names=["sample", "group"],
            dtype=str,
        )

        group_args = []
        for group, sub in groups.groupby("group"):
            samples = sorted(sub["sample"].dropna().unique())
            if samples:
                group_args += ["-g", group] + list(samples)

        cmd = [
            "python",
            params.exe,
            "-t", input.treefile,
            "-w", output.weights,
            "-D", output.distances,
            "--outputTopos", output.topologies,
            "--method", params.method,
            "--abortCutoff", str(params.abort_cutoff),
            "--silent",
        ] + group_args

        if params.method == "fixed":
            cmd += ["--iterations", str(params.iterations)]

        if params.outgroup:
            cmd += ["--outgroup", clean_pop_name(params.outgroup)]

        with open(output.log, "w") as log:
            log.write("COMMAND:\n")
            log.write(" ".join(cmd) + "\n\n")
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


rule plot_twisst_topologies:
    input:
        weights = TWISST_WEIGHTS,
        metadata = TWISST_TREE_METADATA,
        topologies = TWISST_TOPOLOGIES
    output:
        plots = TWISST_PLOT_OUTPUTS
    threads: 1
    resources:
        mem_mb = 4000,
        walltime = 1
    params:
        script = TWISST_PLOT_SCRIPT,
        out_prefix = TWISST_PLOT_PREFIX,
        formats = ",".join(TWISST_PLOT_FORMATS),
        width = TWISST_PLOT_WIDTH,
        height = TWISST_PLOT_HEIGHT,
        top_n = TWISST_PLOT_TOP_N,
        smooth_k = TWISST_PLOT_SMOOTH_K
    shell:
        r"""
        mkdir -p {TWISST_PLOTS_ROOT}

        Rscript {params.script} \
          --weights {input.weights} \
          --metadata {input.metadata} \
          --topologies {input.topologies} \
          --out-prefix {params.out_prefix} \
          --formats {params.formats} \
          --width {params.width} \
          --height {params.height} \
          --top-n {params.top_n} \
          --smooth-k {params.smooth_k}
        """


rule twisst:
    input:
        treefile = TWISST_TREEFILE,
        groups = TWISST_GROUPS_FILE,
        metadata = TWISST_TREE_METADATA,
        weights = TWISST_WEIGHTS,
        distances = TWISST_DISTANCES,
        topologies = TWISST_TOPOLOGIES,
        log = TWISST_LOG,
        plots = TWISST_PLOT_OUTPUTS
