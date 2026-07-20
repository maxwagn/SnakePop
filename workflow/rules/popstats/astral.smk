###############################################################################
# ASTRAL summary tree from SNP window trees
###############################################################################

ASTRAL_CFG = config["popstats"].get("astral", {})
ASTRAL_ROOT = f"{POP_ROOT}/astral"

ASTRAL_EXE = ASTRAL_CFG.get("executable", "astral")
ASTRAL_MAPPING_COLUMN = ASTRAL_CFG.get(
    "mapping_column",
    config["popstats"].get("population_column", "morphology"),
)
ASTRAL_BRANCH_ANNOTATE = ASTRAL_CFG.get("branch_annotate", 3)

ASTRAL_ROOT_CFG = ASTRAL_CFG.get("root", {})
ASTRAL_DO_ROOT = ASTRAL_ROOT_CFG.get("enabled", True)
ASTRAL_OUTGROUPS = ASTRAL_ROOT_CFG.get("outgroups", [])

ASTRAL_RES = RES.get("astral", {})

ASTRAL_GENE_TREES = f"{ASTRAL_ROOT}/window_trees.newick"
ASTRAL_MAPPING = f"{ASTRAL_ROOT}/astral_mapping.tsv"
ASTRAL_TREE = f"{ASTRAL_ROOT}/astral.tree"
ASTRAL_TOPOLOGY = f"{ASTRAL_ROOT}/astral.topology.tree"
ASTRAL_LOG = f"{ASTRAL_ROOT}/astral.log"


rule prepare_astral_input:
    input:
        trees = SNPTREE_ALL_TREES,
        metadata = SAMPLE_TABLE
    output:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING
    run:
        import os
        import pandas as pd

        os.makedirs(ASTRAL_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if ASTRAL_MAPPING_COLUMN not in meta.columns:
            raise ValueError(
                f"Missing ASTRAL mapping column in metadata: {ASTRAL_MAPPING_COLUMN}"
            )

        meta = meta[[SAMPLE_COL, ASTRAL_MAPPING_COLUMN]].dropna().copy()
        meta[ASTRAL_MAPPING_COLUMN] = meta[ASTRAL_MAPPING_COLUMN].map(clean_pop_name)

        # ASTRAL mapping format:
        # species: sample1,sample2,sample3
        with open(output.mapping, "w") as out:
            for group, sub in meta.groupby(ASTRAL_MAPPING_COLUMN):
                group_samples = sorted(sub[SAMPLE_COL].dropna().unique())
                if group_samples:
                    out.write(f"{group}: {','.join(group_samples)}\n")

        df = pd.read_csv(input.trees, sep="\t", compression="gzip")

        if "tree" not in df.columns:
            raise ValueError("Missing tree column in SNP-tree table.")

        with open(output.gene_trees, "w") as out:
            for tree in df["tree"].dropna():
                tree = str(tree).strip()
                if tree:
                    out.write(tree + "\n")


rule astral_tree:
    input:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING
    output:
        tree = ASTRAL_TREE,
        log = ASTRAL_LOG
    threads: ASTRAL_RES.get("threads", 1)
    resources:
        mem_mb = ASTRAL_RES.get("mem_mb", 16000),
        walltime = ASTRAL_RES.get("walltime", 4)
    params:
        exe = ASTRAL_EXE,
        branch_annotate = ASTRAL_BRANCH_ANNOTATE,
        do_root = ASTRAL_DO_ROOT,
        outgroups = ASTRAL_OUTGROUPS
    run:
        import os
        import subprocess

        os.makedirs(ASTRAL_ROOT, exist_ok=True)

        cmd = [
            params.exe,
            "-i", input.gene_trees,
            "-a", input.mapping,
            "-o", output.tree,
            "-t", str(params.branch_annotate),
        ]

        if params.do_root and params.outgroups:
            cmd += ["--outgroup", ",".join(params.outgroups)]

        with open(output.log, "w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)



rule astral_topology_tree:
    input:
        tree = ASTRAL_TREE
    output:
        topology = ASTRAL_TOPOLOGY
    params:
        do_root = ASTRAL_DO_ROOT,
        outgroups = ASTRAL_OUTGROUPS
    run:
        import re
        from ete3 import Tree

        with open(input.tree) as f:
            newick = f.read().strip()

        # Remove ASTRAL bracket annotations, e.g. [q1=...;q2=...]
        newick = re.sub(r"\[[^\[\]]*\]", "", newick)

        # Parse cleaned Newick
        t = Tree(newick, format=1)

        if params.do_root and params.outgroups:
            outgroup = params.outgroups[0]

            try:
                t.set_outgroup(t & outgroup)
            except Exception:
                tips = sorted([leaf.name for leaf in t.iter_leaves()])
                raise ValueError(
                    f"Outgroup '{outgroup}' not found in ASTRAL tree. "
                    f"Available tips: {tips}"
                )

        # Write plain topology only: no branch lengths, no support values
        topology = t.write(format=9).strip()

        # Optional: normalize order so outgroup is visibly first
        if params.do_root and params.outgroups:
            outgroup = params.outgroups[0]
            if not topology.startswith(f"({outgroup},"):
                # topology is still correctly rooted, this only affects display order
                pass

        with open(output.topology, "w") as out:
            out.write(topology + "\n")

rule astral:
    input:
        gene_trees = ASTRAL_GENE_TREES,
        mapping = ASTRAL_MAPPING,
        tree = ASTRAL_TREE,
        topology = ASTRAL_TOPOLOGY,
        log = ASTRAL_LOG
