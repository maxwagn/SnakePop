###############################################################################
# SNP window trees with IQ-TREE
###############################################################################

SNPTREE_CFG = config["popstats"].get("snptrees_iqtree", {})

SNPTREE_ROOT = f"{POP_ROOT}/snptrees_iqtree"

SNPTREE_SAMPLE_LIST = f"{SNPTREE_ROOT}/samples.txt"
SNPTREE_WINDOWS = f"{SNPTREE_ROOT}/windows.tsv"

SNPTREE_WINDOW_SIZE = SNPTREE_CFG.get("window_size", 1000000)
SNPTREE_INCREMENT = SNPTREE_CFG.get("increment", 1000000)
SNPTREE_MIN_SNPS = SNPTREE_CFG.get("min_snps", 500)

SNPTREE_VCF2PHYLIP = SNPTREE_CFG.get(
    "vcf2phylip",
    "bin/vcf2phylip/vcf2phylip.py",
)
SNPTREE_IQTREE = SNPTREE_CFG.get("iqtree", "iqtree")
SNPTREE_MODEL = SNPTREE_CFG.get("model", "GTR")
SNPTREE_UFBOOTS = SNPTREE_CFG.get("UFBoots", 0)
SNPTREE_FAST = SNPTREE_CFG.get("fast", True)
SNPTREE_CLEANUP = SNPTREE_CFG.get("cleanup", True)
SNPTREE_TMPDIR = SNPTREE_CFG.get("tmpdir", None)

SNPTREE_FILTER = SNPTREE_CFG.get("sample_filter", {})
SNPTREE_FILTER_IDS = SNPTREE_FILTER.get("ids", [])
SNPTREE_FILTER_COLUMN = SNPTREE_FILTER.get("column", None)
SNPTREE_FILTER_VALUES = SNPTREE_FILTER.get("values", [])

SNPTREE_RES = RES.get("snptrees_iqtree", {})

SNPTREE_CHROM_TREES = expand(
    f"{SNPTREE_ROOT}/{{chrom}}.window_trees.tsv.gz",
    chrom=CHROMOSOMES,
)

SNPTREE_ALL_TREES = f"{SNPTREE_ROOT}/all.window_trees.tsv.gz"
SNPTREE_SUMMARY = f"{SNPTREE_ROOT}/snptrees_iqtree.summary.tsv"


rule make_snptree_sample_list:
    input:
        metadata = SAMPLE_TABLE
    output:
        samples = SNPTREE_SAMPLE_LIST
    run:
        import os
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        meta = pd.read_csv(input.metadata, sep="\t", dtype=str)

        if SAMPLE_COL not in meta.columns:
            raise ValueError(f"Missing sample column in metadata: {SAMPLE_COL}")

        if SNPTREE_FILTER_IDS:
            samples = pd.Series(SNPTREE_FILTER_IDS).dropna().drop_duplicates()

        elif SNPTREE_FILTER_COLUMN is not None:
            if SNPTREE_FILTER_COLUMN not in meta.columns:
                raise ValueError(
                    f"Missing snptrees_iqtree sample_filter column: {SNPTREE_FILTER_COLUMN}"
                )

            if not SNPTREE_FILTER_VALUES:
                raise ValueError(
                    "snptrees_iqtree sample_filter.column is set, but values is empty."
                )

            samples = (
                meta.loc[
                    meta[SNPTREE_FILTER_COLUMN].isin(SNPTREE_FILTER_VALUES),
                    SAMPLE_COL,
                ]
                .dropna()
                .drop_duplicates()
            )

        else:
            samples = meta[SAMPLE_COL].dropna().drop_duplicates()

        if len(samples) < 4:
            raise ValueError(
                f"SNP tree sample list has only {len(samples)} samples. "
                "Use at least 4 samples for tree inference."
            )

        samples.to_csv(output.samples, index=False, header=False)


rule make_snptree_windows:
    input:
        fai = config["ref"]["fasta"] + ".fai"
    output:
        windows = SNPTREE_WINDOWS
    run:
        import os
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        chrom_set = set(CHROMOSOMES)
        rows = []

        with open(input.fai) as f:
            for line in f:
                chrom, length = line.split("\t")[:2]
                length = int(length)

                if chrom not in chrom_set:
                    continue

                start = 1
                while start <= length:
                    end = min(start + SNPTREE_WINDOW_SIZE - 1, length)

                    rows.append({
                        "tree_id": f"{chrom}:{start}-{end}",
                        "scaffold": chrom,
                        "start": start,
                        "end": end,
                    })

                    start += SNPTREE_INCREMENT

        pd.DataFrame(rows).to_csv(output.windows, sep="\t", index=False)


rule snptrees_iqtree_chrom:
    input:
        windows = SNPTREE_WINDOWS,
        samples = SNPTREE_SAMPLE_LIST,
        vcf = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi"
    output:
        trees = f"{SNPTREE_ROOT}/{{chrom}}.window_trees.tsv.gz"
    threads: SNPTREE_RES.get("threads", 4)
    resources:
        mem_mb = SNPTREE_RES.get("mem_mb", 16000),
        walltime = SNPTREE_RES.get("walltime", 24)
    params:
        model = SNPTREE_MODEL,
        min_snps = SNPTREE_MIN_SNPS,
        ufboots = SNPTREE_UFBOOTS,
        fast = SNPTREE_FAST,
        cleanup = SNPTREE_CLEANUP,
        tmpdir = SNPTREE_TMPDIR
    run:
        import os
        import gzip
        import glob
        import shutil
        import tempfile
        import subprocess
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        chrom = wildcards.chrom

        tmp_base = params.tmpdir
        if tmp_base is not None:
            os.makedirs(tmp_base, exist_ok=True)

        workdir = tempfile.mkdtemp(
            prefix=f"snakepop_snptrees_{chrom}_",
            dir=tmp_base,
        )

        windows = pd.read_csv(input.windows, sep="\t")
        windows = windows[windows["scaffold"].astype(str) == str(chrom)].copy()

        tmp_vcf = os.path.join(workdir, "window.vcf.gz")
        phy_prefix = "window"

        tmp_out = str(output.trees) + ".tmp"

        def clean_window_files():
            for fn in glob.glob(os.path.join(workdir, "window*")):
                try:
                    os.remove(fn)
                except FileNotFoundError:
                    pass

        def count_vcf_snps(vcf_file):
            cmd = ["bcftools", "view", "-H", vcf_file]
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
            p2 = subprocess.Popen(
                ["wc", "-l"],
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                text=True,
            )
            p1.stdout.close()
            out, _ = p2.communicate()
            return int(out.strip())

        def phylip_n_sites(phy_file):
            with open(phy_file) as f:
                first = f.readline().strip().split()

            if len(first) < 2:
                return 0

            try:
                return int(first[1])
            except ValueError:
                return 0

        try:
            with gzip.open(tmp_out, "wt") as out:
                out.write(
                    "tree_id\tscaffold\tstart\tend\tn_vcf_snps\tn_phy_sites\ttree\n"
                )

                for _, row in windows.iterrows():
                    tree_id = row["tree_id"]
                    start = int(row["start"])
                    end = int(row["end"])
                    region = f"{chrom}:{start}-{end}"

                    clean_window_files()

                    subprocess.run(
                        [
                            "bcftools", "view",
                            "-S", input.samples,
                            "-r", region,
                            "-m2", "-M2",
                            "-v", "snps",
                            "-Oz",
                            "-o", tmp_vcf,
                            input.vcf,
                        ],
                        check=True,
                    )

                    subprocess.run(
                        ["bcftools", "index", "-t", "-f", tmp_vcf],
                        check=True,
                    )

                    n_vcf_snps = count_vcf_snps(tmp_vcf)

                    if n_vcf_snps == 0:
                        clean_window_files()
                        continue

                    subprocess.run(
                        [
                            "python",
                            SNPTREE_VCF2PHYLIP,
                            "-i",
                            tmp_vcf,
                            "--output-folder",
                            workdir,
                            "--output-prefix",
                            phy_prefix,
                        ],
                        check=True,
                    )

                    phy_files = glob.glob(
                        os.path.join(workdir, f"{phy_prefix}*.phy")
                    )

                    if not phy_files:
                        clean_window_files()
                        continue

                    phy = phy_files[0]
                    n_phy_sites = phylip_n_sites(phy)

                    if n_phy_sites < params.min_snps:
                        clean_window_files()
                        continue

                    iq_prefix = os.path.join(workdir, "window.iqtree")

                    iq_cmd = [
                        SNPTREE_IQTREE,
                        "-s", phy,
                        "--seqtype", "DNA",
                        "-m", params.model,
                        "-T", "1",
                        "--prefix", iq_prefix,
                        "-quiet",
                        "-redo",
                    ]

                    if params.fast:
                        iq_cmd.append("--fast")

                    if int(params.ufboots) > 0:
                        iq_cmd += ["-B", str(params.ufboots)]

                    subprocess.run(iq_cmd, check=True)

                    treefile = f"{iq_prefix}.treefile"

                    if not os.path.exists(treefile):
                        raise ValueError(
                            f"IQ-TREE produced no treefile for {tree_id}"
                        )

                    with open(treefile) as tf:
                        tree = tf.read().strip()

                    out.write(
                        f"{tree_id}\t{chrom}\t{start}\t{end}\t"
                        f"{n_vcf_snps}\t{n_phy_sites}\t{tree}\n"
                    )

                    out.flush()
                    clean_window_files()

            os.replace(tmp_out, output.trees)

        finally:
            if params.cleanup:
                shutil.rmtree(workdir, ignore_errors=True)
                try:
                    os.remove(tmp_out)
                except FileNotFoundError:
                    pass
            else:
                print(f"Temporary SNP-tree files kept in: {workdir}")
                print(f"Temporary chromosome output kept as: {tmp_out}")


rule merge_snptrees_iqtree:
    input:
        trees = SNPTREE_CHROM_TREES
    output:
        all_trees = SNPTREE_ALL_TREES,
        summary = SNPTREE_SUMMARY
    run:
        import os
        import gzip
        import pandas as pd

        os.makedirs(SNPTREE_ROOT, exist_ok=True)

        summaries = []

        with gzip.open(output.all_trees, "wt") as out:
            out.write(
                "tree_id\tscaffold\tstart\tend\tn_vcf_snps\tn_phy_sites\ttree\n"
            )

            for fn in input.trees:
                df = pd.read_csv(fn, sep="\t", compression="gzip")

                summaries.append({
                    "file": fn,
                    "n_trees": len(df),
                    "mean_phy_sites": df["n_phy_sites"].mean() if len(df) else 0,
                    "min_phy_sites": df["n_phy_sites"].min() if len(df) else 0,
                    "max_phy_sites": df["n_phy_sites"].max() if len(df) else 0,
                    "mean_vcf_snps": df["n_vcf_snps"].mean() if len(df) else 0,
                    "min_vcf_snps": df["n_vcf_snps"].min() if len(df) else 0,
                    "max_vcf_snps": df["n_vcf_snps"].max() if len(df) else 0,
                })

                for _, row in df.iterrows():
                    out.write(
                        f"{row['tree_id']}\t{row['scaffold']}\t{row['start']}\t"
                        f"{row['end']}\t{row['n_vcf_snps']}\t"
                        f"{row['n_phy_sites']}\t{row['tree']}\n"
                    )

        pd.DataFrame(summaries).to_csv(output.summary, sep="\t", index=False)


rule snptrees_iqtree:
    input:
        samples = SNPTREE_SAMPLE_LIST,
        windows = SNPTREE_WINDOWS,
        all_trees = SNPTREE_ALL_TREES,
        summary = SNPTREE_SUMMARY
