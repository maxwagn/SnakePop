from pathlib import Path

REF_NAME = config["ref"]["name"]
CALLSET_ID = config["callset"]["id"]

IND_FILTER_ID = config["variant_calling"]["ind_filter_id"]
SITE_FILTER_ID = config["variant_calling"]["site_filter_id"]

VC_ROOT = f"results/variants/{CALLSET_ID}_{REF_NAME}"
POP_ROOT = f"results/popstats/{CALLSET_ID}_{REF_NAME}"
PCA_ROOT = f"{POP_ROOT}/pca"
REPORT_ROOT = f"reports/popstats/{CALLSET_ID}_{REF_NAME}/pca"

CHROM_LIST = config["ref"]["fasta"] + ".chromosomes.txt"

if Path(CHROM_LIST).exists():
    with open(CHROM_LIST) as handle:
        CHROMOSOMES = [
            line.strip()
            for line in handle
            if line.strip()
        ]
else:
    CHROMOSOMES = []

SAMPLE_TABLE = config["sample_table"]
SAMPLE_COL = config.get("sample_id_column", "id")

RES = config.get("resources", {})
PCA_RES = RES.get("pca", {})
PCA_CFG = config.get("popstats", {}).get("pca", {})

N_PCS = PCA_CFG.get("n_pcs", 10)
MAF = PCA_CFG.get("maf", 0.05)
MIND = PCA_CFG.get("mind", 0.2)
GENO = PCA_CFG.get("geno", 0.2)
COLOR_BY = PCA_CFG.get("color_by", "species")
LABEL_SAMPLES = PCA_CFG.get("label_samples", False)
LD_PRUNE = PCA_CFG.get("ld_prune", False)

LD = PCA_CFG.get("indep_pairwise", {})
LD_WINDOW = LD.get("window_kb", 50)
LD_STEP = LD.get("step", 5)
LD_R2 = LD.get("r2", 0.2)

PLINK_CHR_SET = PCA_CFG.get("plink_chr_set", "23 no-xy no-mt")

BIALLELIC_SNPS = expand(
    f"{VC_ROOT}/vcf/biallelic_snps.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz",
    chrom=CHROMOSOMES,
)
