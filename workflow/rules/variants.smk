###############################################################################
# SnakePop: variant calling, filtering, and final callsets
###############################################################################

import copy
import yaml
import pysam
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt

###############################################################################
# Config
###############################################################################

REF_NAME = config["ref"]["name"]
REF_FASTA = config["ref"]["fasta"]
REF_FAI = REF_FASTA + ".fai"
CHROM_LIST = REF_FASTA + ".chromosomes.txt"

CALLSET_ID = config["callset"]["id"]
VC_ROOT = f"results/variants/{CALLSET_ID}_{REF_NAME}"
REPORT_ROOT = f"reports/filtering_stats/{CALLSET_ID}_{REF_NAME}"

ALIGN_ROOT = f"results/alignment/{REF_NAME}"

SAMPLE_TABLE = config["sample_table"]
SAMPLE_COL = config.get("sample_id_column", "id")

sample_mt = pd.read_csv(SAMPLE_TABLE, dtype=str, sep="\t")
SAMPLES = sample_mt[SAMPLE_COL].tolist()
N_SAMPLES = len(SAMPLES)

CHROMOSOMES = [line.strip() for line in open(CHROM_LIST) if line.strip()]

CHUNK_SIZE = int(config["variant_calling"].get("chunk_size", 5000000))
IND_FILTER_ID = config["variant_calling"]["ind_filter_id"]
SITE_FILTER_ID = config["variant_calling"]["site_filter_id"]
KEEP_INTERMEDIATES = config["variant_calling"].get("keep_intermediates", False)

RES = config.get("resources", {})

chrom_length = pd.read_csv(
    REF_FAI,
    sep="\t",
    usecols=[0, 1],
    names=["chrom", "len"],
    index_col=0,
).squeeze()


def maybe_temp(path):
    return path if KEEP_INTERMEDIATES else temp(path)


def get_chunks(chrom):
    return range(int(np.ceil(chrom_length[chrom] / CHUNK_SIZE)))


def region_for_chunk(wildcards):
    start = int(int(wildcards.chunk) * CHUNK_SIZE + 1)
    end = int((int(wildcards.chunk) + 1) * CHUNK_SIZE)
    return f"{wildcards.chrom}:{start}-{end}"


def bam_files(wildcards):
    return expand(
        f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam",
        sample=SAMPLES,
    )


def get_total_dp_files(wildcards):
    return [
        f"{VC_ROOT}/coverage/total_coverage_{chrom}.{chunk}.txt"
        for chrom in CHROMOSOMES
        for chunk in get_chunks(chrom)
    ]


def read_dp(fn):
    dp_dic = {}
    with open(fn) as f:
        for line in f:
            sample, dp = line.strip().split()
            dp_dic[sample] = float(dp)
    return dp_dic


def write_dp(fn, samples, dps):
    with open(fn, "w") as f:
        for sample, dp in zip(samples, dps):
            f.write(f"{sample}\t{dp}\n")


def get_filter_command(wildcards):
    commands = []
    filters = config["site_filter_sets"][wildcards.site_filter_id]["filters"]

    for filter_name, fd in filters.items():
        assert fd["threshold_type"] == "absolute"

        commands.append(
            "bcftools filter "
            f"--soft-filter {filter_name} "
            "--mode + "
            "-O u "
            f"--exclude '{fd['tag']} {fd['operator']} {fd['threshold']}'"
        )

    return " | ".join(commands)


RAW_BCF_TARGETS = [
    f"{VC_ROOT}/bcf/raw/all_sites.{chrom}.{chunk}.bcf.csi"
    for chrom in CHROMOSOMES
    for chunk in get_chunks(chrom)
]

FILTERED_CHROM_TARGETS = [
    f"{VC_ROOT}/bcf/final/all_sites.{IND_FILTER_ID}.{SITE_FILTER_ID}.{chrom}.bcf.csi"
    for chrom in CHROMOSOMES
]

FINAL_VCF_TARGETS = (
    expand(
        f"{VC_ROOT}/vcf/all_sites.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi",
        chrom=CHROMOSOMES,
    )
    + expand(
        f"{VC_ROOT}/vcf/variants.{IND_FILTER_ID}.{SITE_FILTER_ID}.{{chrom}}.vcf.gz.tbi",
        chrom=CHROMOSOMES,
    )
    + expand(
        f"{VC_ROOT}/vcf/variants.{IND_FILTER_ID}.{SITE_FILTER_ID}.biallelic_snps.{{chrom}}.vcf.gz.tbi",
        chrom=CHROMOSOMES,
    )
)

###############################################################################
# Public targets
###############################################################################

rule raw_calling:
    input:
        RAW_BCF_TARGETS


rule filtering:
    input:
        FILTERED_CHROM_TARGETS,
        f"{REPORT_ROOT}/genotype_filter_summary.tsv",
        f"{REPORT_ROOT}/genotype_filter_summary.png",
        f"{REPORT_ROOT}/site_filter_summary.tsv",
        f"{REPORT_ROOT}/site_filter_summary.png"


rule final_callset:
    input:
        FINAL_VCF_TARGETS


rule variants:
    input:
        rules.final_callset.input


###############################################################################
# Raw all-sites calling
###############################################################################

rule call_all_sites_bcf:
    input:
        bams = bam_files,
        ref = REF_FASTA,
        fai = REF_FAI
    output:
        bcf = maybe_temp(f"{VC_ROOT}/bcf/raw/all_sites.{{chrom}}.{{chunk}}.bcf"),
        csi = maybe_temp(f"{VC_ROOT}/bcf/raw/all_sites.{{chrom}}.{{chunk}}.bcf.csi"),
        tot_dp = f"{VC_ROOT}/coverage/total_coverage_{{chrom}}.{{chunk}}.txt"
    threads: RES.get("call_all_sites_bcf", {}).get("threads", 2)
    resources:
        mem_mb = RES.get("call_all_sites_bcf", {}).get("mem_mb", 18000),
        walltime = RES.get("call_all_sites_bcf", {}).get("walltime", 24)
    params:
        region = region_for_chunk,
        max_depth = config["variant_calling"].get("max_depth", 3000)
    shell:
        r"""
        bcftools mpileup \
          -a FORMAT/AD,FORMAT/DP \
          -d {params.max_depth} \
          --threads {threads} \
          -Ou \
          --regions {params.region} \
          -f {input.ref} \
          {input.bams} \
        | bcftools call \
          --threads {threads} \
          -f GQ \
          -m \
          -Ou \
        | bcftools +fill-tags \
          --threads {threads} \
          -Ob \
        | tee {output.bcf} \
        | bcftools query -f '[%DP\t]\n' \
        | awk '{{for(i=1; i<=NF; i++) dp_tot[i]+=$i}} END{{for(i=1; i<=NF; i++) {{printf dp_tot[i]"\t"}}; printf "\n"}}' \
          > {output.tot_dp}

        bcftools index -f {output.bcf}
        """


###############################################################################
# Mean coverage and depth thresholds
###############################################################################

rule get_mean_coverage:
    input:
        total_dps = get_total_dp_files
    output:
        cov = f"{VC_ROOT}/coverage/mean_coverage.txt"
    run:
        mean_dps = np.zeros(N_SAMPLES)

        for fn in input.total_dps:
            with open(fn) as f:
                mean_dps += np.array([float(x) for x in f.readline().strip().split()])

        mean_dps = mean_dps / chrom_length.loc[CHROMOSOMES].sum()
        write_dp(output.cov, SAMPLES, mean_dps)


rule min_max_dp:
    input:
        cov = f"{VC_ROOT}/coverage/mean_coverage.txt"
    output:
        min_dp = f"{VC_ROOT}/coverage/low_coverage.{{ind_filter_id}}.txt",
        max_dp = f"{VC_ROOT}/coverage/excess_coverage.{{ind_filter_id}}.txt"
    params:
        filter_thresh = lambda wc: config["individual_filter_sets"][wc.ind_filter_id]
    run:
        mean_dp_dic = read_dp(input.cov)

        min_dp = []
        max_dp = []

        for dp in mean_dp_dic.values():
            min_dp.append(stats.poisson.ppf(params.filter_thresh["min_dp_to_missing_pval"], dp))
            max_dp.append(stats.poisson.isf(params.filter_thresh["max_dp_to_missing_pval"], dp))

        write_dp(output.min_dp, mean_dp_dic.keys(), min_dp)
        write_dp(output.max_dp, mean_dp_dic.keys(), max_dp)


###############################################################################
# Genotype-level filtering
###############################################################################

rule genotype_filter_all_sites:
    input:
        bcf = f"{VC_ROOT}/bcf/raw/all_sites.{{chrom}}.{{chunk}}.bcf",
        cov = f"{VC_ROOT}/coverage/mean_coverage.txt",
        min_dp = f"{VC_ROOT}/coverage/low_coverage.{{ind_filter_id}}.txt",
        max_dp = f"{VC_ROOT}/coverage/excess_coverage.{{ind_filter_id}}.txt"
    output:
        bcf = maybe_temp(f"{VC_ROOT}/bcf/genotype_filtered/all_sites.{{ind_filter_id}}.{{chrom}}.{{chunk}}.bcf"),
        csi = maybe_temp(f"{VC_ROOT}/bcf/genotype_filtered/all_sites.{{ind_filter_id}}.{{chrom}}.{{chunk}}.bcf.csi"),
        filter_per_sample = f"{REPORT_ROOT}/per_chunk/filter_per_sample.{{ind_filter_id}}.{{chrom}}.{{chunk}}.tsv",
        total_per_sample = f"{REPORT_ROOT}/per_chunk/total_per_sample.{{ind_filter_id}}.{{chrom}}.{{chunk}}.tsv"
    threads: RES.get("genotype_filter_all_sites", {}).get("threads", 1)
    resources:
        mem_mb = RES.get("genotype_filter_all_sites", {}).get("mem_mb", 36000),
        walltime = RES.get("genotype_filter_all_sites", {}).get("walltime", 12)
    params:
        filter_thresh = lambda wc: config["individual_filter_sets"][wc.ind_filter_id]
    run:
        mean_dp_dic = read_dp(input.cov)
        min_dp_dic = read_dp(input.min_dp)
        max_dp_dic = read_dp(input.max_dp)

        bcf_in = pysam.VariantFile(input.bcf)

        samples_in_vcf = list(bcf_in.header.samples)
        assert samples_in_vcf == SAMPLES, "Sample order mismatch between metadata and BCF."

        for tag, number, typ, desc in [
            ("DD", 1, "Float", "Average normalized deviation from mean sample coverage."),
            ("ExcHetOrig", "A", "Float", "Original ExcessHet value before filtering."),
            ("AB_Het", 1, "Float", "Phred-scaled allele-balance test statistic."),
            ("MF", 1, "Float", "Fraction of missing genotypes after filtering."),
        ]:
            if tag not in bcf_in.header.info:
                bcf_in.header.info.add(tag, number, typ, desc)

        filltag_stream = subprocess.Popen(
            ["bcftools", "+fill-tags", "-Ob"],
            stdin=subprocess.PIPE,
            stdout=open(output.bcf, "wb"),
            stderr=subprocess.PIPE,
        )

        bcf_out = pysam.VariantFile(filltag_stream.stdin, "wu", header=bcf_in.header)

        filter_per_sample = {
            s: {"allele_balance": 0, "low_depth": 0, "high_depth": 0, "low_gq": 0}
            for s in SAMPLES
        }

        total_per_sample = copy.deepcopy(filter_per_sample)

        for rec in bcf_in:
            deviations = []
            ads_hets = np.array([0, 0], dtype=int)
            missing = 0

            for sample_name, sample_call in rec.samples.items():
                dp = sample_call.get("DP")

                if dp is not None:
                    if dp <= min_dp_dic[sample_name]:
                        sample_call["GT"] = (None, None)
                        filter_per_sample[sample_name]["low_depth"] += 1
                    total_per_sample[sample_name]["low_depth"] += 1

                    if dp >= max_dp_dic[sample_name]:
                        sample_call["GT"] = (None, None)
                        filter_per_sample[sample_name]["high_depth"] += 1
                    total_per_sample[sample_name]["high_depth"] += 1

                    if mean_dp_dic[sample_name] > 0:
                        deviations.append(abs(dp - mean_dp_dic[sample_name]) / np.sqrt(mean_dp_dic[sample_name]))

                gq = sample_call.get("GQ")

                if gq is not None:
                    if gq < params.filter_thresh["min_GQ"]:
                        sample_call["GT"] = (None, None)
                        filter_per_sample[sample_name]["low_gq"] += 1
                    total_per_sample[sample_name]["low_gq"] += 1

                gt = sample_call.get("GT")

                if gt is not None and None not in gt and len(gt) == 2 and gt[0] != gt[1]:
                    ad_raw = sample_call.get("AD")

                    if ad_raw is not None:
                        try:
                            ad = np.array(ad_raw, dtype=int)[list(gt)]
                            ad_sum = int(np.sum(ad))

                            if ad_sum > 0:
                                ads_hets += ad

                                ab_pval = stats.binomtest(int(ad[0]), n=ad_sum, p=0.5).pvalue

                                if ab_pval <= params.filter_thresh["ab_test_thresh"]:
                                    sample_call["GT"] = (None, None)
                                    filter_per_sample[sample_name]["allele_balance"] += 1

                            total_per_sample[sample_name]["allele_balance"] += 1

                        except Exception:
                            pass

                if sample_call.get("GT") == (None, None):
                    missing += 1

            rec.info["MF"] = missing / len(SAMPLES)

            if deviations:
                rec.info["DD"] = float(np.mean(deviations))

            if int(np.sum(ads_hets)) > 0:
                ab_pval = stats.binomtest(
                    int(ads_hets[0]),
                    n=int(np.sum(ads_hets)),
                    p=0.5,
                ).pvalue

                if ab_pval > 0:
                    rec.info["AB_Het"] = float(-10 * np.log10(ab_pval))

                if "ExcHet" in rec.info:
                    rec.info["ExcHetOrig"] = rec.info["ExcHet"]

            bcf_out.write(rec)

        bcf_out.close()
        _, stderr = filltag_stream.communicate()

        if filltag_stream.returncode != 0:
            raise RuntimeError(stderr.decode())

        shell(f"bcftools index -f --csi {output.bcf}")

        pd.DataFrame(filter_per_sample).to_csv(output.filter_per_sample, sep="\t")
        pd.DataFrame(total_per_sample).to_csv(output.total_per_sample, sep="\t")


###############################################################################
# Site-level filtering
###############################################################################

rule add_site_filters:
    input:
        bcf = f"{VC_ROOT}/bcf/genotype_filtered/all_sites.{{ind_filter_id}}.{{chrom}}.{{chunk}}.bcf"
    output:
        bcf = maybe_temp(f"{VC_ROOT}/bcf/site_filtered/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.{{chunk}}.bcf"),
        csi = maybe_temp(f"{VC_ROOT}/bcf/site_filtered/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.{{chunk}}.bcf.csi")
    params:
        filter_command = get_filter_command
    shell:
        r"""
        bcftools view -O u {input.bcf} \
        | {params.filter_command} \
        | bcftools view -O b -o {output.bcf}

        bcftools index -f {output.bcf}
        """


rule get_filter_info:
    input:
        bcf = f"{VC_ROOT}/bcf/site_filtered/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.{{chunk}}.bcf"
    output:
        info = f"{REPORT_ROOT}/per_chunk/site_filters.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.{{chunk}}.tsv"
    shell:
        r"""
        bcftools query \
          --include 'FILTER != "PASS"' \
          --format '%CHROM\t%POS\t%FILTER\n' \
          {input.bcf} \
          > {output.info}
        """


###############################################################################
# Chromosome-level final BCFs
###############################################################################

rule make_chrom_all_sites:
    input:
        bcfs = lambda wc: expand(
            f"{VC_ROOT}/bcf/site_filtered/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.{{chunk}}.bcf",
            ind_filter_id=wc.ind_filter_id,
            site_filter_id=wc.site_filter_id,
            chrom=wc.chrom,
            chunk=get_chunks(wc.chrom),
        )
    output:
        bcf = f"{VC_ROOT}/bcf/final/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf",
        csi = f"{VC_ROOT}/bcf/final/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf.csi"
    shell:
        r"""
        bcftools concat -Ou {input.bcfs} \
        | bcftools view --apply-filters PASS -Ob -o {output.bcf}

        bcftools index -f {output.bcf}
        """


rule make_chrom_variants:
    input:
        bcf = f"{VC_ROOT}/bcf/final/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf"
    output:
        bcf = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf",
        csi = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf.csi"
    shell:
        r"""
        bcftools view \
          --exclude 'ALT="."' \
          -Ob \
          -o {output.bcf} \
          {input.bcf}

        bcftools index -f {output.bcf}
        """


rule make_chrom_biallelic_snps:
    input:
        bcf = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf"
    output:
        bcf = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.biallelic_snps.{{chrom}}.bcf",
        csi = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.biallelic_snps.{{chrom}}.bcf.csi"
    shell:
        r"""
        bcftools view \
          --max-alleles 2 \
          --types snps \
          -Ob \
          -o {output.bcf} \
          {input.bcf}

        bcftools index -f {output.bcf}
        """


###############################################################################
# Export final VCFs
###############################################################################

rule export_all_sites_vcf:
    input:
        bcf = f"{VC_ROOT}/bcf/final/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf"
    output:
        vcf = f"{VC_ROOT}/vcf/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/all_sites.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.vcf.gz.tbi"
    shell:
        r"""
        bcftools view -Oz -o {output.vcf} {input.bcf}
        bcftools index --tbi -f {output.vcf}
        """


rule export_variants_vcf:
    input:
        bcf = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.bcf"
    output:
        vcf = f"{VC_ROOT}/vcf/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/variants.{{ind_filter_id}}.{{site_filter_id}}.{{chrom}}.vcf.gz.tbi"
    shell:
        r"""
        bcftools view -Oz -o {output.vcf} {input.bcf}
        bcftools index --tbi -f {output.vcf}
        """


rule export_biallelic_snps_vcf:
    input:
        bcf = f"{VC_ROOT}/bcf/final/variants.{{ind_filter_id}}.{{site_filter_id}}.biallelic_snps.{{chrom}}.bcf"
    output:
        vcf = f"{VC_ROOT}/vcf/variants.{{ind_filter_id}}.{{site_filter_id}}.biallelic_snps.{{chrom}}.vcf.gz",
        tbi = f"{VC_ROOT}/vcf/variants.{{ind_filter_id}}.{{site_filter_id}}.biallelic_snps.{{chrom}}.vcf.gz.tbi"
    shell:
        r"""
        bcftools view -Oz -o {output.vcf} {input.bcf}
        bcftools index --tbi -f {output.vcf}
        """


###############################################################################
# Filtering reports
###############################################################################

rule aggregate_genotype_filter_stats:
    input:
        filter_files = [
            f"{REPORT_ROOT}/per_chunk/filter_per_sample.{IND_FILTER_ID}.{chrom}.{chunk}.tsv"
            for chrom in CHROMOSOMES
            for chunk in get_chunks(chrom)
        ]
    output:
        tsv = f"{REPORT_ROOT}/genotype_filter_summary.tsv",
        png = f"{REPORT_ROOT}/genotype_filter_summary.png"
    run:
        totals = None

        for fn in input.filter_files:
            df = pd.read_csv(fn, sep="\t", index_col=0)
            totals = df if totals is None else totals.add(df, fill_value=0)

        totals.to_csv(output.tsv, sep="\t")

        by_filter = totals.sum(axis=1).sort_values(ascending=False)

        plt.figure(figsize=(8, 5))
        by_filter.plot(kind="bar")
        plt.ylabel("Number of genotypes masked")
        plt.tight_layout()
        plt.savefig(output.png, dpi=200)
        plt.close()


rule aggregate_site_filter_stats:
    input:
        info_files = [
            f"{REPORT_ROOT}/per_chunk/site_filters.{IND_FILTER_ID}.{SITE_FILTER_ID}.{chrom}.{chunk}.tsv"
            for chrom in CHROMOSOMES
            for chunk in get_chunks(chrom)
        ]
    output:
        tsv = f"{REPORT_ROOT}/site_filter_summary.tsv",
        png = f"{REPORT_ROOT}/site_filter_summary.png"
    run:
        counts = {}

        for fn in input.info_files:
            with open(fn) as f:
                for line in f:
                    fields = line.strip().split("\t")
                    if len(fields) < 3:
                        continue

                    filters = fields[2].split(";")
                    for flt in filters:
                        counts[flt] = counts.get(flt, 0) + 1

        df = pd.DataFrame(
            sorted(counts.items(), key=lambda x: x[1], reverse=True),
            columns=["filter", "count"],
        )

        df.to_csv(output.tsv, sep="\t", index=False)

        if len(df) > 0:
            plt.figure(figsize=(8, 5))
            plt.bar(df["filter"], df["count"])
            plt.ylabel("Number of filtered sites")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.savefig(output.png, dpi=200)
            plt.close()
        else:
            plt.figure(figsize=(6, 4))
            plt.text(0.5, 0.5, "No filtered sites", ha="center", va="center")
            plt.axis("off")
            plt.savefig(output.png, dpi=200)
            plt.close()

###############################################################################
# Cleanup
###############################################################################

rule clean_intermediates:
    message:
        "Removing intermediate variant-calling BCF folders."
    shell:
        r"""
        rm -rf {VC_ROOT}/bcf/raw
        rm -rf {VC_ROOT}/bcf/genotype_filtered
        rm -rf {VC_ROOT}/bcf/site_filtered
        """
