###############################################################################
# SnakePop: reference preparation and alignment rules
###############################################################################

import gzip
import pandas as pd

###############################################################################
# Config
###############################################################################

REF_NAME = config["ref"]["name"]
REF_ORIGINAL = config["ref"]["original_fasta"]
REF_FASTA = config["ref"]["fasta"]
ASSEMBLY_REPORT = config["ref"]["assembly_report"]

ALIGN_ROOT = f"results/alignment/{REF_NAME}"

SAMPLE_TABLE = config["sample_table"]
SAMPLE_COL = config.get("sample_id_column", "id")
READ_COL = config.get("read_files_column", "read_files")

RES = config.get("resources", {})

###############################################################################
# Sample metadata
###############################################################################

sample_mt = pd.read_csv(SAMPLE_TABLE, dtype=str, sep="\t")
SAMPLES = sample_mt[SAMPLE_COL].tolist()

READS = {}
for _, row in sample_mt.iterrows():
    reads = row[READ_COL].split(";")
    if len(reads) != 2:
        raise ValueError(
            f"Expected two read files separated by ';' for sample {row[SAMPLE_COL]}"
        )

    READS[row[SAMPLE_COL]] = {
        "r1": reads[0],
        "r2": reads[1],
    }


def read1(wildcards):
    return READS[wildcards.sample]["r1"]


def read2(wildcards):
    return READS[wildcards.sample]["r2"]


###############################################################################
# Reference preparation
###############################################################################

rule rename_reference:
    input:
        ref = REF_ORIGINAL,
        report = ASSEMBLY_REPORT
    output:
        ref = REF_FASTA,
        chrom_map = REF_FASTA + ".chrom_map.tsv",
        chrom_list = REF_FASTA + ".chromosomes.txt"
    threads: RES.get("rename_reference", {}).get("threads", 1)
    resources:
        mem_mb = RES.get("rename_reference", {}).get("mem_mb", 4000),
        walltime = RES.get("rename_reference", {}).get("walltime", 1)
    params:
        chrom_prefix = config["ref"].get("chrom_prefix", "chr"),
        mt_name = config["ref"].get("mt_name", "mtDNA")
    run:
        report = pd.read_csv(input.report, sep="\t", dtype=str).fillna("")

        rename = {}
        chromosomes = []

        for _, row in report.iterrows():
            old = row["GenBank seq accession"]
            role = row["Role"]
            mol_type = row["Molecule type"]
            chrom = row["Chromosome name"]
            seq_name = row["Sequence name"]

            if not old:
                continue

            if role == "assembled-molecule" and mol_type == "Chromosome":
                new = f"{params.chrom_prefix}{chrom}"
                chromosomes.append(new)

            elif mol_type.lower().startswith("mitochond"):
                new = params.mt_name

            else:
                new = seq_name

            rename[old] = new

        with open(output.chrom_map, "w") as f:
            for old, new in rename.items():
                f.write(f"{old}\t{new}\n")

        with open(output.chrom_list, "w") as f:
            for chrom in chromosomes:
                f.write(f"{chrom}\n")

        opener = gzip.open if str(input.ref).endswith(".gz") else open

        with opener(input.ref, "rt") as fin, open(output.ref, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    old = line[1:].split()[0]

                    if old not in rename:
                        raise ValueError(
                            f"No sequence-report mapping found for: {old}"
                        )

                    fout.write(f">{rename[old]}\n")
                else:
                    fout.write(line)


rule bwa_index:
    input:
        ref = rules.rename_reference.output.ref
    output:
        amb  = REF_FASTA + ".amb",
        ann  = REF_FASTA + ".ann",
        bwt  = REF_FASTA + ".bwt",
        pac  = REF_FASTA + ".pac",
        sa   = REF_FASTA + ".sa",
        fai  = REF_FASTA + ".fai",
        dict = REF_FASTA.replace(".fa", "").replace(".fasta", "") + ".dict"
    threads: RES.get("bwa_index", {}).get("threads", 8)
    resources:
        mem_mb = RES.get("bwa_index", {}).get("mem_mb", 8000),
        walltime = RES.get("bwa_index", {}).get("walltime", 4)
    params:
        name = REF_NAME,
        species = config["ref"].get("species", "unknown")
    shell:
        r"""
        samtools faidx {input.ref}

        samtools dict \
          -a {params.name} \
          -s "{params.species}" \
          {input.ref} \
          -o {output.dict}

        bwa index {input.ref}
        """


###############################################################################
# Alignment
###############################################################################

rule align_reads:
    input:
        r1 = read1,
        r2 = read2,
        ref = rules.rename_reference.output.ref,
        amb = rules.bwa_index.output.amb,
        ann = rules.bwa_index.output.ann,
        bwt = rules.bwa_index.output.bwt,
        pac = rules.bwa_index.output.pac,
        sa  = rules.bwa_index.output.sa,
        fai = rules.bwa_index.output.fai
    output:
        bam = f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam",
        bai = f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam.bai"
    threads: RES.get("align_reads", {}).get("threads", 8)
    resources:
        mem_mb = RES.get("align_reads", {}).get("mem_mb", 18000),
        walltime = RES.get("align_reads", {}).get("walltime", 48)
    params:
        tmp = lambda wc: f"{ALIGN_ROOT}/{wc.sample}",
        rg = lambda wc: f"@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA",
        samtools_threads = lambda wc, threads: max(1, threads - 1)
    shell:
        r"""
        bwa mem -t {threads} \
          -R "{params.rg}" \
          {input.ref} \
          <(gzip -dc {input.r1}) \
          <(gzip -dc {input.r2}) \
        | samtools fixmate -@ {params.samtools_threads} -m - - \
        | samtools sort -@ {params.samtools_threads} -T {params.tmp}.sort.tmp - \
        | samtools markdup -@ {params.samtools_threads} -T {params.tmp}.markdup.tmp - {output.bam}

        samtools index -@ {params.samtools_threads} {output.bam}
        """


rule flagstat:
    input:
        bam = f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam"
    output:
        flagstat = f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam.flagstat"
    threads: RES.get("flagstat", {}).get("threads", 1)
    resources:
        mem_mb = RES.get("flagstat", {}).get("mem_mb", 4000),
        walltime = RES.get("flagstat", {}).get("walltime", 1)
    shell:
        r"""
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        """


###############################################################################
# Alignment target
###############################################################################

rule alignment:
    input:
        expand(
            f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam.bai",
            sample=SAMPLES,
        ),
        expand(
            f"{ALIGN_ROOT}/{{sample}}.fixmate.sort.markdup.rg.bam.flagstat",
            sample=SAMPLES,
        )
