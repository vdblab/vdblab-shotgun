import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


validate(config, os.path.join(str(workflow.basedir), "../../config/config.schema.yaml"))


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,

metaphlan = f"metaphlan/{config['sample']}_metaphlan3_profile.txt"
metaphlan_sam = f"metaphlan/{config['sample']}.sam.bz2"
krona = f"reports/{config['sample']}_metaphlan3_profile.txt.krona.html"

all_inputs = [
    metaphlan,
    krona,
    metaphlan_sam,
]


rule all:
    input:
        all_inputs,


# the cat paired end reads and metaphlan and humann3 part
rule cat_pair:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        joined=temp("kneaddata/{sample}_knead_cat.fastq.gz"),
    conda:
        "../envs/base.yaml"
    log:
        e="logs/cat_pair_{sample}.e",
    shell:
        "cat {input.R1} {input.R2} > {output.joined} 2> {log.e}"

# Metaphlan3 and Strainphlan
rule metaphlan_run:
    input:
        fastq="kneaddata/{sample}_knead_cat.fastq.gz",
        db=config["metaphlan_db"],
    output:
        outfile="metaphlan/{sample}_metaphlan3_profile.txt",
        sam="metaphlan/{sample}.sam.bz2",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/metaphlan.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: 30 * 1024 * attempt,
        runtime=lambda wc, attempt: 2 * 60 * attempt,
    threads: 64
    log:
        e="logs/metaphlan_{sample}.e",
    shell:
        """
        # the presense of this file causes an error from metaphlan
        # which makes rerunning irritating
        if [ -f "{input.fastq}.bowtie2out.txt" ]
        then
            rm {input.fastq}.bowtie2out.txt
        fi
        export METAPHLAN_BOWTIE2_DB={input.db}
        metaphlan {input.fastq} \
            --bowtie2db {input.db} \
            --index mpa_vJan21_CHOCOPhlAnSGB_202103 \
            --input_type fastq \
            --sample_id  {wildcards.sample} \
            -s {output.sam} \
            --add_viruses \
            --unclassified_estimation \
            --nproc {threads} \
            -t rel_ab_w_read_stats \
            -o {output.outfile} \
            2> {log.e}
        """


rule metaphlan2_krona:
    input:
        infile="metaphlan/{sample}_metaphlan3_profile.txt",
    output:
        outfile="metaphlan/{sample}_metaphlan3_profile.txt.krona",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/metaphlan.yaml"
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/metaphlan2krona_{sample}.e",
    shell:
        """
        cat {input.infile} | cut -f 1-4 > {wildcards.sample}_tmp_metaphlan.tsv
        metaphlan2krona.py -p {wildcards.sample}_tmp_metaphlan.tsv -k {output.outfile} 2> {log.e}
        """


rule krona:
    input:
        infile="metaphlan/{sample}_metaphlan3_profile.txt.krona",
    output:
        outfile="reports/{sample}_metaphlan3_profile.txt.krona.html",
    container:
        config["docker_krona"]
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/krona_{sample}.e",
    shell:
        """
        ktImportText {input.infile} -o {output.outfile}
        """
