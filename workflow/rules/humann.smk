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


pabun_cpm = f"humann/{config['sample']}_humann3_pathabundance_cpm.tsv"
ko_cpm = f"humann/{config['sample']}_humann3_KO_cpm.tsv"

all_inputs = [
    ko_cpm,
    pabun_cpm,
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


rule humann3_run_uniref90:
    input:
        fastq="kneaddata/{sample}_knead_cat.fastq.gz",
        metaphlan_profile=config["metaphlan_ab"],
        choco_db=config["choco_db"],
        uniref90_db=config["uniref90_db"],
    output:
        ab="humann/{sample}_humann3_pathabundance.tsv",
        genefam="humann/{sample}_humann3_genefamilies.tsv",
        stats="humann/{sample}_humann3_stats.txt",
    params:
        out_prefix="{sample}_humann3",
        out_dir=lambda w, output: os.path.dirname(output[0]),
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * 1024
        * max(input.fastq.size // 1000000000, 1)
        * 5 * attempt,
        runtime=lambda wc, attempt: 8 * 60 * attempt,
    threads: 64
    # we have an extra log in case there is an error with humann.  Cause
    # we skip the built in logging because
    # they dont actually log errors to their --o-log :(
    log:
        e="logs/humann_{sample}.e",
        o="logs/humann_{sample}.o",
    shell:
        """
        # see https://forum.biobakery.org/t/metaphlan-v4-0-2-and-huma-3-6-metaphlan-taxonomic-profile-provided-was-not-generated-with-the-expected-database/4296/8
        cat {input.metaphlan_profile} | cut -f 1-4 > {wildcards.sample}_tmp_metaphlan.tsv
        humann \
            --input {input.fastq} \
            --output {params.out_dir} \
            --output-basename {params.out_prefix} \
            --search-mode uniref90 \
            --remove-column-description-output \
            --protein-database {input.uniref90_db} \
            --nucleotide-database {input.choco_db} \
            --taxonomic-profile {wildcards.sample}_tmp_metaphlan.tsv  \
            --remove-temp-output \
            --output-max-decimals 5 \
            --threads {threads} \
            > {log.o} 2> {log.e}

        cat {log.o} | \
            grep "reads" | \
            grep "\: [0-9]" \
            > {output.stats}
        """


# Renormalize gene family and pathway abundance from RPK to relative abundance(CPM)
# for input of lefse
rule renormalize_pabun:
    input:
        path="humann/{sample}_humann3_pathabundance.tsv",
    output:
        path_out="humann/{sample}_humann3_pathabundance_cpm.tsv",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    log:
        e="logs/renomalize_pabun_{sample}.e",
    shell:
        """
        humann_renorm_table \
            --input {input.path} \
            --units cpm \
            -s n \
            --output {output.path_out} \
            2> {log.e}
        """


rule regroup_genefam_2_KO:
    input:
        genefamilies="humann/{sample}_humann3_genefamilies.tsv",
        utility_mapping_db=os.path.join(
            config["utility_mapping_db"], "map_ko_uniref90.txt.gz"
        ),
    output:
        out="humann/{sample}_humann3_KO.tsv",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    log:
        e="logs/regroup_genefam_2_ko_{sample}.e",
    shell:
        """
        # since humann_config doesn't work with docker we provide the uniref to ko as a
        # custom regrouping file
        # Dispite the file is  provided as ko's to uniref90 we DON'T add the
        # --reversed flag to convert our uniref90 genefamilies to KOs :(
        humann_regroup_table \
            --input {input.genefamilies} \
            --custom {input.utility_mapping_db} \
            --output {output.out} \
            2> {log.e}
        """


rule renormalize_KO:
    input:
        path="humann/{sample}_humann3_KO.tsv",
    output:
        path_out="humann/{sample}_humann3_KO_cpm.tsv",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    log:
        e="logs/renormalize_ko_{sample}.e",
    shell:
        """
        humann_renorm_table \
            --input {input.path} \
            --units cpm \
            -s n \
            --output {output.path_out} \
            2> {log.e}
        """


# Join humann output per sample into one table
rule join_table:
    input:
        res_dir="humann",
    output:
        pathabun="humann/humann3_pathabundance_cpm_joined_{sample}.tsv",
        ko_cpm="humann3_KO_cpm_joined_{sample}.tsv",
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    log:
        e="logs/join_table_{sample}.e",
    shell:
        """
        humann_join_tables \
            -s \
            --input {input.res_dir} \
            --file_name humann3_pathabundance_cpm \
            --output {output.pathabun} \
            2> {log.e}
        humann_join_tables \
            -s \
            --input {input.res_dir} \
            --file_name humann3_KO_cpm \
            --output {output.ko_cpm} \
            2>> {log.e}
        """


rule split_stratified:
    input:
        pabun="humann3_pathabundance_cpm_joined_{sample}.tsv",
        ko_cpm="humann3_KO_cpm_joined_{sample}.tsv",
    output:
        touch("split_stratified_{sample}.done"),
    params:
        output_dir=directory("humann3_final_out"),
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    log:
        e="logs/split_stratified_{sample}.e",
    shell:
        """
        humann_split_stratified_table \
            --input {input.pabun} \
            --output {params.output_dir} \
            2> {log.e}
        humann_split_stratified_table \
            --input {input.ko_cpm} \
            --output {params.output_dir} \
            2>> {log.e}
        """


