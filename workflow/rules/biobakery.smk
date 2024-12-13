import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


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
metaphlan = f"metaphlan/{config['sample']}_metaphlan3_profile.txt"
ko_cpm = f"humann/{config['sample']}_humann3_KO_cpm.tsv"
metaphlan_sam = f"metaphlan/{config['sample']}.sam.bz2"
krona = f"reports/{config['sample']}_metaphlan3_profile.txt.krona.html"

all_inputs = [
    ko_cpm,
    metaphlan,
    pabun_cpm,
    krona,
    metaphlan_sam,
]


rule all:
    input:
        all_inputs,


# the cat paired end reads and metaphlan and humann3 part
rule cat_pair:
    input:
        unpack(get_config_inputs),
    output:
        joined=temp("kneaddata/{sample}_knead_cat.fastq.gz"),
    conda:
        "../envs/base.yaml"
    log:
        e="logs/cat_pair_{sample}.e",
    shell:
        "cat {input} > {output.joined} 2> {log.e}"


rule humann3_stage1_nucleotide:
    input:
        fastq="kneaddata/{sample}_knead_cat.fastq.gz",
        metaphlan_profile="metaphlan/{sample}_metaphlan3_profile.txt",
        choco_db=config["choco_db"],
    output:
        # in the case where there is no Chocodb, we only generate the log:
        log="humann/{sample}_humann3_nucleotide_align_humann_temp/{sample}_humann3_nucleotide_align.log"
        #aligned_tsv="humann/{sample}_humann3_nucleotide_align_humann_temp/{sample}_humann3_nucleotide_align_bowtie2_aligned.tsv",
        #unaligned_fa="humann/{sample}_humann3_nucleotide_align_humann_temp/{sample}_humann3_nucleotide_align_bowtie2_unaligned.fa",
    params:
        out_prefix="{sample}_humann3_nucleotide_align",
        out_dir=lambda w, output: os.path.dirname(os.path.dirname(output[0])),
    container:
        config["docker_biobakery"]
    conda:
        "../envs/humann.yaml"
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * 1024
        * max(input.fastq.size // 1000000000, 1)
        * 5,
        runtime=lambda wc, attempt: 8 * 60 * attempt,
    threads: 30
    # we have an extra log in case there is an error with humann.  Cause
    # we skip the built in logging because
    # they dont actually log errors to their --o-log :(
    log:
        e="logs/humann_nucleotide_align_{sample}.e",
        o="logs/humann_nucleotide_align_{sample}.o",
    shell:
        """
        # see https://forum.biobakery.org/t/metaphlan-v4-0-2-and-huma-3-6-metaphlan-taxonomic-profile-provided-was-not-generated-with-the-expected-database/4296/8
        cat {input.metaphlan_profile} | cut -f 1-4 > {wildcards.sample}_tmp_metaphlan.tsv
        humann \
            --input {input.fastq} \
            --output {params.out_dir} \
            --output-basename {params.out_prefix} \
            --bypass-translated-search \
            --remove-column-description-output \
            --nucleotide-database {input.choco_db} \
            --taxonomic-profile {wildcards.sample}_tmp_metaphlan.tsv  \
            --output-max-decimals 5 \
            --threads {threads} \
            > {log.o} 2> {log.e}
        """

rule humann3_stage2_shard_unaligned_nucleotide_reads:
    input:
        log="humann/{sample}_humann3_nucleotide_align_humann_temp/{sample}_humann3_nucleotide_align.log",
    output:
        unaligned_shards=temp(
            expand(
                "humann/split_unaligned_fasta/{{sample}}_humann3_nucleotide_align_bowtie2_unaligned.part_{shard}.fa",
                shard=SHARDS,
            )
        ),
    params:
        num_shards=len(SHARDS),
        out_dir=lambda w, output: os.path.dirname(output[0]),
    container:
        config["docker_seqkit"]
    conda:
        "../envs/humann.yaml"
    resources:
        mem_mb=lambda wc, attempt: attempt * 1024 * 5,
        runtime=lambda wc, attempt: 60 * attempt,
    threads: 4  # see their docs (apparently!)
    log:
        e="logs/humann_shard_{sample}.e",
        o="logs/humann_shard_{sample}.o",
    shell:
        """
        NUC_ALIGN_OUT=humann/{wildcards.sample}_humann3_nucleotide_align_humann_temp/{wildcards.sample}_humann3_nucleotide_align_bowtie2_unaligned.fa
        ORIGINAL_FASTQ=kneaddata/{wildcards.sample}_knead_cat.fastq.gz
        if [ ! -f $NUC_ALIGN_OUT ]; then zcat $ORIGINAL_FASTQ > $NUC_ALIGN_OUT; fi
        seqkit split2 \
            --threads {threads} \
            --read1 $NUC_ALIGN_OUT \
            --by-part {params.num_shards} \
            --force \
            --out-dir {params.out_dir}
        """

rule humann3_stage3_translated_alignment_of_shards:
    input:
        unaligned_shard = "humann/split_unaligned_fasta/{sample}_humann3_nucleotide_align_bowtie2_unaligned.part_{shard}.fa",
    output:
        translated_aligned_reads =temp("humann/translated_aligned_reads/{sample}_translated_aligned.part_{shard}.tsv"),
    params:
        diamond_db=config["diamond_db_file"],
        evalue_threshold=1.0, #from humann defaults
        outdir=lambda w, output: os.path.dirname(output[0]),
    container:
        config["docker_diamond"]
    resources:
        mem_mb=lambda wc, attempt: 1024 * 5 * attempt,
        runtime=lambda wc, attempt: 5 * 60 * attempt,
    threads: 4
    log:
        e="logs/humann_translation_align_diamond_{sample}_{shard}.e",
        o="logs/humann_translation_align_diamond_{sample}_{shard}.o",
    shell:
        """
        diamond blastx \
            --query {input.unaligned_shard} \
            --evalue {params.evalue_threshold} \
            --threads {threads} \
            --db {params.diamond_db} \
            --out {output.translated_aligned_reads} \
            --tmpdir {params.outdir}"
        """

rule humann3_stage4_combine_all_tsvs:
    input:
        log="humann/{sample}_humann3_nucleotide_align_humann_temp/{sample}_humann3_nucleotide_align.log",
        translation_aligned_shards=expand("humann/translated_aligned_reads/{{sample}}_translated_aligned.part_{shard}.tsv",
                shard=SHARDS),
    output:
        all_aligned_tsv="humann/{sample}_humann3_all_aligned.tsv",
    resources:
        mem_mb=lambda wc, attempt: 1024 * 5 * attempt,
        runtime=lambda wc, attempt: 5 * 60 * attempt,
    threads: 4
    log:
        e="logs/humann_combine_shards_{sample}.e",
        o="logs/humann_combine_shards_{sample}.o",
    shell:
        """
        NUC_ALIGN_OUT=humann/{sample}_humann3_nucleotide_align_temp/{sample}_humann3_nucleotide_align_bowtie2_aligned.tsv
        if [ -f $NUC_ALIGN_OUT ]; then cat $NUC_ALIGN_OUT > {output.all_aligned_tsv}; fi
        cat {input.translated_aligned_reads} >> {output.all_aligned_tsv}
        """

rule humann3_stage5_rerun_humann_on_all_alignments:
    input:
        all_aligned_tsv="humann/{sample}_humann3_all_aligned.tsv",
    output:
        ab="humann/{sample}_humann3_pathabundance.tsv",
        genefam="humann/{sample}_humann3_genefamilies.tsv",
        stats="humann/{sample}_humann3_stats.txt",
    params:
        out_prefix="{sample}_humann3",
        out_dir=lambda w, output: os.path.dirname(output[0]),
    resources:
        mem_mb=lambda wc, attempt: 1024 * 5 * attempt,
        runtime=lambda wc, attempt: 5 * 60 * attempt,
    threads: 4
    log:
        e="logs/humann_combine_shards_{sample}.e",
        o="logs/humann_combine_shards_{sample}.o",
    shell:
        """
        humann \
            --input {input.all_aligned_tsv} \
            --input-format blastm8 \
            --output {params.out_dir} \
            --output-basename {params.out_prefix} \
            --bypass-nucleotide-search \
            --bypass-translated-search \
            --remove-column-description-output \
            --output-max-decimals 5 \
            --threads {threads}
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
        # first submission is given 30GB, then 45,
        mem_mb=lambda wildcards, attempt: 30 * 1024 * attempt,
        runtime=lambda wc, attempt: 2 * 60 * attempt,
    threads: 32
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
