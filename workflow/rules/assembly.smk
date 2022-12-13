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
    with open("config_used.yml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


assemblies = f"{config['sample']}.assembly.fasta"
quast_outputs = [
    "quast/quast_{sample}/report.pdf".format(sample=config["sample"]),
    "quast/quast_{sample}/transposed_report.tsv".format(sample=config["sample"]),
]


rule all:
    input:
        assemblies,
        quast_outputs,


rule SPAdes_run:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        assembly="{sample}.assembly.fasta",
    container:
        config["docker_spades"]
    conda:
        "../envs/spades.yaml"
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 20),
        runtime="48:00",
    threads: 64
    log:
        e="logs/spades_{sample}.log",
    shell:
        """
        spades.py \
            -1 {input.R1} \
            -2 {input.R2} \
            -t {threads} \
            --meta \
            -o spades_{wildcards.sample} \
            -m $(({resources.mem_mb}/1024)) \
            2> {log.e}
        mv spades_{wildcards.sample}/scaffolds.fasta {output.assembly}
        """


rule quast_run:
    """ see http://quast.sourceforge.net/docs/manual.html
    NOTE: giving --blast-db {input.blast_16s_db_nsq} as the NCBI 16s db
    results in silent errors because their automatic genome retrieval
    uses names not accessions :(
    so instead we allow it to use its own SILVA one
    I tried downloading their silva and preprocessing it but
    ran into additional errors, presumably due to the OLD version
    of BLAST used by quast...
    TODO: make --reference work with chocophlan db somehow?

    As of 2022-11-09 we switched to non-meta quast for runtime issues
    with pulling genomes
    """
    input:
        assembly="{sample}.assembly.fasta",
        blast_16s_db_nsq=config["blast_16s_db_nsq"],
    output:
        report_tsv="quast/quast_{sample}/transposed_report.tsv",
        report="quast/quast_{sample}/report.pdf",
    params:
        dir="quast/quast_{sample}/",
    threads: 16
    container:
        config["docker_quast"]
    conda:
        "../envs/spades.yaml"
    log:
        e="logs/quast_{sample}.e",
        o="logs/quast_{sample}.o",
    resources:
        mem_mb=8 * 1024,
        runtime="3:00",
    shell:
        """
        quast.py \
          -o {params.dir} \
          --split-scaffolds \
          --threads {threads} \
          {input.assembly} \
          --ambiguity-usage all \
          --min-contig 100 \
          --no-snps \
          --no-icarus \
          --labels {wildcards.sample} \
          > {log.o} 2> {log.e}
        """
