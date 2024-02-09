import os
import glob
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


# the str is needed as when running from github workflow.current_basedir is a Githubfile, not a string or a path so os.path objects
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")



envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])
MAGs = glob.glob(os.path.join(config["bindir"] + "*.fa"))
if not MAGs:
    raise ValueError(f'No MAGs found with pattern {os.path.join(config["bindir"] + "*.fa")}')

onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,




rule all:
    input:
        f"gtdbtk-{config['sample']}/gtdbtk.json",
        f"kraken2/{config['sample']}_kraken2.report"

rule gtdb_classify:
    """ Cant have the tabular output because if a sample lacks bacteria or archaea it will omit that domains output file
    """
    input:
        bindir=config["bindir"],
    output:
        "gtdbtk-{sample}/gtdbtk.json"
    container:
        config["docker_gtdbtk"]
    params:
        db=config["gtdb_db"],
        outdir=lambda wc, output: os.path.dirname(output[0])
    threads: 32
    resources:
        # ~64 is recommended by their docs
        mem_mb=68*1024
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf --genome_dir {input.bindir} --cpus {threads} -x fa --out_dir {params.outdir} --mash_db {params.db}/mash
        """

module kraken:
    snakefile:
        "kraken.smk"
    config:
        config
    skip_validation:
        True

use rule kraken_standard_run from kraken with:
    input:
        R1=MAGs,
        db=config["kraken2_db"],
    output:
        out="kraken2/{sample}_kraken2.out",
        unclass_R1=temp("kraken2/{sample}_kraken2_unclassified_1.fastq"),
        unclass_R2=temp("kraken2/{sample}_kraken2_unclassified_2.fastq"),
        report="kraken2/{sample}_kraken2.report"
