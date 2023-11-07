import os
import glob
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


# the str is needed as when running from github workflow.current_basedir is a Githubfile, not a string or a path so os.path objects
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")


validate(
    config,
    os.path.join(str(workflow.current_basedir), "../../config/config.schema.yaml"),
)


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
        "gtdb-tk/results.txt",
        f"kraken2/{config['sample']}_kraken2.report"

rule gtdb_classify:
    input:
        bindir=config["bindir"],
    output:
        "gtdb-tk/results.txt"
    container:
        config["docker_gtdbtk"]
    params:
        db="/data/brinkvd/resources/dbs/GTDB/214/release214/"
    threads: 16
    resources:
        # ~64 is recommended by their docs
        mem_mb=68*1024
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf --genome_dir {input.bindir} --cpus {threads} -x fa --out_dir results/ --mash_db {params.db}/mash
        touch {output}
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
