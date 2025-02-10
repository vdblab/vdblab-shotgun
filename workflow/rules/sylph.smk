import glob
import os
import sys
from collections import defaultdict


include: "common.smk"


sylphdbpath = config["sylphdbpath"]

addn_profiles = []
if config["addn_dbs"]:
    addn_dbs = {}
    for _db in  config["addn_dbs"]:
        dbname = os.path.basename(_db).replace(".syldb", "")
        metadata_path = glob.glob(_db.replace("syldb", "metadata*"))
        if not metadata_path or len(metadata_path) > 1:
            raise ValueError("Additional databases must have corresponding metadata, eg path/to/db.syldb must had path/to/db.metadata.tsv")
        addn_dbs[dbname] = {
            "db": _db,
            "metadata": metadata_path[0]
            }
    addn_profiles = expand(f"addn_taxprofiles/{{db}}/{config['sample']}.sylphmpa", db=addn_dbs.keys())

dbs = {
    "fungi": {
        "db": os.path.join(sylphdbpath, "fungi-refseq-2024-07-25-c200-v0.3.syldb"),
        "metadata": os.path.join(
            sylphdbpath, "fungi_refseq_2024-07-25_metadata.tsv.gz"
        ),
    },
    "viruses": {
        "db": os.path.join(sylphdbpath, "imgvr_c200_v0.3.0.syldb"),
        "metadata": os.path.join(sylphdbpath, "IMGVR_4.1_metadata.tsv.gz"),
    },
    "prok": {
        "db": os.path.join(sylphdbpath, "gtdb-r220-c200-dbv1.syldb"),
        "metadata": os.path.join(sylphdbpath, "gtdb_r220_metadata.tsv.gz"),
    },
}

merged_anno_db = os.path.join(config["sylphdbpath"], "3kingdom_metadata.tsv")
if "R2" not in config:
    raise ValueError("This is not configured to work on single-end data")

rule all:
    input:
        f"taxprofiles/{config['sample']}.sylphmpa",
        addn_profiles,


def get_sylph_input(wc):
    if len(config["R1"]) == 1:
        input_R1 = config["R1"]
        input_R2 = config["R2"]
    else:
        input_R1 = f"concatenated/{config['sample']}_R1.fastq.gz"
        input_R2 = f"concatenated/{config['sample']}_R2.fastq.gz"
    return {"R1": input_R1, "R2": input_R2}


# Utils Module
module utils:
    snakefile:
        "utils.smk"
    config:
        config
    skip_validation:
        True


use rule concat_lanes_fix_names from utils as utils_sylph_concat_lanes_fix_names with:
    input:
        fq=get_sylph_input,
    output:
        fq=temp("concatenated/{sample}_R{rd}.fastq.gz"),
    log:
        e="logs/concat_lanes_fix_names_{sample}_R{rd}.e",


rule sketch:
    input:
        unpack(get_sylph_input),
    output:
        sketch="sketches/{sampleid}.paired.sylsp",
    params:
        inpstr=lambda wc, input: (
            f"-1 {input.R1} -2 {input.R2}"
            if hasattr(input, "R2")
            else f"-1 {input.R1}"
        ),
        outdir=lambda wc, output: os.path.dirname(output.sketch),
    threads: 3
    container:
        "docker://ghcr.io/vdblab/sylph:0.6.1a"
    shell:
        "sylph sketch {params.inpstr} --sample-names {wildcards.sampleid} --sample-output-directory {params.outdir} -t {threads}"


rule profile:
    input:
        sketch="sketches/{sampleid}.paired.sylsp",
        db=[y["db"] for x, y in dbs.items()],
    output:
        profile="profiles/{sampleid}.tsv",
    threads: 3
    container:
        "docker://ghcr.io/vdblab/sylph:0.6.1a"
    resources:
        mem_mb=lambda wc, attempt: 1024 * attempt * 22,
    shell:
        "sylph profile -t {threads}  {input.db} {input.sketch} -o {output.profile}"


rule create_taxa_profile:
    input:
        profile="profiles/{sampleid}.tsv",
        dbmeta=merged_anno_db,
    output:
        profile="taxprofiles/{sampleid}.sylphmpa",
    resources:
        mem_mb=lambda wc, attempt: 1024 * attempt * 8,
    params:
        prefix=lambda wc, output: os.path.dirname(output.profile) + "/",
    container:
        "docker://ghcr.io/vdblab/sylph:0.6.1a"
    shell:
        "sylph_to_taxprof.py -s {input.profile} -m  {input.dbmeta}  -o {params.prefix}"


def get_addn_db(wildcards):
    return(addn_dbs[wildcards.db]["db"])
def get_addn_metadata(wildcards):
    return(addn_dbs[wildcards.db]["metadata"])


use rule profile as profile_addn with:
    input:
        sketch="sketches/{sampleid}.paired.sylsp",
        db=get_addn_db,
    output:
        profile="addn_db_profiles/{db}/{sampleid}.tsv",

use rule create_taxa_profile as create_taxa_profile_addn with:
    input:
        profile="addn_db_profiles/{db}/{sampleid}.tsv",
        dbmeta=get_addn_metadata,
    output:
        profile="addn_taxprofiles/{db}/{sampleid}.sylphmpa",
