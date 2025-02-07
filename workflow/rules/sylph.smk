import os
import sys
from collections import defaultdict

include: "common.smk"


sylphdbpath=config["sylphdbpath"]


dbs = {
    "fungi": {
        "db": os.path.join(sylphdbpath, "fungi-refseq-2024-07-25-c200-v0.3.syldb"),
        "metadata": os.path.join(sylphdbpath, "fungi_refseq_2024-07-25_metadata.tsv.gz")},
    "viruses": {
        "db": os.path.join(sylphdbpath, "imgvr_c200_v0.3.0.syldb"),
        "metadata": os.path.join(sylphdbpath, "IMGVR_4.1_metadata.tsv.gz")},
    "prok": {
        "db": os.path.join(sylphdbpath, "gtdb-r220-c200-dbv1.syldb"),
        "metadata": os.path.join(sylphdbpath, "gtdb_r220_metadata.tsv.gz")}
    }

merged_anno_db = os.path.join(config["sylphdbpath"], "3kingdom_metadata.tsv")
if "R2" not in config:
    raise ValueError("This is not configured to work on single-end data")

rule all:
    input:
        f"taxprofiles/{config['sample']}.sylphmpa",


# Utils Module
def get_sylph_input(wc):
    if len(config["R1"]) == 1:
        input_R1 = config["R1"]
        input_R2 = config["R2"]
    else:
        input_R1 = f"concatenated/{config['sample']}_R1.fastq.gz"
        input_R2 = f"concatenated/{config['sample']}_R2.fastq.gz"
    return {"R1": input_R1, "R2": input_R2}


module utils:
    snakefile:
        "utils.smk"
    config:
        config
    skip_validation:
        True


use rule concat_lanes_fix_names from utils as utils_sylph_concat_lanes_fix_names with:
    input:
        fq=get_concat_input,
    output:
        fq=temp("concatenated/{sample}_R{rd}.fastq.gz"),
    log:
        e="logs/concat_lanes_fix_names_{sample}_R{rd}.e",

rule sketch:
    input:
        unpack(get_sylph_input),
    output:
        sketch="sketches/{sampleid}.paired.sylsp"
    params:
        inpstr=lambda wc, input: (
            f"-1 {input.R1} -2 {input.R2}" if hasattr(input, "R2") else f"-1 {input.R1}"
        ),
        outdir=lambda wc, output: os.path.dirname(output.sketch),

    threads: 3
    container: "docker://ghcr.io/vdblab/sylph:0.6.1a"
    shell:"sylph sketch {params.inpstr} --sample-names {wildcards.sampleid} --sample-output-directory {params.outdir} -t {threads}"

rule profile:
    input:
        sketch="sketches/{sampleid}.paired.sylsp",
        db=[y["db"] for x, y in dbs.items()],
    output:
        profile="profiles/{sampleid}.tsv"
    threads: 3
    container: "docker://ghcr.io/vdblab/sylph:0.6.1a"
    resources:
        mem_mb = lambda wc, attempt: 1024 * attempt * 22,
    shell:"sylph profile {input.db} {input.sketch} -o {output.profile}"


# rule merge_anno:
#     input:
#         db=[y["metadata"] for x, y in dbs.items()],
#     output:
#         "3kingdommeta.tsv.gz"
#     shell:
#         "cat {input.db} > {output}"

rule create_taxa_profile:
    input:
        profile="profiles/{sampleid}.tsv",
        dbmeta=merged_anno_db,
    output:
        profile="taxprofiles/{sampleid}.sylphmpa"
    resources:
        mem_mb = lambda wc, attempt: 1024 * attempt * 8,
    params:
        prefix = lambda wc, output: os.path.dirname(output.profile) + "/",
    container: "docker://ghcr.io/vdblab/sylph:0.6.1a"
    shell: "sylph_to_taxprof.py -s {input.profile} -m  {input.dbmeta}  -o {params.prefix}"

# rule micom_agg:
#     input:
#         profile=expand("taxprofiles/{sampleid}.sylphmpa", sampleid = [f"{x}" for x, y in inputs.items()]),
#     output:
#         "micom_species.csv"
#     shell:""" echo -e "sample_id,genus,species,abundance" > {output}
#     for f in {input.profile}
#     do
#        x=$( basename $f | sed "s|.sylphmpa||g")
#        cat $f | grep "s__" | grep -v "t__" | sed "s|^.*g__|g__|g" | sed "s | \t g" | cut -f 1-3 | sed "s|^|$x\t|g" | tr '	' ,
#     done >> {output}
# """
