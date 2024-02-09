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


localrules:
    all,


quast_outputs = [
    "quast/quast_{sample}/report.pdf".format(sample=config["sample"]),
    "quast/quast_{sample}/transposed_report.tsv".format(sample=config["sample"]),
]
all_inputs = [quast_outputs]


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")
    if config["assembler"].lower() == "spades":
        if len(config["R1"]) > 1:
            print("WARNING: Concatenating multiple inputs into a single paired library as metaspades does not support multiple libraries")
        print("Running SPAdes")
    else:
        print("Running megahit")


if config["assembler"].lower() == "spades":
    assemblies = [f"spades_{config['sample']}.assembly.fasta", f"spades_{config['sample']}_metaviral/scaffolds.fasta"]
    assemblies_labels = ",".join([f"metaspades_{config['sample']}",f"metaviralspades_{config['sample']}"])
    all_inputs.append(f"{config['sample']}.cleaned_assembly_files")
    all_inputs.extend(assemblies)
# elif config["assembler"].lower() == "both": # should this be enabled?
#     assemblies = [f"spades_{config['sample']}.assembly.fasta", f"spades_{config['sample']}_metaviral/scaffolds.fasta", f"megahit_{config['sample']}.assembly.fasta"]
#     assemblies_labels = ",".join([f"metaspades_{config['sample']}",f"metaviralspades_{config['sample']}",f"megahit_{config['sample']}"])
#     all_inputs.append(f"{config['sample']}.cleaned_assembly_files")
#     all_inputs.extend(assemblies)
else:
    assemblies = [f"megahit_{config['sample']}.assembly.fasta"]
    assemblies_labels = f"megahit_{config['sample']}"
    all_inputs.extend(assemblies)


rule all:
    input:
        all_inputs,


#
if len(config["R1"]) == 1:
    input_R1 = config["R1"]
    input_R2 = config["R2"]
else:
    input_R1 = [f"concatenated/{config['sample']}_R1.fastq.gz"]
    input_R2 = [f"concatenated/{config['sample']}_R2.fastq.gz"]


# Utils Module
module utils:
    snakefile:
        "utils.smk"
    config:
        config
    skip_validation:
        True


use rule concat_R1_R2 from utils as utils_concat_R1_R2 with:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp("concatenated/{sample}_R1.fastq.gz"),
        R2=temp("concatenated/{sample}_R2.fastq.gz"),
    log:
        e="logs/concat_r1_r2_{sample}.e",

rule megahit:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        outdir=directory("megahit_{sample}"),
        assembly="megahit_{sample}.assembly.fasta",
    container:
        config["docker_megahit"]
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 20),
        runtime=24 * 60,
    threads: 64
    params:
        input_string=lambda wildcards, input: str(
            "-1 " + ",".join(input.R1) + " -2 " + ",".join(input.R2)
        ),
    shell:
        """
    mkdir -p ${{TMPDIR}}/megahit_{wildcards.sample}/
    megahit {params.input_string} --out-dir megahit_{wildcards.sample}/ --out-prefix {wildcards.sample} --tmp-dir ${{TMPDIR}}/megahit_{wildcards.sample}/ --memory $(({resources.mem_mb} * 1024 ))  --num-cpu-threads {threads}
    rm -r ${{TMPDIR}}/megahit_{wildcards.sample}/
    mv megahit_{wildcards.sample}/{wildcards.sample}.contigs.fa {output.assembly}
    """


rule SPAdes_run:
    # TODO: add in params for read length to experiment with larger kmers than default
    input:
        R1=input_R1,
        R2=input_R2,
    output:
        assembly="spades_{sample}.assembly.fasta",
        graph="spades_{sample}.assembly_graph.gfa",
        spades_log="spades_{sample}/spades.log",
    container:
        config["docker_spades"]
    conda:
        "../envs/spades.yaml"
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 20),
        runtime=48 * 60,
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
        mv spades_{wildcards.sample}/assembly_graph_with_scaffolds.gfa {output.graph}
        """

rule viral_SPAdes_run:
    # max_kmer https://github.com/ablab/spades/discussions/1188
    # --onlyassembler seems to be neccessary when using assembly graph input
    input:
        R1=input_R1,
        R2=input_R2,
        assembly_graph=rules.SPAdes_run.output.graph,
    output:
        assembly="spades_{sample}_metaviral/scaffolds.fasta",
        spades_log="spades_{sample}_metaviral/spades.log",
    container:
        config["docker_spades"]
    params:
        max_kmer=55
    conda:
        "../envs/spades.yaml"
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 20),
        runtime=48 * 60,
    threads: 64
    log:
        e="logs/spades_{sample}_metaviral.log",
    shell:
        """
        spades.py \
            --metaviral \
            -1 {input.R1} \
            -2 {input.R2} \
            --assembly-graph {input.assembly_graph} \
            --only-assembler \
            -t {threads} \
            -o spades_{wildcards.sample}_metaviral/ \
            -m $(({resources.mem_mb}/1024)) \
            -k {params.max_kmer}  \
            2> {log.e}
        # deal with missing scaffolds file if no viruses recovered
        if grep -q "No complete extrachromosomal contigs assembled" "spades_{wildcards.sample}_metaviral/spades.log"; then
            touch {output.assembly}
        fi
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
        assembly=assemblies,
        blast_16s_db_nsq=config["blast_16s_db_nsq"],
    output:
        report_tsv="quast/quast_{sample}/transposed_report.tsv",
        report="quast/quast_{sample}/report.pdf",
    params:
        dir="quast/quast_{sample}/",
        labels=assemblies_labels,
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
        runtime=3 * 60,
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
          --labels {params.labels} \
          > {log.o} 2> {log.e}
        """


rule clean_up:
    """"{sample}_metaerg.gff" is used as an input to ensure
    that step is done before we clean.
    """
    input:
        agg_files=[
            "quast/quast_{sample}/report.pdf".format(sample=config["sample"]),
        ],
        spades_log="spades_{sample}/spades.log",
    output:
        touch("{sample}.cleaned_assembly_files"),
    shell:
        """
        # We dont need any of the spades intermediate files (corrected reads, temp, per-kmer assembly)
        ls spades_{wildcards.sample}/
        find spades_{wildcards.sample}/  -type f | xargs --no-run-if-empty rm
        """
