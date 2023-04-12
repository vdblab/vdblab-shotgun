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


kraken_outputs = f"kraken2/{config['sample']}_kraken2.report"
kraken_unclassified_outputs = expand(
    "kraken2/{sample}_kraken2_unclassified_{readdir}.fastq.gz",
    sample=config["sample"],
    readdir=[1, 2],
)
brackenreport = expand(
    "kraken2/{sample}_kraken2.bracken.{taxlevel}.report",
    sample=config["sample"],
    taxlevel=["G", "S"],
)
krona = expand(
    "reports/{sample}_kraken2.bracken.S.report.krona.html", sample=config["sample"]
)


rule all:
    input:
        kraken_outputs,
        kraken_unclassified_outputs,
        brackenreport,

#
# Utils Module
if len(config["R1"]) == 1:
    input_R1 = config["R1"]
    input_R2 = config["R2"]
else:
    input_R1 = "concatenated/{sample}_R1.fastq.gz"
    input_R2 = "concatenated/{sample}_R2.fastq.gz"

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


rule kraken_standard_run:
    """ Profile the microbiome with Kraken2.

    We run this on the raw reads rather than the host-depleted trimmed reads because
    - its nice to have validation of the host detection percentages
    - we want uniform read lengths for bracken
    - its what the authors seem to recommend, although
      Nick has asked for confirmation: https://github.com/DerrickWood/kraken2/issues/646

    We set the confidence threshold after reading
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full

    We wanted to
      - minimize the L1 distance to the known community composition
      - maximize the alpha diversity similarity to the actual diversity

    We don't care as much about reads assigned so long as the assignements we
    do get are high quality and representative
    """
    input:
        R1=input_R1,
        R2=input_R2,
        db=config["kraken2_db"],
    output:
        out="kraken2/{sample}_kraken2.out",
        unclass_R1="kraken2/{sample}_kraken2_unclassified_1.fastq",
        unclass_R2="kraken2/{sample}_kraken2_unclassified_2.fastq",
        report="kraken2/{sample}_kraken2.report",
    params:
        unclass_template=lambda wildcards, output: output["unclass_R1"].replace(
            "_1.fastq", "#.fastq"
        ),
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/kraken_{sample}.log",
    resources:
        # use 8GB mem if using a minikraken db to make testing easier
        mem_mb=lambda wildcards, attempt: 8 * 1024
        if "mini" in config["kraken2_db"]
        else (64 * 1024 * attempt),
    threads: 16
    shell:
        """
        kraken2 \
            --threads {threads} \
            --use-names \
            --confidence 0.2 \
            --unclassified-out {params.unclass_template} \
            --db {input.db} \
            --report {output.report} \
            --paired {input.R1} {input.R2} \
            > {output.out} 2> {log.e}

        # some (mock) datasets are perfect but we still need these files
        if [ ! -f "{output.unclass_R1}" ]
        then
            touch {output.unclass_R1}
            touch {output.unclass_R2}
        fi
        """


checkpoint get_read_len:
    # get the read length needed by bracken, and name a file in its honor
    input:
        R1=config["R1"],
    output:
        readlen_dir=directory("kraken2/readlens/"),
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/get_read_len.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8 * 1024,
    threads: 1
    shell:
        """
        # thanks https://www.biostars.org/p/72433/ for the readlen calculation with awk
        # however, there is a strange error with awk here where it returns the proper result but returns an error code of 1.
        # we allow pipefails here to account for that
        set +o pipefail
        READLEN=$(zcat {input.R1} | head -n 100 | awk "NR%4 == 2 {{lengths[length(\$0)]++}} END {{for (l in lengths) {{print l}}}}")
        mkdir -p kraken2/readlens
        touch kraken2/readlens/len${{READLEN}}
        set -o pipefail
        """


rule bracken:
    """Run bracken using the probabilities based on the closest read length pre-computed.
    Note that this depreciates the need for the get_db_mers function,
    but we are leaving it here in case we decide that pre-computed
    probabilities are inadequate
    """
    input:
        #db_mers=get_dbs_needed,
        readlen_file=parse_read_lengths,
        report="kraken2/{sample}_kraken2.report",
        db=config["kraken2_db"],
    output:
        report="kraken2/{sample}_kraken2.bracken.{taxlevel}.report",
        out="kraken2/{sample}_kraken2.bracken.{taxlevel}.out",
    params:
        threshold=0,
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/bracken_{sample}.{taxlevel}.log",
    resources:
        mem_mb=8 * 1024,
    threads: 4
    shell:
        """
        # trim off the approximate and the len part of the read len signal
        READLEN=$(basename {input.readlen_file} | sed 's|len||' | sed  's|approx||')
        bracken -d {input.db} -i {input.report} -o {output.out} -w {output.report} -r $READLEN -l {wildcards.taxlevel} -t {params.threshold} 2> {log.e}
        find .
        """


rule kraken2krona:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        report="kraken2/{sample}_kraken2.bracken.S.report",
    output:
        report="kraken2/{sample}_kraken2.bracken.S.report.krona",
    conda:
        "../envs/kraken2.yaml"
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken2krona_{sample}.e",
        o="logs/kraken2krona_{sample}.o",
    shell:
        """
        kreport2krona.py --report-file {input.report} \
        --output {output.report} > {log.o} 2> {log.e}
        """


rule krona:
    input:
        report="kraken2/{sample}_kraken2.bracken.S.report.krona",
    output:
        report="reports/{sample}_kraken2.bracken.S.report.krona.html",
    container:
        config["docker_krona"]
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/krona_{sample}.e",
    shell:
        """
        ktImportText {input.report} -o {output.report}
        """


rule compress_kraken_unclassified:
    input:
        R1="kraken2/{sample}_kraken2_unclassified_1.fastq",
        R2="kraken2/{sample}_kraken2_unclassified_2.fastq",
    output:
        R1="kraken2/{sample}_kraken2_unclassified_1.fastq.gz",
        R2="kraken2/{sample}_kraken2_unclassified_2.fastq.gz",
    conda:
        "../envs/pigz.yaml"
    # this container is used in the 16S pipeline
    container:
        config["docker_cutadapt"]
    threads: 8
    log:
        e="logs/kraken_compressed_unclassified_{sample}.e",
    shell:
        """
        cat {input.R1} | pigz -p {threads} -9 > {output.R1} 2> {log.e}
        cat {input.R2} | pigz -p {threads} -9 > {output.R2} 2>> {log.e}
        """


rule kraken_merge_shards:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        reports=[
            f"kraken2/{config['sample']}_shard{shard}_kraken2.report"
            for shard in SHARDS
        ],
    output:
        out="kraken2/{sample}_kraken2_merged.report",
    conda:
        "../envs/kraken2.yaml"
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken_merge_{sample}.e",
        o="logs/kraken_merge_{sample}.o",
    shell:
        """
        combine_kreports.py \
            --only-combined \
            --no-headers \
            -r {input.reports} \
            -o {output.out} \
            > {log.o} 2> {log.e}
        """
