import os

from pathlib import Path


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


if not os.path.exists("logs"):
    os.makedirs("logs")


localrules:
    all,


depths = config["depths"]
reps = config["reps"]
samples = [config["sample"]]


wildcard_constraints:
    rep="\d+",


downsampled_fastqs_r1 = expand(
    "ds{depth}_rep{rep}/ds{depth}_rep{rep}_{sample}_R1_001.fastq.gz",
    depth=depths,
    sample=samples,
    rep=reps,
)
downsampled_fastqs_r2 = expand(
    "ds{depth}_rep{rep}/ds{depth}_rep{rep}_{sample}_R2_001.fastq.gz",
    depth=depths,
    sample=samples,
    rep=reps,
)

downsampling_reports = expand("reports/downsample_stats_{sample}.tab", sample=samples)


rule all:
    input:
        downsampled_fastqs_r1,
        downsampled_fastqs_r2,
        downsampling_reports,


if len(config["R1"]) == 1:
    ds_input_R1 = config["R1"]
    ds_input_R2 = config["R2"]
else:
    ds_input_R1 = "concatenated/{sample}_R1.fastq.gz"
    ds_input_R2 = "concatenated/{sample}_R2.fastq.gz"


module utils:
    snakefile:
        "utils.smk"
    config:
        config


use rule concat_R1_R2 from utils as utils_concat_R1_R2 with:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp("concatenated/{sample}_R1.fastq.gz"),
        R2=temp("concatenated/{sample}_R2.fastq.gz"),
    log:
        e="logs/concat_r1_r2_{sample}.e",


rule downsample_fastq:
    input:
        R1=ds_input_R1,
        R2=ds_input_R2,
    output:
        R1="ds{depth}_rep{rep}/ds{depth}_rep{rep}_{sample}_R1_001.fastq.gz",
        R2="ds{depth}_rep{rep}/ds{depth}_rep{rep}_{sample}_R2_001.fastq.gz",
    params:
        tmpdir="ds{depth}_rep{rep}",
    container:
        config["docker_seqkit"]
    threads: 4  # see their docs
    resources:
        mem_mb=8000,
    shell:
        """
        set -euxo pipefail
        # assume 4 line fastqs.
        nreads=$(pigz -dc {input.R1} | awk 'NR % 4== 2' | wc -l )


        # reads_per_file is the target depth/2. eg 100 reads is 50 read pairs.
        reads_per_file=$( echo "print(round({wildcards.depth}/2))" | python3)

        if [ "$nreads" -lt "$reads_per_file" ]
        then
        echo "Insufficient depth to downsample to requested depth"
        exit 1
        fi

        # calculate the proportion needed, but avoid rounding issues by doing the *1.2 factor
        generous_proportion=$( echo "print($reads_per_file/$nreads * 1.2)" | python3)
        # get min of 1 in case the proportion *1.2 exceeds 1
        generous_proportion=$( echo "print(min($generous_proportion, 1.0))" | python3)

        # sample -p + head is seqkit's recommended approach for downsampling,
        # as it avoids having to read the whole thing into memory, which happens
        #  when you sample directly to a number of reads using -n
        # Don't ask why piping isn't working >:(
        fq_tmpfile={params.tmpdir}/tmp.fq.gz
        seqkit sample --rand-seed {wildcards.rep} --threads {threads} {input.R1} \
          -p $generous_proportion -o $fq_tmpfile
        seqkit head  -n $reads_per_file -o {output.R1} $fq_tmpfile
        rm $fq_tmpfile

        seqkit sample --rand-seed {wildcards.rep} --threads {threads} {input.R2} \
          -p $generous_proportion -o $fq_tmpfile
        seqkit head -n $reads_per_file -o {output.R2} $fq_tmpfile
        rm $fq_tmpfile

        """


rule ds_stats:
    input:
        R1=config["R1"],
        R2=config["R2"],
        R1ds=downsampled_fastqs_r1,
        R2ds=downsampled_fastqs_r2,
    output:
        stats="reports/downsample_stats_{sample}.tab",
    container:
        config["docker_seqkit"]
    shell:
        """seqkit stats --basename --tabular  --out-file {output.stats} {input}"""
