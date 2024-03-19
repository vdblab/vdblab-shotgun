import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)
    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


index_base = "path_to_chm13.draft_v1.0_plusY.1.bt2"
bowtie_outputs = f"bowtie/{config['sample']}.{index_base}.bam"


rule all:
    input:
        bowtie_outputs,


rule bowtie2:
    container:
        config["docker_bowtie2"]
    input:
        R1=config["R1"],
        R2=config["R2"],
        idx=multiext(
            index_base,
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bam="bowtie/{sample}.bam",
        unmapped_R1=f"bowtie/{{sample}}.without_{os.path.basename(index_base)}.R1.fastq.gz",
        unmapped_R2=f"bowtie/{{sample}}.without_{os.path.basename(index_base)}.R2.fastq.gz",
    log:
        # unfortunately we cant do lambdas on the log names
        e=f"logs/bowtie2_{{sample}}.{os.path.basename(index_base)}.e",
        o=f"logs/bowtie2_{{sample}}.{os.path.basename(index_base)}.o",
    params:
        umapped_template=lambda wildcards, output: output["unmapped_R1"].replace(
            ".R1.fastq.gz", ".R%.fastq.gz"
        ),
        db_prefix=lambda wildcards, input: os.path.splitext(
            os.path.splitext(input.idx[0])[0]
        )[0],
        extra="",  # optional parameters
    threads: 24
    shell:
        """
          (bowtie2 \
         --threads {threads} \
         -1 {input.R1} -2 {input.R2} \
         --un-conc-gz {params.umapped_template} \
         -x {params.db_prefix} \
        | samtools view -bh --threads $(({threads} - 1)) - \
        ) > {output.bam} 2> {log.e}
        """
