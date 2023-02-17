import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import validate

HTTP = HTTPRemoteProvider()
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")

validate(config, os.path.join(str(workflow.basedir), "../../config/config.schema.yaml"))

rule all:
    input:
        expand("reports/{sample}_result.tsv", sample=config['sample'])


rule OptiType_subset_fastq:
    input:
        R1=config["R1"],
        R2=config["R2"],
        alleles=HTTP.remote("https://raw.githubusercontent.com/FRED-2/OptiType/master/data/hla_reference_dna.fasta")
    output:
        R1="tmp_{sample}_R1.bam",
        R2="tmp_{sample}_R2.bam",
    container:
        config["docker_optitype"]
    threads: 16
    shell:
        """
        razers3 -i 95 -m 1 -dr 0 -o {output.R1} {input.alleles} {input.R1}
        razers3 -i 95 -m 1 -dr 0 -o {output.R2} {input.alleles} {input.R2}
        """

rule OptiType_bam2fastq:
    input:
        R1="tmp_{sample}_R1.bam",
        R2="tmp_{sample}_R2.bam",
    output:
        R1="tmp_{sample}.fastq",
    container:
        config["docker_bowtie2"]
    threads: 16
    shell:
        """
        samtools bam2fq {input.R1} > {output.R1}
        samtools bam2fq {input.R2} >> {output.R1}
        """

rule OptiType_fastq:
    # TODO move ref to resources repo
    input:
        R1="tmp_{sample}.fastq",
    output:
        results="reports/{sample}_result.tsv",
    container:
        config["docker_optitype"]
    threads: 16
    shell:
        """
        OptiTypePipeline.py -i {input.R1}  --dna --prefix {wildcards.sample} --outdir reports/ --verbose
        find .
        """
