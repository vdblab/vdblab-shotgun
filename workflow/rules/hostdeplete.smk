import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


configfile: os.path.join(workflow.basedir, "../runconfig.yaml")


config["pipeline_version"] = get_pipeline_version()

# validate(config, os.path.join(workflow.current_basedir, "../runconfig.schema.yaml"))


onstart:
    with open("config_used.yml", "w") as outfile:
        yaml.dump(config, outfile)
    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    hostdeplete_all,


bowtie2_human_db_name = os.path.basename(config["bowtie2_human_index_base"])
bowtie2_mouse_db_name = os.path.basename(config["bowtie2_mouse_index_base"])
snap_human_db_name = os.path.basename(os.path.dirname(config["snap_human_index_dir"]))
snap_mouse_db_name = os.path.basename(os.path.dirname(config["snap_mouse_index_dir"]))


human_bowtie_outputs = f"01-bowtie/{config['sample']}.{bowtie2_human_db_name}.bam"
human_snap_outputs = f"02-snap/{config['sample']}.{snap_human_db_name}.bam"
mouse_bowtie_outputs = f"04-bowtie/{config['sample']}.{bowtie2_mouse_db_name}.bam"
mouse_snap_outputs = f"05-snap/{config['sample']}.{snap_mouse_db_name}.bam"
cleaned_reads_R1 = f"06-nohuman-nomouse/{config['sample']}.R1.fastq.gz"
cleaned_reads_R2 = f"06-nohuman-nomouse/{config['sample']}.R2.fastq.gz"
table = f"{config['sample']}_depletion.stats"


rule hostdeplete_all:
    input:
        human_bowtie_outputs,
        mouse_bowtie_outputs,
        human_snap_outputs,
        mouse_snap_outputs,
        cleaned_reads_R2,
        table,


module bowtie2:
    snakefile:
        "bowtie2.smk"
    config:
        config


module snap:
    snakefile:
        "snap.smk"
    config:
        config


use rule bowtie2 from bowtie2 as bowtie2_human with:
    input:
        R1=config["R1"],
        R2=config["R2"],
        idx=multiext(
            config["bowtie2_human_index_base"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bam=f"01-bowtie/{config['sample']}.{bowtie2_human_db_name}.bam",
        unmapped_R1=f"01-bowtie/{config['sample']}.without_{bowtie2_human_db_name}.R1.fastq.gz",
        unmapped_R2=f"01-bowtie/{config['sample']}.without_{bowtie2_human_db_name}.R2.fastq.gz",
    log:
        e=f"logs/bowtie2_{config['sample']}.{bowtie2_human_db_name}.e",
        o=f"logs/bowtie2_{config['sample']}.{bowtie2_human_db_name}.o",


use rule snapalign from snap as snapalign_human with:
    input:
        # using the direct syntax breaks when importing module
        # presumably some namespace clashing with bowtie2 when we redefine it
        # the others below seem to be fine
        #        R1=rules.bowtie2_human.output.unmapped_R1,
        #        R2=rules.bowtie2_human.output.unmapped_R2,
        R1=f"01-bowtie/{{sample}}.without_{bowtie2_human_db_name}.R1.fastq.gz",
        R2=f"01-bowtie/{{sample}}.without_{bowtie2_human_db_name}.R2.fastq.gz",
        idx_genome=f"{config['snap_human_index_dir']}Genome",
    output:
        bam=f"02-snap/{{sample}}.{snap_human_db_name}.bam",
    log:
        e=f"logs/snap_{{sample}}.{snap_human_db_name}.e",
        o=f"logs/snap_{{sample}}.{snap_human_db_name}.o",


use rule get_unmapped from snap as get_unmapped_human with:
    input:
        bam=rules.snapalign_human.output.bam,
    output:
        unmapped_R1=f"03-nohuman/{{sample}}.without_{snap_human_db_name}.R1.fastq",
        unmapped_R2=f"03-nohuman/{{sample}}.without_{snap_human_db_name}.R2.fastq",
        flagstat=f"{rules.snapalign_human.output.bam}.flagstat",
    log:
        e=f"logs/get_unmapped_human_{{sample}}.e",
        o=f"logs/get_unmapped_human_{{sample}}.o",


use rule bowtie2 from bowtie2 as bowtie2_mouse with:
    input:
        R1=rules.get_unmapped_human.output.unmapped_R1,
        R2=rules.get_unmapped_human.output.unmapped_R2,
        idx=multiext(
            config["bowtie2_mouse_index_base"],
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bam=f"04-bowtie/{{sample}}.{bowtie2_mouse_db_name}.bam",
        unmapped_R1=f"04-bowtie/{{sample}}.without_{bowtie2_mouse_db_name}.R1.fastq.gz",
        unmapped_R2=f"04-bowtie/{{sample}}.without_{bowtie2_mouse_db_name}.R2.fastq.gz",
    log:
        e=f"logs/bowtie2_{{sample}}.{bowtie2_mouse_db_name}.e",
        o=f"logs/bowtie2_{{sample}}.{bowtie2_mouse_db_name}.o",


use rule snapalign from snap as snapalign_mouse with:
    input:
        R1=rules.bowtie2_mouse.output.unmapped_R1,
        R2=rules.bowtie2_mouse.output.unmapped_R2,
        idx_genome=f"{config['snap_mouse_index_dir']}Genome",
    output:
        bam=f"05-snap/{{sample}}.{snap_mouse_db_name}.bam",
    log:
        e=f"logs/snap_{{sample}}.{snap_mouse_db_name}.e",
        o=f"logs/snap_{{sample}}.{snap_mouse_db_name}.o",


use rule get_unmapped from snap as get_unmapped_human_mouse with:
    input:
        bam=rules.snapalign_mouse.output.bam,
    output:
        unmapped_R1=f"06-nohuman-nomouse/{{sample}}.R1.fastq",
        unmapped_R2=f"06-nohuman-nomouse/{{sample}}.R2.fastq",
        flagstat=f"{rules.snapalign_mouse.output.bam}.flagstat",
    log:
        e=f"logs/get_unmapped_human_mouse_{{sample}}.e",
        o=f"logs/get_unmapped_human_mouse_{{sample}}.o",


rule tally_depletion:
    input:
        bam01=f"01-bowtie/{{sample}}.{bowtie2_human_db_name}.bam",
        bam02=f"02-snap/{{sample}}.{snap_human_db_name}.bam",
        bam04=f"04-bowtie/{{sample}}.{bowtie2_mouse_db_name}.bam",
        bam05=f"05-snap/{{sample}}.{snap_mouse_db_name}.bam",
    output:
        table="hostdepleted/{sample}_hostdepletion.stats",
    container:
        config["docker_bowtie2"]
    threads: 1
    shell:
        """
        human_bowtie=$(samtools view  -c  {input.bam01})
        human_snap=$(samtools view  -c  {input.bam02})
        mouse_bowtie=$(samtools view  -c  {input.bam04})
        mouse_snap=$(samtools view  -c  {input.bam05})
        human_aligned_bowtie=$(samtools view -F 0x04 -c  {input.bam01})
        human_aligned_snap=$(samtools view -F 0x04 -c  {input.bam02})
        mouse_aligned_bowtie=$(samtools view -F 0x04 -c  {input.bam04})
        mouse_aligned_snap=$(samtools view -F 0x04 -c  {input.bam05})
        echo -e "sample\tbowtie2_human\tbowtie2_human_aligned\tsnap_human\tsnap_human_aligned\tbowtie2_mouse\tbowtie2_mouse_aligned\tsnap_mouse\tsnap_mouse_aligned" > {output.table}
        echo -e "{wildcards.sample}\t$human_bowtie\t$human_aligned_bowtie\t$human_snap\t$human_aligned_snap\t$mouse_bowtie\t$mouse_aligned_bowtie\t$mouse_snap\t$mouse_aligned_snap" >> {output.table}
        """


rule xHLA:
    input:
        bam=f"01-bowtie/{{sample}}.{bowtie2_human_db_name}.bam",
    output:
        results="xHLA/{sample}.json",
    container:
        config["docker_hla"]
    threads: 1
    shell:
        """ run.py \
        --sample_id {wildcards.sample} --input_bam_path {input.bam} \
        --output_path xHLA
        find .
        """
