import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")


validate(
    config,
    os.path.join(str(workflow.current_basedir), "../../config/config.schema.yaml"),
)


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)
    if not os.path.exists("logs"):
        os.makedirs("logs")


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


rule s01_bowtie2:
    container:
        config["docker_bowtie2"]
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
        bam=temp(f"01-bowtie/{config['sample']}.{bowtie2_human_db_name}.bam"),
        unmapped_R1=temp(
            f"01-bowtie/{config['sample']}.without_{bowtie2_human_db_name}.R1.fastq.gz"
        ),
        unmapped_R2=temp(
            f"01-bowtie/{config['sample']}.without_{bowtie2_human_db_name}.R2.fastq.gz"
        ),
    log:
        e=f"logs/bowtie2_{config['sample']}.{bowtie2_human_db_name}.e",
        o=f"logs/bowtie2_{config['sample']}.{bowtie2_human_db_name}.o",
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

rule s02_snapalign:
    container:
        config["docker_snap"]
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
        bam=temp(f"02-snap/{{sample}}.{snap_human_db_name}.bam"),
    log:
        e=f"logs/snap_{{sample}}.{snap_human_db_name}.e",
        o=f"logs/snap_{{sample}}.{snap_human_db_name}.o",
    params:
        db_prefix=lambda wildcards, input: os.path.dirname(input.idx_genome),
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 10),
        runtime=48 * 60,
    threads: 24  # Use at least two threads
    shell:
        """
        snap-aligner paired {params.db_prefix} \
        {input.R1} {input.R2}  -o {output.bam} -t {threads} -xf 2.0
        """


rule s03_get_unmapped:
    """ see mgen/10.1099/mgen.0.000393
       > This two-stage approach first classified, and then discarded,
         ‘human’ reads using one method, and then performed a second round of
          classification using a second method. In this way, a
         method with high precision could be supplemented by a
        method with high sensitivity, maximizing the utility of both.

    # uses logic from https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd

    """
    container:
        config["docker_bowtie2"]
    input:
        bam=f"02-snap/{{sample}}.{snap_human_db_name}.bam",
    output:
        unmapped_R1=temp(f"03-nohuman/{{sample}}.without_{snap_human_db_name}.R1.fastq"),
        unmapped_R2=temp(f"03-nohuman/{{sample}}.without_{snap_human_db_name}.R2.fastq"),
        flagstat=temp(f"02-snap/{{sample}}.{snap_human_db_name}.bam.flagstat"),
    log:
        e=f"logs/get_unmapped_human_{{sample}}.e",
        o=f"logs/get_unmapped_human_{{sample}}.o",
    threads: 8
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        # R1 unmapped, R2 mapped
        samtools view -u -f 4 -F 264 {input.bam} > tmp_unmap_map_{wildcards.sample}.bam
        # R1 mapped, R2 unmapped
        samtools view -u -f 8 -F 260 {input.bam} > tmp_map_unmap_{wildcards.sample}.bam
        # R1 & R2 unmapped
        samtools view -u -f 12 -F 256 {input.bam} > tmp_unmap_unmap_{wildcards.sample}.bam

        samtools merge -u tmp_unmapped_{wildcards.sample}.bam tmp_unmap_map_{wildcards.sample}.bam tmp_map_unmap_{wildcards.sample}.bam tmp_unmap_unmap_{wildcards.sample}.bam
        samtools flagstat tmp_unmapped_{wildcards.sample}.bam

        # note this outputs uncompressed only
        bamToFastq -i tmp_unmapped_{wildcards.sample}.bam -fq {output.unmapped_R1} -fq2 {output.unmapped_R2}
        rm tmp*_{wildcards.sample}.bam
        """






use rule s01_bowtie2  as s04_bowtie2_mouse with:
    input:
        R1=rules.s03_get_unmapped.output.unmapped_R1,
        R2=rules.s03_get_unmapped.output.unmapped_R2,
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
        bam=temp(f"04-bowtie/{{sample}}.{bowtie2_mouse_db_name}.bam"),
        unmapped_R1=temp(
            f"04-bowtie/{{sample}}.without_{bowtie2_mouse_db_name}.R1.fastq.gz"
        ),
        unmapped_R2=temp(
            f"04-bowtie/{{sample}}.without_{bowtie2_mouse_db_name}.R2.fastq.gz"
        ),
    log:
        e=f"logs/bowtie2_{{sample}}.{bowtie2_mouse_db_name}.e",
        o=f"logs/bowtie2_{{sample}}.{bowtie2_mouse_db_name}.o",


use rule s02_snapalign as s04_snapalign_mouse with:
    input:
        R1=rules.s04_bowtie2_mouse.output.unmapped_R1,
        R2=rules.s04_bowtie2_mouse.output.unmapped_R2,
        idx_genome=f"{config['snap_mouse_index_dir']}Genome",
    output:
        bam=temp(f"05-snap/{{sample}}.{snap_mouse_db_name}.bam"),
    log:
        e=f"logs/snap_{{sample}}.{snap_mouse_db_name}.e",
        o=f"logs/snap_{{sample}}.{snap_mouse_db_name}.o",


use rule s03_get_unmapped  as s06_get_unmapped_human_mouse with:
    input:
        bam=f"05-snap/{{sample}}.{snap_mouse_db_name}.bam",
    output:
        unmapped_R1=temp(f"06-nohuman-nomouse/{{sample}}.R1.fastq"),
        unmapped_R2=temp(f"06-nohuman-nomouse/{{sample}}.R2.fastq"),
        flagstat=f"05-snap/{{sample}}.{snap_mouse_db_name}.bam.flagstat",
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
