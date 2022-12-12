index_base = "path/to/referece"
# rule all:
#     input:
#         unmapped_R1=f"snap/{config['sample']}.without_{index_base}.R1.fastq.gz",
#         unmapped_R2=f"snap/{config['sample']}.without_{index_base}.R2.fastq.gz",


rule snapalign:
    container:
        config["docker_snap"]
    input:
        R1=config["R1"],
        R2=config["R2"],
        idx_genome=f"path/to/reference/Genome",
    output:
        bam="snap/{sample}.bam",
    log:
        e=f"logs/snap_{{sample}}.{index_base}.e",
        o=f"logs/snap_{{sample}}.{index_base}.o",
    params:
        db_prefix=lambda wildcards, input: os.path.dirname(input.idx_genome),
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * (max(input.size // 1000000, 1024) * 10),
        runtime="48:00",
    threads: 24  # Use at least two threads
    shell:
        """
        snap-aligner paired {params.db_prefix} \
        {input.R1} {input.R2}  -o {output.bam} -t {threads} -xf 2.0
        """


rule get_unmapped:
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
        bam="snap/{sample}.bam",
    output:
        unmapped_R1=f"snap/{{sample}}.without_{index_base}.R1.fastq",
        unmapped_R2=f"snap/{{sample}}.without_{index_base}.R2.fastq",
        flagstat="snap/{sample}.bam.flagstat",
    log:
        e=f"logs/snap_{{sample}}.{os.path.basename(index_base)}.e",
        o=f"logs/snap_{{sample}}.{os.path.basename(index_base)}.o",
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
