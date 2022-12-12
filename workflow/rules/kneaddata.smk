
include: "common.smk"


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])
TMPDIR = Path(os.environ["TMPDIR"])


# Remove human reads contamination using Kneaddata
rule kneaddata_run:
    input:
        R1=config["R1"],  #"trimmed/{sample}_shard{shard}_trim_R1.fastq.gz",
        R2=config["R2"],  #"trimmed/{sample}_shard{shard}_trim_R2.fastq.gz",
        bmtagger_human=config["bmtagger_human"],
        bmtagger_mouse=config["bmtagger_mouse"],
    output:
        zip1=temp("kneaddata/{sample}_shard{shard}_knead_paired_1.fastq"),
        zip2=temp("kneaddata/{sample}_shard{shard}_knead_paired_2.fastq"),
    params:
        out_prefix="{sample}_shard{shard}_knead",
        out_dir=os.path.join(TMPDIR, "{sample}_shard{shard}"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000,
        runtime="12:00",
    container:
        config["docker_kneaddata"]
    conda:
        "../envs/kneaddata.yaml"
    log:
        e="logs/{sample}_shard{shard}_knead.e",
        o="kneaddata/{sample}_shard{shard}_knead.log",
    threads: 24
    shell:
        """
        echo START > {log.o}
        kneaddata \
            -i {input.R1} -i {input.R2} \
            -o {params.out_dir} \
            -db {input.bmtagger_human} \
            -db {input.bmtagger_mouse} \
            --output-prefix {params.out_prefix} \
            --bypass-trim \
            --run-bmtagger \
            -t {threads} \
            --log {log.o} 2> {log.e}
        echo END >> {log.o}
        mv {params.out_dir}/{params.out_prefix}_paired_1.fastq kneaddata
        mv {params.out_dir}/{params.out_prefix}_paired_2.fastq kneaddata
        """


rule kneaddata_stats:
    # input fastqs are just used as a trigger, cause the count tool
    # takes an input directory not the actual log files
    input:
        trigger=expand(
            "kneaddata/{{sample}}_shard{shard}_knead_paired_2.fastq", shard=SHARDS
        ),
    output:
        table="hostdepleted/{sample}.hostdepletion.stats",
    resources:
        mem_mb=2000,
    container:
        config["docker_kneaddata"]
    conda:
        "../envs/kneaddata.yaml"
    threads: 1
    log:
        e="logs/kneaddata_{sample}.e",
    shell:
        """
        kneaddata_read_count_table \
            --input kneaddata \
            --output {output.table} \
            2> {log.e}
        """
