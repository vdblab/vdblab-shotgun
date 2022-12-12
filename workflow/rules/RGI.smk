include: "common.smk"


configfile: os.path.join(workflow.basedir, "../../config/config.yaml")


allele_mappings = f"rgi/{config['sample']}.allele_mapping_data.txt"


rule all:
    input:
        allele_mappings,


rule RGI:
    input:
        R1=config["R1"],
        R2=config["R2"],
        db=config["CARD_db_json"],
    output:
        allele_mapping="rgi/{sample}.allele_mapping_data.txt",
        artifacts="rgi/{sample}.artifacts_mapping_stats.txt",
        mapstats="rgi/{sample}.overall_mapping_stats.txt",
        refmapstats="rgi/{sample}.reference_mapping_stats.txt",
    params:
        outpre=lambda wildcards, output: output.allele_mapping.split(
            ".allele_mapping"
        )[0],
        CARD_VER="3.2.4",
        CARD_DB_DIR=lambda wildcards, input: os.path.dirname(input["db"]),
    conda:
        "../envs/RGI.yaml"
    container:
        config["docker_rgi"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64 * 1024,
        runtime="6:00",
    threads: 24
    log:
        e="rgi/{sample}.e",
        o="rgi/{sample}.o",
    shell:
        """
        #rgi load -i {input.db} --card_annotation {params.CARD_DB_DIR}/card_database_v{params.CARD_VER}.fasta --local

        rgi load -i {input.db} --wildcard_annotation {params.CARD_DB_DIR}/wildcard_database_v{params.CARD_VER}.fasta  \
            --wildcard_index {params.CARD_DB_DIR}/wildcard/index-for-model-sequences.txt \
            --card_annotation {params.CARD_DB_DIR}/card_database_v{params.CARD_VER}.fasta --local
        echo "running bwt"
        rgi bwt --read_one {input.R1} \
            --read_two {input.R2} \
            --aligner bowtie2 \
            --debug \
            --output_file {params.outpre} --threads {threads} --include_wildcard \
            --local --clean
        ls rgi
        """
