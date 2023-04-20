include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",


allele_mappings = f"rgi/{config['sample']}.allele_mapping_data.txt"


rule all:
    input:
        allele_mappings,
        f"rgi/{config['sample']}.allele_mapping_mqc.png",


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
        CARD_VER=os.path.basename(os.path.dirname(config["CARD_db_json"])),
        CARD_DB_DIR=lambda wildcards, input: os.path.dirname(input["db"]),
        tmpdir="rgitmp",
    conda:
        "../envs/RGI.yaml"
    container:
        config["docker_rgi"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64 * 1024,
        runtime=12 * 60,
    threads: 32
    log:
        e="rgi/{sample}.e",
        o="rgi/{sample}.o",
    shell:
        """
        # RGI uses samtools which will use TMPDIR. here we make sure that exists
        # on our cluster, certain nodes appear to have
        mkdir -p {params.tmpdir}
        # --wildcard_index {params.CARD_DB_DIR}/wildcard/index-for-model-sequences.txt
        # --wildcard_annotation {params.CARD_DB_DIR}/wildcard_database_{params.CARD_VER}.fasta

        # # --include_wildcard
        rgi load -i {input.db} \
            --card_annotation {params.CARD_DB_DIR}/card_database_{params.CARD_VER}.fasta --local
        echo "running bwt"
        rgi bwt --read_one {input.R1} \
            --read_two {input.R2} \
            --aligner bowtie2 \
            --debug \
            --output_file {params.outpre} --threads {threads} \
            --local --clean
        ls rgi
        rm -r {params.tmpdir} localDB/
        """


rule plotRGI:
    input:
        allele_mapping="rgi/{sample}.allele_mapping_data.txt",
    #        plotscript=config["rgi_plotscript"],
    output:
        heatmap="rgi/{sample}.allele_mapping_mqc.png",
    params:
        plotscript_url="https://raw.githubusercontent.com/vdblab/vdblab-pipelines/fix/annotate-binning-isabl/vdb_shotgun/scripts/plot_RGI_heatmap.R",
        plotscript_base="plot_RGI_heatmap.R",
        cov_percent=50,
    log:
        e="logs/plotRGI_{sample}.e",
        o="logs/plotRGI_{sample}.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    shell:
        """
        wget {params.plotscript_url}
        Rscript {params.plotscript_base} {input.allele_mapping} {output.heatmap} {params.cov_percent} > {log.o} 2> {log.e}
        rm {params.plotscript_base}
        """
