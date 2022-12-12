import os
import sys


configfile: os.path.join(workflow.basedir, "runconfig.yaml")


assembly_base = os.path.basename(config["assembly"])


localrules:
    localize,


if not os.path.exists("logs"):
    os.makedirs("logs")

abricates = expand(
    "{sample}_abricate_{tool}.tab",
    sample=config["sample"],
    tool=[
        "argannot",
        "card",
        "ecoh",
        "ecoli_vf",
        "megares",
        "ncbi",
        "plasmidfinder",
        "resfinder",
        "vfdb",
    ],
)


rule all:
    input:
        f"{config['sample']}_metaerg.gff",
        f"{config['sample']}_antismash.gbk",
        f"{config['sample']}_amrfinder.tab",
        "cazi_db_scan/done",
        abricates,


rule localize:
    input:
        config["assembly"],
    output:
        assembly_base,
    shell:
        """
        cp {input} {output}
        """


rule annotate_orfs:
    container:
        config["docker_metaerg"]
    input:
        assembly=assembly_base,
    output:
        gff="{sample}_metaerg.gff",
    resources:
        mem_mb=16 * 1024,
    threads: 32
    params:
        metaerg_db_dir=config["metaerg_db_dir"],
        check_contigs=1 if config["check_contigs"] else 0,
        contig_annotation_thresh=config["contig_annotation_thresh"],
    shell:
        """
        set -x
        # get the first (longest) contig name that contains the length
        maxlength=$(cat {input.assembly} |head -n1 | cut -f 4 -d _)
        if [ "$maxlength" -lt "{params.contig_annotation_thresh}" ] && [ "{params.check_contigs}" -eq "1" ]
        then
            echo "longest contig less than 1000bp; skipping annotation"
            touch {output.gff}
        else
            rm -r  annotation
            metaerg.pl --cpus {threads} --dbdir {params.metaerg_db_dir} --outdir annotation {input.assembly}
            mv annotation/data/all.gff {output.gff}
        fi
        """


rule antismash:
    container:
        config["docker_antismash"]
    input:
        assembly=assembly_base,
        gff=rules.annotate_orfs.output.gff,
    resources:
        mem_mb=16000,
    threads: 16
    output:
        gbk="{sample}_antismash.gbk",
    shell:
        """
        set +e -x
        antismash --cpus {threads}  --output-dir antismash {input.assembly} --genefinding-gff {input.gff}
        exitcode=$?
        if [ ! $exitcode -eq 0 ]
        then
            echo "Error running antismash; this can occur if the assembly is very fragmented/poor"
            touch {output.gbk}

        else
            cat antismash/*.gbk > {output.gbk}
        fi
        """


rule tabulate_antismash:
    input:
        gbk="{sample}_antismash.gbk",
    output:
        tab="{sample}_antismash.tab",
    script:
        "scripts/parse_antismash_gbk.py"


rule annotate_abricate:
    input:
        assembly=assembly_base,
    output:
        out="{sample}_abricate_{tool}.tab",
    container:
        config["docker_abricate"]
    resources:
        mem_mb=4000,
    shell:
        """
        abricate --db {wildcards.tool} {input.assembly} > {output.out}
        """


rule annotate_AMR:
    input:
        assembly=assembly_base,
    output:
        amr="{sample}_amrfinder.tab",
    resources:
        mem_mb=4000,
    container:
        config["docker_amrfinder"]
    shell:
        """
        amrfinder -n {input.assembly}  --plus  > {output.amr}
        # /tmp/amrfinder/data/2022-04-04.1/
        """


rule annotate_CAZI:
    input:
        assembly=assembly_base,
    output:
        done="cazi_db_scan/done",
    params:
        cazi_db=config["cazi_db"],
        check_contigs=1 if config["check_contigs"] else 0,
        contig_annotation_thresh=config["contig_annotation_thresh"],
    resources:
        mem_mb=8000,
    container:
        config["docker_dbcan"]
    threads: 16
    shell:
        """
        set -x
        maxlength=$(cat {input.assembly} |head -n1 | cut -f 4 -d _)
        if [ "$maxlength" -lt "{params.contig_annotation_thresh}" ] && [ "{params.check_contigs}" -eq "1" ]
        then
            echo "longest contig less than 1000bp; skipping CAZI annotation"
        else
            run_dbcan {input.assembly} meta --out_dir cazi_db_scan/ -t all --db_dir /app/db --dia_cpu {threads} --hmm_cpu {threads} --eCAMI_jobs {threads} --tf_cpu {threads} --stp_cpu {threads}
        fi
        touch cazi_db_scan/done
        """
