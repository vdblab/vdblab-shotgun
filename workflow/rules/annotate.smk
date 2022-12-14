import os
import sys
from shutil import rmtree
import glob


configfile: os.path.join(workflow.basedir, "../../config/config.yaml")


localrules:
    all,


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

outputs = [
    f"{config['sample']}_antismash.gbk",
    f"{config['sample']}_antismash.tab",
    f"{config['sample']}_amrfinder.tab",
    f"{config['sample']}_cazi_overview.txt",
    abricates,
    f"{config['sample']}.cleaned_dirs",
]
if config["check_contigs"]:
    outputs.extend(
        [f"{config['sample']}_metaerg.gff", f"{config['sample']}_antismash.gbk"]
    )


rule all:
    input:
        outputs,


rule annotate_orfs:
    container:
        config["docker_metaerg"]
    input:
        assembly="tmp/{batch}.fasta",
    output:
        gff="annotation/annotation_{batch}/data/either_all_or_master.gff",
    resources:
        mem_mb=8 * 1024,
        runtime=2 * 60,
    threads: 4
    params:
        metaerg_db_dir=config["metaerg_db_dir"],
    shell:
        """
        # turn off strict so we don't fail even if we have the gff file.
        # currently the output_report.pl script is failing
        # see issues https://github.com/xiaoli-dong/metaerg/pull/38 and
        # https://github.com/xiaoli-dong/metaerg/issues/12
        set -e
        metaerg.pl --cpus {threads} --dbdir {params.metaerg_db_dir} --outdir annotation/annotation_{wildcards.batch} {input.assembly} --force || echo "Finished running Metaerg"
        # if metaerg successfully packaged everything up
        if [ -f "annotation/annotation_{wildcards.batch}/data/master.gff.txt" ]
        then
            mv annotation/annotation_{wildcards.batch}/data/master.gff.txt {output.gff}
        else
            # if it successed but failed at output_report.pl, no need to do anything
            echo "sample likely failed at output_report.pl but gff should be present"
            mv annotation/annotation_{wildcards.batch}/data/all.gff {output.gff}
        fi
        """


rule antismash:
    # note that we don't require the web index as antismash can fail on small samples.
    container:
        config["docker_antismash"]
    input:
        assembly=config["assembly"],
        gff="{sample}_metaerg.gff",
    resources:
        mem_mb=16 * 1024,
        runtime=3 * 60,
    threads: 16
    output:
        gbk="{sample}_antismash.gbk",
    shell:
        """
        set +e -x
        antismash --cpus {threads} --allow-long-headers \
            --output-dir antismash_{wildcards.sample} {input.assembly} \
            --genefinding-gff {input.gff}
        exitcode=$?
        if [ ! $exitcode -eq 0 ]
        then
            echo "Error running antismash; this can occur if the assembly is very fragmented/poor"
            touch {output.gbk}

        else
            cat antismash_{wildcards.sample}/*.gbk > {output.gbk}
        fi
        """


rule tabulate_antismash:
    input:
        gbk="{sample}_antismash.gbk",
    output:
        tab="{sample}_antismash.tab",
    container:
        config["docker_biopython"]
    script:
        "../scripts/parse_antismash_gbk.py"


rule annotate_abricate:
    input:
        assembly=config["assembly"],
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
        assembly=config["assembly"],
    output:
        amr="{sample}_amrfinder.tab",
    resources:
        mem_mb=4000,
        runtime=3 * 60,
    container:
        config["docker_amrfinder"]
    shell:
        """
        amrfinder -n {input.assembly}  --plus  > {output.amr}
        """


checkpoint split_assembly:
    """
    These dummy inputs are intended to be overwritten when importing the rule
    """
    input:
        R1=config["assembly"],
    output:
        directory("tmp"),
    params:
        outdir="tmp/",
        nseqs=200,
        minlen=config["contig_annotation_thresh"],
    container:
        "docker://pegi3s/seqkit:2.3.0"
    threads: 4  # see their docs
    resources:
        mem_mb=8000,
    log:
        e="logs/split_assembly.e",
        o="logs/split_assembly.o",
    shell:
        """
        seqkit shuffle {input.R1} --two-pass  |
        seqkit seq --min-len {params.minlen} --threads {threads} | \
            seqkit split2 --by-size {params.nseqs} --out-dir {params.outdir}  > {log.o} 2>> {log.e}
        """


rule annotate_CAZI_split:
    input:
        assembly="tmp/{batch}.fasta",
    output:
        overview="cazi_db_scan/{batch}/overview.txt",
    params:
        cazi_db=config["cazi_db"],
        contig_annotation_thresh=config["contig_annotation_thresh"],
    resources:
        mem_mb=4 * 1024,
        runtime=90,
    container:
        config["docker_dbcan"]
    threads: 2
    shell:
        """
        # turn off strict so we don't fail on empty hmmer outputs within run_dbcan
        set -e
        # touch this file so it exists even if this fails
        # if you can figure out checkpoints that can deal with missing files,
        # you win!
        touch {output.overview}
        run_dbcan {input.assembly} meta --out_dir cazi_db_scan/{wildcards.batch}/ -t all --db_dir /app/db --dia_cpu {threads} --hmm_cpu {threads} --eCAMI_jobs {threads} --tf_cpu {threads} --stp_cpu {threads}   ||  echo "WARNING: no output from this split!"
        """


def aggregate_cazi_results(wildcards):
    """
    aggregate the file names of the random number of files
    generated at the scatter step
    """
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output[0]
    return expand(
        "cazi_db_scan/{batch}/overview.txt",
        batch=glob_wildcards(os.path.join(checkpoint_output, "{batch}.fasta")).batch,
    )


def aggregate_metaerg_results(wildcards):
    """
    aggregate the file names of the random number of files
    generated at the scatter step
    """
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output[0]
    return expand(
        "annotation/annotation_{batch}/data/either_all_or_master.gff",
        batch=glob_wildcards(os.path.join(checkpoint_output, "{batch}.fasta")).batch,
    )


rule join_CAZI:
    input:
        aggregate_cazi_results,
    output:
        f"{config['sample']}_cazi_overview.txt",
    shell:
        """
        # deal with header
        head -n 1 {input[0]} > {output}
        for f in {input}
        do
            tail -n+2 $f >> {output}
        done
        """


rule join_gffs:
    input:
        aggregate_metaerg_results,
    output:
        "{sample}_metaerg.gff",
    shell:
        """
        # deal with header
        head -n 1 {input[0]} > {output}
        for f in {input}
        do
            tail -n+2 $f >> {output}
        done
        """


rule clean_up:
    """"{sample}_metaerg.gff" is used as an input to ensure
    that step is done before we clean.
    """
    input:
        agg_file="{sample}_metaerg.gff",
    output:
        touch("{sample}.cleaned_dirs"),
    shell:
        """
        find annotation/ -name "annotation_stdin.part_*" -type d | xargs --no-run-if-empty rm -r
        rm -r tmp/
        """
