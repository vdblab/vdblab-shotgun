import os
import sys
from shutil import rmtree
import glob
from math import ceil


include: "common.smk"


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
    f"{config['sample']}_cazi_substrate.out",
    f"{config['sample']}_annotated_cazymes_RPM.tsv",
    abricates,
    f"{config['sample']}.cleaned_dirs",
]
if config["check_contigs"]:
    outputs.extend(
        [f"{config['sample']}_metaerg.gff", f"{config['sample']}_antismash.gbk"]
    )


nseqs = 200
nparts = 10  # will be overriden on workflow start


onstart:
    ncontigs = 0
    with open(config["assembly"], "r") as inf:
        for line in inf:
            if line.startswith(">"):
                ncontigs = ncontigs + 1
    nparts = ceil(ncontigs / nseqs)
    logger.info(f"Breaking assembly into {nparts} {nseqs}-contig chunks")


BATCHES = [f"stdin.part_{x}" for x in make_assembly_split_names(nparts)]


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
        ffn="annotation/annotation_{batch}/data/cds.ffn",
        faa="annotation/annotation_{batch}/data/cds.faa",
    resources:
        mem_mb=8 * 1024,
        runtime=45,
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
        metaerg.pl --cpus {threads} --dbdir {params.metaerg_db_dir} --outdir annotation/annotation_{wildcards.batch} --locustag {wildcards.batch} {input.assembly} --force || echo "Finished running Metaerg"
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
        runtime=6 * 60,
    threads: 16
    log:
        o="logs/antismash_{sample}.log",
    output:
        gbk="{sample}_antismash.gbk",
        outdir=directory("antismash_{sample}"),
    shell:
        """
        set +e -x
        antismash --cpus {threads} --allow-long-headers \
            --output-dir antismash_{wildcards.sample} {input.assembly} \
            --genefinding-gff {input.gff} --verbose --logfile {log.o} \
            --genefinding-tool none  # ignore contigs without genes: https://www.biostars.org/p/9539337/
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


rule split_assembly:
    """
    These dummy inputs are intended to be overwritten when importing the rule
    """
    input:
        assembly=config["assembly"],
    output:
        directory("tmp"),
        chunks=expand("tmp/{batch}.fasta", batch=BATCHES),
        assembly=temp("tmp-" + os.path.basename(config["assembly"])),
    params:
        outdir="tmp/",
        nbatches=len(BATCHES),
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
        # If you get a weird issue with missing contigs, check if its related to
        # https://github.com/shenwei356/seqkit/issues/364

        #We have to copy the input here because if we run on isabl samples,
        # the permissions don't allow the creation of the index file
        cp {input.assembly} {output.assembly}
        seqkit shuffle {output.assembly} --two-pass  |
        seqkit seq --min-len {params.minlen} --threads {threads} | \
            seqkit split2 --by-part {params.nbatches} --out-dir {params.outdir}  > {log.o} 2>> {log.e}
        """


def get_annotate_cazi_runtime(wildcards, attempt):
    return attempt * 6 * 60

def get_annotate_cazi_memory(wildcards, attempt):
    return attempt * 6 * 1024


rule annotate_CAZI_split:
    input:
        faa="annotation/annotation_{batch}/data/cds.faa",
        gff="annotation/annotation_{batch}/data/either_all_or_master.gff",
    output:
        overview="cazi_db_scan/{batch}/overview.txt",
        substrate="cazi_db_scan/{batch}/substrate.out",
        cgc="cazi_db_scan/{batch}/cgc.out",
    params:
        cazi_db=config["cazi_db"],
        contig_annotation_thresh=config["contig_annotation_thresh"],
    resources:
        mem_mb=get_annotate_cazi_memory,
        runtime=get_annotate_cazi_runtime,
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
        touch {output.substrate}
        touch {output.cgc}
        run_dbcan {input.faa} protein --out_dir cazi_db_scan/{wildcards.batch}/ -t all --db_dir /app/db -c {input.gff} --cgc_substrate --dia_cpu {threads} --hmm_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} --dbcan_thread {threads} ||  echo "WARNING: no output from this split!"
        """


rule join_metaerg_outputs:
    input:
        gff=expand(
            "annotation/annotation_{batch}/data/either_all_or_master.gff",
            batch=BATCHES,
        ),
        ffn=expand(
            "annotation/annotation_{batch}/data/cds.ffn",
            batch=BATCHES,
        ),
    output:
        gff=f"{config['sample']}_metaerg.gff",
        ffn=f"{config['sample']}_metaerg.ffn",
    container:
        config["docker_seqkit"]
    shell:
        """
        # deal with header
        head -n 1 {input.gff[0]} > {output.gff}
        for f in {input.gff}
        do
            tail -n+2 $f >> {output.gff}
        done
        for f in {input.ffn}
        do
            cat $f >> {output.ffn}
        done
        """


rule align_annotated_genes:
    input:
        ffn="annotation/annotation_{batch}/data/cds.ffn",
        r1=config["R1"],
        r2=config["R2"],
    output:
        bamfile="tmp/{batch}_aligned_reads.bam",
    container:
        config["docker_bowtie2"]
    resources:
        mem_mb=16 * 1024,
        runtime=get_annotate_cazi_runtime,
        threads=16,
        cores=16,
    params:
        bowtie_dir="./tmp_bowtie_indices_{batch}",
        bowtie_index="./tmp_bowtie_indices_{batch}/{batch}_bowtie2_index",
    shell:
        """
        mkdir -p {params.bowtie_dir}
        bowtie2-build \
            --threads {resources.threads} \
            {input.ffn} \
            {params.bowtie_index}
        bowtie2 --threads {resources.threads} -1 {input.r1} -2 {input.r2} -x {params.bowtie_index}  | samtools view -@ {resources.threads} -Sb | samtools sort -o {output.bamfile} -@ {resources.threads} 
        rm -rf {params.bowtie_dir}
        """


rule bedtools_coverage:
    input:
        gff="annotation/annotation_{batch}/data/either_all_or_master.gff",
        bamfile="tmp/{batch}_aligned_reads.bam",
    output:
        coverage="tmp/{batch}_annotated_gene_coverage.txt",
    container:
        config["docker_bedtools"]
    shell:
        """
        bedtools coverage -a {input.gff} -b {input.bamfile} > {output.coverage}
        """


rule create_RPM_counts:
    input:
        coverage="tmp/{batch}_annotated_gene_coverage.txt",
        overview="cazi_db_scan/{batch}/overview.txt",
        substrate="cazi_db_scan/{batch}/substrate.out",
        cgc="cazi_db_scan/{batch}/cgc.out",
        r1=config["R1"],
    output:
        rpm_file="tmp/{batch}_annoted_cazymes_RPM.tsv",
    conda:
        "../envs/annotate_output_parse.yaml"
    log:
        e="tmp/{batch}_create_RPM.err",
        o="tmp/{batch}_create_RPM.out",
    script:
        "../scripts/generate_RPM_annotation_files.py"


rule join_CAZI:
    input:
        overview=expand("cazi_db_scan/{batch}/overview.txt", batch=BATCHES),
        substrate=expand("cazi_db_scan/{batch}/substrate.out", batch=BATCHES),
        cgc=expand("cazi_db_scan/{batch}/cgc.out", batch=BATCHES),
        rpm=expand("tmp/{batch}_annoted_cazymes_RPM.tsv", batch=BATCHES),
    output:
        overview=f"{config['sample']}_cazi_overview.txt",
        substrate=f"{config['sample']}_cazi_substrate.out",
        cgc=f"{config['sample']}_cazi_cgc.out",
        rpm=f"{config['sample']}_annotated_cazymes_RPM.tsv",
    shell:
        """
        join_files(){{
            output_file=$1
            echo $output_file
            # deal with header
            head -n 1 $2 > $output_file
            # This weird for loop handles the fact we pass the infiles as
            # a list - (which doesn't work well in bash):
            for (( i=2; i <= "$#"; i++ )); do
                tail -n+2 ${{!i}} >> $output_file
            done
        }}
        
        join_files {output.overview} {input.overview}
        join_files {output.substrate} {input.substrate}
        join_files {output.cgc} {input.cgc}
        join_files {output.rpm} {input.rpm}
        """


rule clean_up:
    """"{sample}_metaerg.gff"  and the cazi merged output is used as an input to ensure
    this is done last.
    """
    input:
        agg_file="{sample}_metaerg.gff",
        cazi=f"{config['sample']}_cazi_overview.txt",
        rpm_file=f"{config['sample']}_annotated_cazymes_RPM.tsv",
    output:
        touch("{sample}.cleaned_dirs"),
    shell:
        """
        #find annotation/ -name "annotation_stdin.part_*" -type d | xargs --no-run-if-empty rm -r
        #rm -r tmp/
        echo "commented out stuff should add it back!!"
        """
