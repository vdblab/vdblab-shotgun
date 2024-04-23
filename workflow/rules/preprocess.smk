import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])
TMPDIR = Path(os.environ["TMPDIR"])


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


sortmerna_outputs = f"reports/{config['sample']}_sortmerna.merged.log"
readdirs = [1,2] if config["lib_layout"] == "paired" else [1]
cleaned_fastqs = expand(
    "hostdepleted/{sample}_R{read_dir}.fastq.gz",
    sample=config["sample"],
    read_dir=readdirs,
)


wildcard_constraints:
    sample="[^/]+",


rule all:
    input:
        clean_fastqs=cleaned_fastqs,
        hostdeplete_stats_mqc=f"reports/{config['sample']}_hostdeplete.stats.summary_mqc.tsv",
        fastqcs_mqc=expand(f"reports/{config['sample']}_R{{rd}}_fastqc.html", rd=readdirs),
        sortmerna_blast=f"sortmerna/{config['sample']}_sortmerna.blast.gz",
        host_reads=expand(f"host/{config['sample']}_all_host_reads_R{{rd}}.fastq.gz",  rd=readdirs),


# note:  we could make the concatenation conditional on how many libraries we
# have, but a nice side-effect of creating a temporary file is that the
# resulting filename is uniform.  For tools like fastqc, this means we can
# predict the outputs much easier
module utils:
    snakefile:
        "utils.smk"
    config:
        config


use rule concat_lanes_fix_names from utils as utils_concat_lanes_fix_names with:
    input:
        unpack(get_input_fastqs),
    output:
        R1=temp("concatenated/{sample}_R1.fastq.gz"),
        R2=temp("concatenated/{sample}_R2.fastq.gz"),
    log:
        e="logs/concat_r1_r2_{sample}.e",


# Initial FastQC to see the quality
rule initial_fastqc_run:
    input:
        rules.utils_concat_lanes_fix_names.output
        # R1="concatenated/{sample}_R1.fastq.gz",
        # R2="concatenated/{sample}_R2.fastq.gz",
    output:
        reports=expand("reports/{{sample}}_R{rd}_fastqc.html", rd=readdirs)
    threads: 4
    container:
        config["docker_fastqc"]
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=4000,
    log:
        e="logs/fastqc_{sample}.e",
        o="logs/fastqc_{sample}.o",
    shell:
        """
        fastqc \
            --outdir reports/ \
            --threads {threads} \
            --noextract {input} \
            > {log.o} 2> {log.e}
        """


rule bbmap_dedup:
    input:
        rules.utils_concat_lanes_fix_names.output,
#        R1="concatenated/{sample}_R1.fastq.gz",
#        R2="concatenated/{sample}_R2.fastq.gz",
    output:
        reads=temp(expand("dedup/{{sample}}_R{rd}.fastq.gz", rd=readdirs)),
        dedup_stats="reports/{sample}_dedup.stats",
    threads: 8
    params:
        allowed_subs=3,
        flags=bbmap_dedup_params_flags,
        inputstring = "in1={input[0]} out2={input[1]}" if is_paired() else "in1={input[0]}",
        outputstring = "out1={output[0]} out2={output[1]}" if is_paired() else "out1={output[0]}",

    resources:
        mem_mb=lambda wildcards, input, attempt: attempt
        * (max(input.size // 1000000, 1024) * 16),
        runtime=24 * 60,
    log:
        # this is annoying but we want to be able to extract the stats from
        # the logs, which we can't do without the logs as a file. Perhaps
        # tee-ing would work, if you can do it with stderr
        "logs/bbmap_dedup_{sample}.log",
    container:
        config["docker_bbtools"]
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        clumpify.sh \
            {parms.inputstring} \
            {parms.outputstring} \
            {params.flags} \
            subs={params.allowed_subs} \
            t={threads} \
            -Xmx{resources.mem_mb}M \
            -eoom \
            &> {log}

        echo -e "BBmap clumpify metric\t{wildcards.sample}" > {output.dedup_stats}
        cat {log} | \
            grep "Reads In\:" -A4 | \
            sed "s|\: *|\t|g" \
            >> {output.dedup_stats}

        # clumpify tends to fail silently. We both catch general cases resulting in an empty file
        # and cases of misspecified platform (eg trying to do optical deduplication on SRA samples without coordinates)
        if [ ! -s "{output.R1}" ]
        then
          echo "output file after dedup is empty"
          exit 1
        else
          nlines=$(zcat {output.R1} | wc -l)
          if [ "$nlines" -lt 4 ]
          then
            echo "output file after dedup is empty"
            exit 1
          fi
        fi
        """


use rule split_fastq from utils as utils_split_fastq with:
    input:
        R1=lambda wildcards: files_to_split(
            wildcards, dedup=config["dedup_reads"], read_dir=1
        ),
        R2=lambda wildcards: files_to_split(
            wildcards, dedup=config["dedup_reads"], read_dir=2
        ),
    output:
        R1=temp(expand("split_fastq/{{sample}}_R1.part_{shard}.fastq.gz", shard=SHARDS)),
        R2=temp(expand("split_fastq/{{sample}}_R2.part_{shard}.fastq.gz", shard=SHARDS)),
    log:
        e="logs/split_fastq_{sample}.e",
        o="logs/split_fastq_{sample}.o",
    params:
        outdir="split_fastq",
        nshards=config["nshards"],


# Trim adapters with BBMap
rule bbmap_run:
    input:
        R1=lambda wildcards: files_to_trim(
            wildcards,
            nshards=config["nshards"],
            dedup=config["dedup_reads"],
            read_dir=1,
        ),
        R2=lambda wildcards: files_to_trim(
            wildcards,
            nshards=config["nshards"],
            dedup=config["dedup_reads"],
            read_dir=2,
        ),
        adapter=config["adapters_fasta"],
    output:
        out_R1=temp("trimmed/{sample}_shard{shard}_trim_R1.fastq.gz"),
        out_R2=temp("trimmed/{sample}_shard{shard}_trim_R2.fastq.gz"),
        rm_R1=temp("trimmed/{sample}_shard{shard}_discard_R1.fastq.gz"),
        rm_R2=temp("trimmed/{sample}_shard{shard}_discard_R2.fastq.gz"),
        stats=temp("trimmed/{sample}_shard{shard}_trimmingAQ.txt"),
    threads: 8
    params:
        filter_params=config["filter_params"],
    resources:
        mem_mb=4000,
    container:
        config["docker_bbtools"]
    conda:
        "../envs/bbmap.yaml"
    log:
        "logs/bbmap_{sample}_shard{shard}.e",
    shell:
        """
        bbduk.sh -Xmx{resources.mem_mb}m \
            ordered \
            in={input.R1} in2={input.R2} \
            out={output.out_R1} out2={output.out_R2} \
            outm={output.rm_R1} outm2={output.rm_R2} \
            ref={input.adapter} \
        {params.filter_params} \
            stats={output.stats} \
            threads={threads} \
            2> {log}
        """


bowtie2_human_db_name = os.path.basename(config["bowtie2_human_index_base"])
bowtie2_mouse_db_name = os.path.basename(config["bowtie2_mouse_index_base"])
snap_human_db_name = os.path.basename(os.path.dirname(config["snap_human_index_dir"]))
snap_mouse_db_name = os.path.basename(os.path.dirname(config["snap_mouse_index_dir"]))

if is_paired():
    module hostdeplete:
        snakefile:
            "hostdeplete.smk"
        config:
             config
else:
    module hostdeplete:
        snakefile:
            "hostdeplete-se.smk"
        config:
            config


module bowtie2:
    snakefile:
        "bowtie2.smk"
    config:
        config


use rule * from hostdeplete as bt_*


use rule bowtie2 from bowtie2 as bowtie_human with:
    input:
        R1="trimmed/{sample}_shard{shard}_trim_R1.fastq.gz",
        R2="trimmed/{sample}_shard{shard}_trim_R2.fastq.gz",
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
        bam=temp(f"01-bowtie/{{sample}}_shard{{shard}}.{bowtie2_human_db_name}.bam"),
        unmapped_R1=temp(
            f"01-bowtie/{{sample}}_shard{{shard}}.without_{bowtie2_human_db_name}.R1.fastq.gz"
        ),
        unmapped_R2=temp(
            f"01-bowtie/{{sample}}_shard{{shard}}.without_{bowtie2_human_db_name}.R2.fastq.gz"
        ),
    log:
        e=f"logs/bowtie2_{{sample}}_shard{{shard}}.{bowtie2_human_db_name}.e",
        o=f"logs/bowtie2_{{sample}}_shard{{shard}}.{bowtie2_human_db_name}.o",


use rule tally_depletion from hostdeplete as bt_tally_depletion with:
    input:
        bam01=f"01-bowtie/{{sample}}.{bowtie2_human_db_name}.bam",
        bam02=f"02-snap/{{sample}}.{snap_human_db_name}.bam",
        bam04=f"04-bowtie/{{sample}}.{bowtie2_mouse_db_name}.bam",
        bam05=f"05-snap/{{sample}}.{snap_mouse_db_name}.bam",
    output:
        table=temp("hostdepleted/{sample}_hostdepletion.stats.tmp"),


rule cat_depletion_stats:
    input:
        table=[
            f"hostdepleted/{{sample}}_shard{shard}_hostdepletion.stats.tmp"
            for shard in SHARDS
        ],
    output:
        table="reports/{sample}_hostdepletion.stats",
    shell:
        """
        head -n 1 {input.table[0]} > {output.table}
        for i in {input.table}
        do
            tail  -n+2 $i >> {output.table}
        done
        """


rule merge_fastq_pair:
    input:
        R1=[f"06-nohuman-nomouse/{{sample}}_shard{shard}.R1.fastq" for shard in SHARDS],
        R2=[f"06-nohuman-nomouse/{{sample}}_shard{shard}.R2.fastq" for shard in SHARDS],
    output:
        R1="hostdepleted/{sample}_R1.fastq.gz",
        R2="hostdepleted/{sample}_R2.fastq.gz",
    container:
        config["docker_cutadapt"]
    threads: 16
    log:
        e="logs/merge_compress_{sample}.e",
    shell:
        """
        cat {input.R1} | pigz -p {threads} -9 > {output.R1} 2> {log.e}
        cat {input.R2} | pigz -p {threads} -9 > {output.R2} 2>> {log.e}
        """


rule aligned_host_reads_to_fastq:
    input:
        bam="{id}/{sample}_shard{shard}.{db}.bam",
    output:
        bam=temp("host/{id}/{sample}_shard{shard}.{db}.bam"),
        R1=temp("host/{id}/{sample}_shard{shard}.{db}.R1.fq"),
        R2=temp("host/{id}/{sample}_shard{shard}.{db}.R2.fq"),
    threads: 8
    resources:
        runtime=8 * 60,
    container:
        config["docker_bowtie2"]
    shell:
        """
        # get the aligned reads
        samtools view -f 2 -F 512 -b -o {output.bam} {input.bam}
        # convert to fastq
        bamToFastq -i {output.bam} -fq {output.R1} -fq2 {output.R2}
        """


rule make_combined_host_reads_fastq:
    """ Get all the host-associated reads and convert back to fastqs
    """
    input:
        R1=expand(
            expand(
                "host/{id}/{{sample}}_shard{{shard}}.{db}.R{{{{readdir}}}}.fq",
                zip,
                id=["01-bowtie", "02-snap", "04-bowtie", "05-snap"],
                db=[
                    bowtie2_human_db_name,
                    snap_human_db_name,
                    bowtie2_mouse_db_name,
                    snap_mouse_db_name,
                ],
            ),
            shard=SHARDS,
            sample=config["sample"],
        ),
    output:
        R1="host/{sample}_all_host_reads_R{readdir}.fastq.gz",
    threads: 8
    resources:
        runtime=8 * 60,
    container:
        config["docker_bowtie2"]
    shell:
        """
        cat {input.R1} | pigz -p {threads} -9 >> {output.R1}
        """


rule sortmerna_run:
    # sortmerna has a bug where they default to writing workdir to your
    # home, and to $HOME/<workdir> if you don't provide an absolute path
    input:
        R1="trimmed/{sample}_shard{shard}_trim_R1.fastq.gz",
        R2="trimmed/{sample}_shard{shard}_trim_R2.fastq.gz",
        db=config["blast_16s_db_nsq"].replace(".nsq", ".fna"),
    output:
        blast=temp("sortmerna/{sample}_shard{shard}_sortmerna.blast.gz"),
        stats=temp("sortmerna/{sample}_shard{shard}_sortmerna.log"),
        work=temp(directory("sortmerna/tmp_{sample}_shard{shard}/")),
    params:
        aligned_prefix=lambda wildcards, output: output["blast"].replace(
            ".blast.gz", ""
        ),
    resources:
        mem_mb=16000,
    threads: 16
    message:
        "Quantify rRNA for qPCR normalization and 16S comparison"
    container:
        config["docker_sortmerna"]
    conda:
        "../envs/sortmerna.yaml"
    log:
        e="logs/sortmerna_{sample}_shard{shard}.e",
        o="logs/sortmerna_{sample}_shard{shard}.o",
    shell:
        """
        sortmerna \
            --ref {input.db} \
            --reads {input.R1} \
            --reads {input.R2} \
            --aligned {params.aligned_prefix} \
            --workdir $PWD/{output.work} \
            -e 0.1 \
            --blast '1 cigar qcov qstrand' \
            --threads {threads}  -v \
            > {log.o} 2> {log.e}
        """


rule merge_sortmerna_blast:
    input:
        blast=expand(
            "sortmerna/{sample}_shard{shard}_sortmerna.blast.gz",
            sample=config["sample"],
            shard=SHARDS,
        ),
    output:
        blast=f"sortmerna/{config['sample']}_sortmerna.blast.gz",
    log:
        e=f"logs/merge_sortmerna_blast_{config['sample']}.e",
    shell:
        """
        cat {input.blast} > {output.blast} 2>> {log.e}
        """


rule merge_logs_for_multiqc:
    input:
        sortmernas=expand(
            "sortmerna/{{sample}}_shard{shard}_sortmerna.log", shard=SHARDS
        ),
        knead="reports/{sample}_hostdepletion.stats",
        bbtrim=expand("trimmed/{{sample}}_shard{shard}_trimmingAQ.txt", shard=SHARDS),
    output:
        sortmerna_report="reports/{sample}_sortmerna.stats",
        sortmerna_log="reports/{sample}_sortmerna.merged.log",
        knead="reports/{sample}_hostdeplete.stats.summary_mqc.tsv",
        bbtrim="reports/{sample}_trimmingAQ_summary.txt",
    params:
        sample_name=lambda wc: wc.sample,
    container:
        "docker://pandas/pandas:pip-all"
    log:
        e="logs/merge_logs_{sample}.e",
        o="logs/merge_logs_{sample}.o",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/merge_logs.py"


rule merge_primary_host_align_bams:
    """ Sort the BAMs in a temp folder, then merge.
    The bowtie alignment has all the reads so we use
    -f 2 to get reads in proper pair and
    -F 512 to remove low-quality reads
    """
    input:
        bams=expand(
            "01-bowtie/{sample}_shard{shard}.{db}.bam",
            sample=config["sample"],
            shard=SHARDS,
            db=bowtie2_human_db_name,
        ),
    output:
        bam=f"host/{{sample}}.{bowtie2_human_db_name}.bam",
        bai=f"host/{{sample}}.{bowtie2_human_db_name}.bam.bai",
    container:
        config["docker_bowtie2"]
    shell:
        """
        mkdir -p tmp_{wildcards.sample}_merge
        for i in {input.bams}
        do
            thisbasename=$(basename $i)
            thisbase=${{thisbasename%.*}}
            echo "getting passing mapped reads $i"
            samtools view  -f 2 -F 512 -b -o tmp_{wildcards.sample}_merge/${{thisbase}}.bam $i
            echo "sorting $i"
            samtools sort -o \
                tmp_{wildcards.sample}_merge/${{thisbase}}.sort.bam \
                tmp_{wildcards.sample}_merge/${{thisbase}}.bam
        done
        echo "merging"
        samtools merge --threads {threads} -o {output.bam} \
            tmp_{wildcards.sample}_merge/*.sort.bam
        samtools index {output.bam}
        rm -r tmp_{wildcards.sample}_merge/
        """
