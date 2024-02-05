import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


validate(config, os.path.join(str(workflow.basedir), "../../config/config.schema.yaml"))


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])


onstart:
    # see https://github.com/snakemake/snakemake/issues/2663
    # for why this isn't in envvars anymore
    try:
        print(os.environ["SNAKEMAKE_PROFILE"])
    except KeyError:
        raise ValueError("Warning: SNAKEMAKE_PROFILE must be set")
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)
    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


kraken_outputs = f"kraken2/{config['sample']}_kraken2.report"
kraken_unclassified_outputs = expand(
    "kraken2/{sample}_kraken2_unclassified_{readdir}.fastq.gz",
    sample=config["sample"],
    readdir=[1, 2],
)
brackenreport = expand(
    "kraken2/{sample}_kraken2.bracken.{taxlevel}.report",
    sample=config["sample"],
    taxlevel=["G", "S"],
)
krona = expand(
    "reports/{sample}_kraken2.bracken.S.report.krona.html", sample=config["sample"]
)
phanta_db_names = ["phanta_unmasked_db", "phanta_masked_db"]
phanta_outputs = expand(
    f'phanta_{config["sample"]}/{{db}}/results/final_merged_outputs/counts.txt',
    db=phanta_db_names,
)

phanta_extra_outputs = expand(
    f'phanta_{config["sample"]}/{{db}}/results/{{result}}',
    db=phanta_db_names,
    result=[
        "counts_by_host.tsv",
        "lifestyle_stats.txt",
        "integrated_prophages_detection_results.txt",
    ],
)

all_outputs = [kraken_outputs, kraken_unclassified_outputs, brackenreport]
if not config["skip_phanta"]:
    all_outputs.extend(phanta_outputs)
    all_outputs.extend(phanta_extra_outputs)


rule all:
    input:
        all_outputs,


#
# Utils Module
if len(config["R1"]) == 1:
    input_R1 = config["R1"]
    input_R2 = config["R2"]
else:
    input_R1 = f"concatenated/{config['sample']}_R1.fastq.gz"
    input_R2 = f"concatenated/{config['sample']}_R2.fastq.gz"


module utils:
    snakefile:
        "utils.smk"
    config:
        config
    skip_validation:
        True


use rule concat_R1_R2 from utils as utils_concat_R1_R2 with:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp("concatenated/{sample}_R1.fastq.gz"),
        R2=temp("concatenated/{sample}_R2.fastq.gz"),
    log:
        e="logs/concat_r1_r2_{sample}.e",


rule kraken_standard_run:
    """ Profile the microbiome with Kraken2.

    We run this on the raw reads rather than the host-depleted trimmed reads because
    - its nice to have validation of the host detection percentages
    - we want uniform read lengths for bracken
    - its what the authors seem to recommend, although
      Nick has asked for confirmation: https://github.com/DerrickWood/kraken2/issues/646

    We set the confidence threshold after reading
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full

    We wanted to
      - minimize the L1 distance to the known community composition
      - maximize the alpha diversity similarity to the actual diversity

    We don't care as much about reads assigned so long as the assignments we
    do get are high quality and representative
    """
    input:
        R1=input_R1,
        R2=input_R2,
        db=config["kraken2_db"],
    output:
        out="kraken2/{sample}_kraken2.out",
        unclass_R1=temp("kraken2/{sample}_kraken2_unclassified_1.fastq"),
        unclass_R2=temp("kraken2/{sample}_kraken2_unclassified_2.fastq"),
        report="kraken2/{sample}_kraken2.report",
    params:
        inpstr=lambda wc, input: f"--paired {input.R1} {input.R2}" if hasattr(input, "R2")  else input.R1,
        unclass_template=lambda wildcards, output: output["unclass_R1"].replace(
            "_1.fastq", "#.fastq"
        ),
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/kraken_{sample}.log",
    resources:
        # use 8GB mem if using a minikraken db to make testing easier
        mem_mb=lambda wildcards, attempt: 8 * 1024
        if "mini" in config["kraken2_db"]
        else (64 * 1024 * attempt),
    threads: 16
    shell:
        """
        kraken2 \
            --threads {threads} \
            --use-names \
            --confidence 0.2 \
            --unclassified-out {params.unclass_template} \
            --db {input.db} \
            --report {output.report} \
            {params.inpstr} \
            > {output.out} 2> {log.e}

        # some (mock) datasets are perfect but we still need these files
        if [ ! -f "{output.unclass_R1}" ]
        then
            touch {output.unclass_R1}
            touch {output.unclass_R2}
        fi
        """


checkpoint get_read_len:
    # get the max read length needed by bracken, and name a file in its honor
    input:
        R1=config["R1"],
    output:
        readlen_dir=directory("kraken2/readlens/"),
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/get_read_len.log",
    resources:
        mem_mb=lambda wildcards, attempt: 8 * 1024,
    threads: 1
    shell:
        """
        # thanks https://www.biostars.org/p/72433/ for the readlen calculation with awk
        # however, there is a strange error with awk here where it returns the proper result but returns an error code of 1.
        # we allow pipefails here to account for that
        set +o pipefail
        READLEN=$(zcat {input.R1} | head -n 100 | awk "NR%4 == 2 {{lengths[length(\$0)]++}} END {{for (l in lengths) {{print l}}}}")
        READLENMAX=$(for i in $READLEN; do echo $i ; done  | awk  'NR==1{{max=$1}} $1>max{{max=$1}} END{{print max+0}}')
        mkdir -p kraken2/readlens
        touch kraken2/readlens/len${{READLENMAX}}
        set -o pipefail
        """


rule bracken:
    """Run bracken using the probabilities based on the closest read length pre-computed.
    Note that this depreciates the need for the get_db_mers function,
    but we are leaving it here in case we decide that pre-computed
    probabilities are inadequate
    """
    input:
        #db_mers=get_dbs_needed,
        readlen_file=parse_read_lengths,
        report="kraken2/{sample}_kraken2.report",
        db=config["kraken2_db"],
    output:
        report="kraken2/{sample}_kraken2.bracken.{taxlevel}.report",
        out="kraken2/{sample}_kraken2.bracken.{taxlevel}.out",
    params:
        threshold=0,
    container:
        config["docker_kraken"]
    conda:
        "../envs/kraken2.yaml"
    log:
        e="logs/bracken_{sample}.{taxlevel}.log",
    resources:
        mem_mb=8 * 1024,
    threads: 4
    shell:
        """
        # trim off the approximate and the len part of the read len signal
        READLEN=$(basename {input.readlen_file} | sed 's|len||' | sed  's|approx||')
        bracken -d {input.db} -i {input.report} -o {output.out} -w {output.report} -r $READLEN -l {wildcards.taxlevel} -t {params.threshold} 2> {log.e}
        find .
        """


rule kraken2krona:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        report="kraken2/{sample}_kraken2.bracken.S.report",
    output:
        report="kraken2/{sample}_kraken2.bracken.S.report.krona",
    conda:
        "../envs/kraken2.yaml"
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken2krona_{sample}.e",
        o="logs/kraken2krona_{sample}.o",
    shell:
        """
        kreport2krona.py --report-file {input.report} \
        --output {output.report} > {log.o} 2> {log.e}
        """


rule krona:
    input:
        report="kraken2/{sample}_kraken2.bracken.S.report.krona",
    output:
        report="reports/{sample}_kraken2.bracken.S.report.krona.html",
    container:
        config["docker_krona"]
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/krona_{sample}.e",
    shell:
        """
        ktImportText {input.report} -o {output.report}
        """


rule compress_kraken_unclassified:
    input:
        R1="kraken2/{sample}_kraken2_unclassified_1.fastq",
        R2="kraken2/{sample}_kraken2_unclassified_2.fastq",
    output:
        R1="kraken2/{sample}_kraken2_unclassified_1.fastq.gz",
        R2="kraken2/{sample}_kraken2_unclassified_2.fastq.gz",
    conda:
        "../envs/pigz.yaml"
    # this container is used in the 16S pipeline
    container:
        config["docker_cutadapt"]
    threads: 8
    log:
        e="logs/kraken_compressed_unclassified_{sample}.e",
    shell:
        """
        cat {input.R1} | pigz -p {threads} -9 > {output.R1} 2> {log.e}
        cat {input.R2} | pigz -p {threads} -9 > {output.R2} 2>> {log.e}
        """


rule kraken_merge_shards:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        reports=[
            f"kraken2/{config['sample']}_shard{shard}_kraken2.report"
            for shard in SHARDS
        ],
    output:
        out="kraken2/{sample}_kraken2_merged.report",
    conda:
        "../envs/kraken2.yaml"
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken_merge_{sample}.e",
        o="logs/kraken_merge_{sample}.o",
    shell:
        """
        combine_kreports.py \
            --only-combined \
            --no-headers \
            -r {input.reports} \
            -o {output.out} \
            > {log.o} 2> {log.e}
        """


rule make_phanta_manifest:
    """ {params.thisdir} is needed to list the absolute path when using concatenated libraries.
    PHANTA needs the full path, and Input_R1 is relative the relative paths if we are using the concatenated lanes
    """
    input:
        R1=input_R1,
        R2=input_R2,
    output:
        manifest=temp("phanta_inputs.tsv"),
    params:
        thisdir=os.getcwd() if isinstance(input_R1, str) else "",
        sample=config["sample"],
    shell:
        """echo -e "{params.sample}\t{params.thisdir}/{input.R1}\t{params.thisdir}/{input.R2}" > {output.manifest}"""


def get_singularity_args(wildcards):
    with open(os.path.join(os.environ["SNAKEMAKE_PROFILE"], "config.yaml"), "r") as f:
        return yaml.safe_load(f)["singularity-args"]


rule phanta:
    # we need the results dir because if we try to initiate two snakemake workflows from the same --directory we get LockErrors
    input:
        R1=input_R1,
        R2=input_R2,
        manifest=rules.make_phanta_manifest.output.manifest,
        readlen_file=parse_read_lengths,
    output:
        counts="phanta_{sample}/{db}/results/final_merged_outputs/counts.txt",
        bracken_failed="phanta_{sample}/{db}/results/classification/samples_that_failed_bracken.txt",
        bracken_report="phanta_{sample}/{db}/results/classification/{sample}.krak.report_bracken_species.filtered",
        bracken_species_abundance="phanta_{sample}/{db}/results/classification/{sample}.krak.report.filtered.bracken.scaled",
        relative_read_abundance="phanta_{sample}/{db}/results/final_merged_outputs/relative_read_abundance.txt",
        relative_taxonomic_abundance="phanta_{sample}/{db}/results/final_merged_outputs/relative_taxonomic_abundance.txt",
        total_reads="phanta_{sample}/{db}/results/final_merged_outputs/total_reads.tsv",
        single_end=temp(
            directory("phanta_{sample}/{db}/results/classification/single_end")
        ),
    container:
        config["docker_phanta"]
    params:
        dbpath=lambda wc: config[wc.db],
        sing_args=lambda wc: get_singularity_args(wc),
        cov_thresh_viral=config["cov_thresh_viral"],
        cov_thresh_arc=config["cov_thresh_arc"],
        cov_thresh_bacterial=config["cov_thresh_bacterial"],
        cov_thresh_euk=config["cov_thresh_euk"],
        minimizer_thresh_viral=config["minimizer_thresh_viral"],
        minimizer_thresh_arc=config["minimizer_thresh_arc"],
        minimizer_thresh_bacterial=config["minimizer_thresh_bacterial"],
        minimizer_thresh_euk=config["minimizer_thresh_euk"],
        single_end_krak=config["single_end_krak"],
    resources:
        mem_mb=lambda wildcards, attempt: 58 * 1024 * attempt,
    threads: 16
    # we have to do the dummy profile to keep the job from inheriting SNAKEMAKE_PROFILE;
    #   we want this job to be submitted locally.  Otherwise, we have this unpleasant
    # situation where we need the docker container for its snakefile, but also
    # need the workflow itself to execute within the same container.  Cue mount point conflicts,
    # the inability to test, etc.
    # The workaround is running it without a --cluster directive so it runs 'locally',
    # which is actually on the node this rule got submitted to. While there is a slight
    # inefficiency, phanta doesn't benefit enough from distributing its workflow rules to
    # warrant going the alternative route of adding it to this repo as a submodule,
    # importing the rules, etc.
    shell:
        """
        READLEN=$(basename {input.readlen_file} | sed 's|len||' | sed  's|approx||')
        echo "cores: {threads}" > config.yaml && \
        export SNAKEMAKE_PROFILE="" && \
        snakemake  --profile $PWD --notemp \
        --singularity-args "{params.sing_args}" \
        --snakefile /home/mambauser/phanta/Snakefile \
        --configfile /home/mambauser/phanta/config.yaml \
        --directory phanta_{wildcards.sample}/{wildcards.db}/ \
        --config \
        outdir="results/" \
        read_length=$READLEN \
        cov_thresh_viral={params.cov_thresh_viral} \
        cov_thresh_bacterial={params.cov_thresh_bacterial} \
        cov_thresh_euk={params.cov_thresh_euk} \
        cov_thresh_arc={params.cov_thresh_arc} \
        minimizer_thresh_viral={params.minimizer_thresh_viral} \
        minimizer_thresh_bacterial={params.minimizer_thresh_bacterial} \
        minimizer_thresh_euk={params.minimizer_thresh_euk} \
        minimizer_thresh_arc={params.minimizer_thresh_arc} \
        single_end_krak={params.single_end_krak} \
        class_mem_mb={resources.mem_mb} \
        database={params.dbpath} \
        pipeline_directory=/home/mambauser/phanta/ \
        sample_file=$PWD/{input.manifest}
        """


rule phanta_postprocess:
    input:
        counts="phanta_{sample}/{db}/results/final_merged_outputs/counts.txt",
        tax_relab="phanta_{sample}/{db}/results/final_merged_outputs/relative_taxonomic_abundance.txt",
        manifest=rules.make_phanta_manifest.output.manifest,
        single_end="phanta_{sample}/{db}/results/classification/single_end",
    output:
        byhost="phanta_{sample}/{db}/results/counts_by_host.tsv",
        lifestyle="phanta_{sample}/{db}/results/lifestyle_stats.txt",
        prophage="phanta_{sample}/{db}/results/integrated_prophages_detection_results.txt",
    container:
        config["docker_phanta"]
    params:
        dbpath=lambda wc: config[wc.db],
        bacphlip_thresh=0.5,
        classif_dir=lambda wc, input: os.path.dirname(input.single_end),
        out_dir=lambda wc, output: os.path.dirname(output.prophage),
    resources:
        mem_mb=8 * 1024,
    threads: 1
    shell:
        """
        python /home/mambauser/phanta/post_pipeline_scripts/collapse_viral_abundances_by_host.py {input.counts} {params.dbpath}/host_prediction_to_genus.tsv {output.byhost}
        # check if any viruses detected; file has header
        if [ $(wc -l < "{output.byhost}") -gt 1 ]
        then
            Rscript /home/mambauser/phanta/post_pipeline_scripts/calculate_lifestyle_stats/lifestyle_stats.R \
              {params.bacphlip_thresh} {params.dbpath}/species_name_to_vir_score.txt \
              {input.counts} \
              {input.tax_relab} \
              {output.lifestyle}
            python /home/mambauser/phanta/post_pipeline_scripts/integrated_prophages_detector.py \
              {input.manifest} \
              {params.dbpath} \
              {params.classif_dir} \
              {params.out_dir}
        else
            touch {output.lifestyle}
            touch {output.prophage}
        fi
        """
