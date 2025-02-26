import os
import math
import traceback
import numpy as np
import pandas as pd
from contextlib import redirect_stderr


wildcard_constraints:
    depth="\d+",


include: "common.smk"


module downsample:
    snakefile:
        "downsample.smk"
    config:
        config
    skip_validation:
        True


module kraken:
    snakefile:
        "kraken.smk"
    config:
        config
    skip_validation:
        True


envvars:
    "TMPDIR",


# ---------------------------------------------------------------------------------------------------------------------

SORTMERNA_PERCENTAGE_THRESHOLD = 20
BBDUK_TRIM_PERCENTAGE_THRESHOLD = 10

# ---------------------------------------------------------------------------------------------------------------------


def get_manifest_notes(manifest_path):
    """Collate any notes about file modification into a note for report."""
    full_manifest = pd.read_csv(manifest_path, sep="\t", header=0)

    mod_message = [""]
    for row in full_manifest[full_manifest.notes.notna()].iterrows():
        mod_message.append(
            f"File {row.file_type} was modified for experiment {row['experiment.identifier']}"
            f" of sample {config['sample']}.  Reason: {row.notes}.\n"
        )
    return "".join(mod_message)


def validate_manifest(manifest):
    """Run some sanity checks on the manifest"""

    # To start with we can check if we have any NAs in the manifest - which would imply not all
    # file_types provided for each experiment.
    if manifest.isnull().values.any():
        raise ValueError(
            "Some values in the manifest were null.  Please modify and resubmit. "
            "Perhaps typo in file_type? See problematic experiments here: \n\n"
            f"{manifest[manifest.isna().any(axis=1)]}"
        )


def load_manifest(manifest_path, validate=True):
    """Load in the manifest stored and pivot so each row is an experiment, and columns are files."""
    manifest = (
        pd.read_csv(manifest_path, sep="\t", header=0)
        .drop(
            columns=[
                "notes"
            ]  # We drop notes from the manifest, but include in the report.
        )
        .pivot(index="experiment.identifier", columns="file_type", values="file_path")
    )
    if validate:
        validate_manifest(manifest)
    return manifest


def get_experiments_from_manifest(manifest):
    return [e for e in manifest.index.tolist()]


def get_hostdepleted_fastqs(wildcards):
    """If paired will return [[R1 file list][R2 file list]], else [[R1 file list]]"""
    fqs = []
    for d in get_readdirs():
        fqs.append(MANIFEST[f"fq{d}"].tolist())
    return fqs


def get_host_fastqs(wildcards):
    """If paired will return [[R1 file list][R2 file list]], else [[R1 file list]]"""
    fqs = []
    for d in get_readdirs():
        fqs.append(MANIFEST[f"host_fq{d}"].tolist())
    return fqs


# ---------------------------------------------------------------------------------------------------------------

LOG_PREFIX = "logs/qc"
MANIFEST = pd.DataFrame(columns=["sortmerna_report", "fq1", "fq2"])
EXPERIMENTS = None


onstart:
    MANIFEST = load_manifest(config["manifest"])
    EXPERIMENTS = get_experiments_from_manifest(MANIFEST)
    with open("config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,


def get_all_outputs(wildcards):
    return [
        f"qc/{config['sample']}_qc_merged_R1.fq.gz",
        f"qc/{config['sample']}_qc_merged_R2.fq.gz",
        f"qc/{config['sample']}_qc_stats.txt",
    ]


rule all:
    input:
        get_all_outputs,
        expand(
            "qc/kraken2/ds{depth}_{sample}_kraken2.report.krona.html",
            sample=config["sample"],
            depth=10000,
        ),
        expand(
            "qc/kraken2/ds{depth}_{sample}_kraken2.braycurtis_G.tab",
            sample=config["sample"],
            depth=10000,
        ),
        expand(
            "qc/kraken2/ds{depth}_{sample}_kraken2.braycurtis_S.tab",
            sample=config["sample"],
            depth=10000,
        ),


rule check_sortmerna_rep:
    input:
        smr_logs=MANIFEST["sortmerna_report"].tolist(),
    output:
        sortmerna_report=f"qc/{config['sample']}_pooled_sortme_rna_report.txt",
    log:
        e=f"{LOG_PREFIX}/qc_sortmerna_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_sortmerna_{config['sample']}.o",
    params:
        experiment_ids=EXPERIMENTS,
        threshold=SORTMERNA_PERCENTAGE_THRESHOLD,
    script:
        "../scripts/qc/check_sortmerna_report.py"


# rule check_bbduk_trimming:
#    input:
#        bbduk_report = list(MANIFEST["bbtrim_report"]),
#    output:
#        bbduk_report = f"qc/{config['sample']}_pooled_bbduk_trim_report.txt",
#    log:
#        e=f"{LOG_PREFIX}/qc_bbduk_{config['sample']}.e",
#        o=f"{LOG_PREFIX}/qc_bbduk_{config['sample']}.o",
#    params:
#        experiment_ids = EXPERIMENTS,
#        threshold=BBDUK_TRIM_PERCENTAGE_THRESHOLD,
#    script:
#        "../scripts/qc/check_bbduk_report.py"


def get_fqs_by_experiment(wildcards):
    paths = {}
    paths["R1"] = MANIFEST[MANIFEST.index == wildcards.exp]["fq1"].tolist()
    paths["R2"] = MANIFEST[MANIFEST.index == wildcards.exp]["fq2"].tolist()
    assert paths[
        "R1"
    ], f"retrieving fq from manifest with experiment {wildcards.exp} has failed!"
    return paths


use rule downsample_fastq from downsample with:
    input:
        unpack(get_fqs_by_experiment),
    output:
        R1=temp("tmp/ds{depth}_{exp}_R1_001.fastq.gz"),
        R2=temp("tmp/ds{depth}_{exp}_R2_001.fastq.gz"),
    params:
        tmpdir=lambda wildcards, output: os.path.dirname(output.R1),
        seed=lambda wildcards: wildcards.rep if "rep" in wildcards else 1,
        return_too_shallow="yes",


use rule kraken_standard_run from kraken as kracken_standard_run_qc with:
    input:
        R1="tmp/ds{depth}_{exp}_R1_001.fastq.gz",
        R2="tmp/ds{depth}_{exp}_R2_001.fastq.gz",
        db=config["kraken2_db"],
    output:
        out=temp("qc/kraken2/ds{depth}_{exp}_kraken2.out"),
        unclass_R1=temp("qc/kraken2/ds{depth}_{exp}_kraken2_unclassified_1.fastq"),
        unclass_R2=temp("qc/kraken2/ds{depth}_{exp}_kraken2_unclassified_2.fastq"),
        kreport="qc/kraken2/ds{depth}_{exp}_kraken2.report",
    log:
        e="logs/kraken_ds{depth}_{exp}.log",


rule kraken2krona:
    """ This is only needed when running the whole workflow's snakefile;
    otherwise we do not assume sample is split. with the --only-combined
    option, krakentools will merge individual reports into a single report
    """
    input:
        kreport="qc/kraken2/ds{depth}_{exp}_kraken2.report",
    output:
        kreport="qc/kraken2/ds{depth}_{exp}_kraken2.report.krona",
    resources:
        mem_mb=2000,
    container:
        config["docker_kraken"]
    threads: 1
    log:
        e="logs/kraken2krona_ds{depth}_{exp}.e",
        o="logs/kraken2krona_ds{depth}_{exp}.o",
    shell:
        """
        kreport2krona.py --report-file {input.kreport} \
        --output {output.kreport} > {log.o} 2> {log.e}
        """


def get_experiment_from_krona_path(path):
    # why regex when split do trick
    x = "_".join(path.split("_")[1:])
    return x.split("_kraken2")[0]


rule merged_krona:
    input:
        kreport=expand(
            "qc/kraken2/ds{{depth}}_{exp}_kraken2.report.krona", exp=EXPERIMENTS
        ),
    output:
        kreport=f"qc/kraken2/ds{{depth}}_{config['sample']}_kraken2.report.krona.html",
    container:
        config["docker_krona"]
    # see https://github.com/marbl/Krona/issues/125
    params:
        input_string=lambda wc, input: " ".join(
            [f"{x},{get_experiment_from_krona_path(x)}" for x in input.kreport]
        ),
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/krona_ds{depth}_merge.e",
    shell:
        """
        ktImportText {params.input_string} -o {output.kreport}
        """


rule beta_diversity:
    input:
        kreports=expand("qc/kraken2/ds{{depth}}_{exp}_kraken2.report", exp=EXPERIMENTS),
    output:
        table=f"qc/kraken2/ds{{depth}}_{config['sample']}_kraken2.braycurtis_G.tab",
        table_S=f"qc/kraken2/ds{{depth}}_{config['sample']}_kraken2.braycurtis_S.tab",
    container:
        config["docker_kraken"]
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e="logs/krakentools_betadiv_ds{depth}.e",
    shell:
        """
        python3 /KrakenTools/DiversityTools/beta_diversity.py  -i {input.kreports} --type kreport --level G  >  {output.table}
        python3 /KrakenTools/DiversityTools/beta_diversity.py  -i {input.kreports} --type kreport --level S  >  {output.table_S}
        """


rule merge_and_check_all_reports:
    input:
        sortmerna_report=f"qc/{config['sample']}_pooled_sortme_rna_report.txt",
    #        bbduk_report = f"qc/{config['sample']}_pooled_bbduk_trim_report.txt",
    output:
        qc_stats=f"qc/{config['sample']}_qc_stats.txt",
    log:
        e=f"{LOG_PREFIX}/qc_stats_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_stats_{config['sample']}.o",
    params:
        sample=config["sample"],
    shell:
        """
        cat {input.sortmerna_report} > {output.qc_stats} 2> {log.e}
        """


rule merge_all_experiments:
    input:
        fqs1=MANIFEST[f"fq1"].tolist(),
        fqs2=MANIFEST[f"fq2"].tolist(),
    output:
        R1=f"qc/{config['sample']}_qc_merged_R1.fq.gz",
        R2=f"qc/{config['sample']}_qc_merged_R2.fq.gz",
    log:
        e=f"{LOG_PREFIX}/qc_merge_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_merge_{config['sample']}.o",
    shell:
        """
        cat {input.fqs1} > {output.R1} 2> {log.e}
        cat {input.fqs2} > {output.R2} 2>> {log.e}
        """
