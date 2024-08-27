import os
import math
import traceback
import numpy as np
import pandas as pd
from contextlib import redirect_stderr

include: "common.smk"

envvars:
    "TMPDIR",

# ---------------------------------------------------------------------------------------------------------------------

SORTMERNA_PERCENTAGE_THRESHOLD = 20
BBDUK_TRIM_PERCENTAGE_THRESHOLD = 10
NOTIFICATION_EMAIL = "zzPDL_SKI_microbiome@mskcc.org"

# ---------------------------------------------------------------------------------------------------------------------

def get_manifest_notes(manifest_path):
    """Collate any notes about file modification into a note for report."""
    full_manifest = pd.read_csv(manifest_path, sep="\t",
    header = 0)
    
    mod_message = ['']
    for row in full_manifest[full_manifest.notes.notna()].iterrows():
        mod_message.append(f"File {row.file_type} was modified for experiment {row['experiment.identifier']}"
            f" of sample {config['sample']}.  Reason: {row.notes}.\n")
    return ''.join(mod_message)

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

def load_manifest(manifest_path, validate = True):
    """Load in the manifest stored and pivot so each row is an experiment, and columns are files."""
    manifest = pd.read_csv(manifest_path, sep="\t", header = 0).drop(
            columns = ['notes']  # We drop notes from the manifest, but include in the report. 
        ).pivot(
            index = "experiment.identifier", columns = "file_type", values = "file_path"
        )
    manifest.to_csv("test.csv")
    print(manifest)
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

MANIFEST = load_manifest(config["manifest"])
EXPERIMENTS = get_experiments_from_manifest(MANIFEST)
LOG_PREFIX = "logs/qc"

onstart:
    with open("config_used.yaml", "w") as f:
        yaml.dump(config, f)

localrules:
    all,

def get_all_outputs(wildcards):
    return [
        f"qc/{config['sample']}_qc_merged_R1.fq",
        f"qc/{config['sample']}_qc_merged_R2.fq",
        f"qc/{config['sample']}_qc_stats.txt"
    ]

rule all:
    input:
        get_all_outputs,

rule check_sortmerna_rep:
    input: 
        smr_logs = MANIFEST["sortmerna_report"].tolist(),
    output:
        sortmerna_report = f"qc/{config['sample']}_pooled_sortme_rna_report.txt",
    log:
        e=f"{LOG_PREFIX}/qc_sortmerna_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_sortmerna_{config['sample']}.o",
    params:
        experiment_ids = EXPERIMENTS,
        threshold=SORTMERNA_PERCENTAGE_THRESHOLD,
    script:
        "../scripts/qc/check_sortmerna_report.py"

#rule check_bbduk_trimming:
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


rule merge_and_check_all_reports:
    input:
        sortmerna_report = f"qc/{config['sample']}_pooled_sortme_rna_report.txt",
#        bbduk_report = f"qc/{config['sample']}_pooled_bbduk_trim_report.txt",
    output:
        qc_stats = f"qc/{config['sample']}_qc_stats.txt"
    log:
        e=f"{LOG_PREFIX}/qc_stats_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_stats_{config['sample']}.o",
    params:
        notification_email=NOTIFICATION_EMAIL,
    shell:
        """
        cat {input.sortmerna_report} > {output.qc_stats} 2> {log.e}
        if grep -q -i 'error' ${output.qc_stats}; then
            echo "Error found in qc for sample {config['sample']}. \nPlease check out the qc_stats file for details." | mail -s "Potential IGO rsync error" {params.notification_email}
        fi
        """

rule merge_all_experiments:
    input:
        fqs = get_hostdepleted_fastqs,
    output:
        R1 = f"qc/{config['sample']}_qc_merged_R1.fq",
        R2 = f"qc/{config['sample']}_qc_merged_R2.fq",
    log:
        e=f"{LOG_PREFIX}/qc_merge_{config['sample']}.e",
        o=f"{LOG_PREFIX}/qc_merge_{config['sample']}.o",
    shell:
        """
        cat {input.fqs[0]} > {output.R1} 2> {log.e}
        cat {input.fqs[1]} > {output.R2} 2>> {log.e}
        """

