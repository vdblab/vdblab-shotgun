import os
import json
import yaml

from pathlib import Path
from glob import glob


configfile: os.path.join(workflow.basedir, "../config/config.yaml")


envvars:
    "TMPDIR",


tmpdir = Path(os.environ["TMPDIR"])


sample_pkls = []
with open(config["manifest"], "r") as manifest_f:
    for line in manifest_f:
        line = line.strip()
        if line:
            sample_pkls.append(line)
logger.info(
    "Running the pipeline on the following Sample pkls:\n - %s"
    % "\n  - ".join(sample_names)
)


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)


if not os.path.exists("logs"):
    os.makedirs("logs")

humann_container = "docker://nickp60/humann:3.1.1"
# ==================================================================================================================================================


localrules:
    all,


if "targets" not in config:
    targets = [
        "s__Blautia_coccoides",
        "s__Blautia_wexlerae",
        "s__Blautia_producta",
        "s__Blautia_obeum",
        "s__Enterococcus_faecalis",
        "s__Enterococcus_faecium",
        "s__Erysipelatoclostridium_ramosum",
    ]
else:
    targets = config["targets"]

for t in targets:
    os.makedirs("markers_" + t, exist_ok=True)


rule all:
    input:
        expand(
            "strainphlan/strainphlan_{sp}_output/"
            "RAxML_bestTree.{sp}.StrainPhlAn3.tre",
            sp=targets,
        ),


rule extract_sp_markers:
    input:
        outdir="markers_{sp}/",
        db=config["metaphlan_pkl"],
    output:
        fasta="markers_{sp}/{sp}.fna",
    params:
        sp="{sp}",
    container:
        humann_container
    conda:
        "envs/metaphlan.yaml"
    shell:
        """
        extract_markers.py --clade {params.sp} --database {input.db} --output_dir {input.outdir}
        if [ ! -s "{output.fasta}" ]
        then
            echo "output fasta empty; double check path to metaphlan database pkl file"
            exit 1
        fi
        """


rule strainphlan_run:
    input:
        sp_markers=rules.extract_sp_markers.output.fasta,  #,"markers_{sp}/markers_{sp}.fna",
        sample_pkls=sample_pkls,
    output:
        "strainphlan/strainphlan_{sp}_output/RAxML_bestTree.{sp}.StrainPhlAn3.tre",
    params:
        chocophlan_db=config["choco_db"],
        marker_in_n_samples=config["marker_in_n_samples"],
        sp="{sp}",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        humann_container
    conda:
        "envs/metaphlan.yaml"
    threads: 16
    resources:
        mem_mb=16000,
    shell:
        """
        strainphlan \
            --database {params.chocophlan_db} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            --samples {input.sample_pkls} \
            --clade_markers {input.sp_markers} \
            --output_dir {params.outdir} \
            --clade {params.sp} \
            --phylophlan_mode fast \
            --nproc {threads}
        """
