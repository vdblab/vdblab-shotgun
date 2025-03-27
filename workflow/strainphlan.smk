import os
import sys
import json
import yaml

from pathlib import Path
from glob import glob


configfile: os.path.join(workflow.basedir, "../config/config.yaml")


if not config["strainphlan_markers_dir"]:
    sys.exit(
        "config 'strainphlan_markers_dir' must be set!  This is where we look for existing marker for a sample, or deposit them if the markers do not exist"
    )
if not os.path.exists(config["strainphlan_markers_dir"]):
    sys.exit(f"{config['strainphlan_markers_dir']} does not exist!")
else:
    # create the samples and species subdirectories if needed
    os.makedirs(
        os.path.join(config["strainphlan_markers_dir"], "samples"), exist_ok=True
    )
    os.makedirs(
        os.path.join(config["strainphlan_markers_dir"], "species"), exist_ok=True
    )

if not config["sams"]:
    sys.exit("config 'sams' is a required input for this workflow")


envvars:
    "TMPDIR",


ref_string = (
    "" if "references" not in config else "-r " + " ".join(config["references"])
)
tmpdir = Path(os.environ["TMPDIR"])


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)


if not os.path.exists("logs"):
    os.makedirs("logs")

# ==================================================================================================================================================


localrules:
    all,


# see https://opendata.lifebit.ai/table/sgb
if "targets" not in config:
    targets = [
        # "s__Blautia_coccoides",
        # "s__Blautia_wexlerae",
        # "s__Blautia_producta",
        # "t__SGB4844", #Blautia_obeum
        "t__SGB7962",  # Enterococcus_faecalis",
        "t__SGB7967",  # Enterococcus_faecium",
        # "t__SGB7968", #also Enterococcus_faecium",
        "t__SGB6746",  # Erysipelatoclostridium_sp",
    ]
else:
    targets = config["targets"]

for t in targets:
    os.makedirs("markers_" + t, exist_ok=True)

SAMPLES = {}
with open(config["sams"], "r") as samin:
    for line in samin:
        SAMPLES[os.path.basename(line.strip()).replace(".sam.bz2", "")] = line.strip()


def get_sam_path(wc):
    return SAMPLES[wc.sample]


rule all:
    input:
        expand(
            os.path.join(
                config["strainphlan_markers_dir"], "samples", "{sample}.json.bz2"
            ),
            sample=SAMPLES.keys(),
        ),
        expand(
            os.path.join(config["strainphlan_markers_dir"], "species", "{sp}.fna"),
            sp=targets,
        ),
        expand(
            "strainphlan/strainphlan_{sp}_output/"
            "RAxML_bestTree.{sp}.StrainPhlAn4.tre",
            sp=targets,
        ),


# Run sample2markers for strainphlan2
rule sample2markers_run:
    input:
        inf=get_sam_path,
        db=os.path.join(
            config["metaphlan_db"],
            os.path.basename(os.path.dirname(config["metaphlan_db"] + "/")) + ".pkl",
        ),
    output:
        os.path.join(config["strainphlan_markers_dir"], "samples", "{sample}.json.bz2"),
    threads: 8
    container:
        config["docker_biobakery"]
    resources:
        mem_mb=lambda wildcards, attempt, input: attempt
        * 1024
        * max(input.inf.size // 1000000000, 1)
        * 10,
        runtime=24 * 60,
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    log:
        e="logs/sample2markers_{sample}.e",
    shell:
        """
        sample2markers.py \
            --database {input.db} \
            -i {input.inf} \
            -o {params.outdir} \
            --database {input.db} \
            --nprocs {threads} \
            2> {log.e}
        """


rule extract_sp_markers:
    input:
        db=os.path.join(
            config["metaphlan_db"],
            os.path.basename(os.path.dirname(config["metaphlan_db"])) + ".pkl",
        ),
    output:
        fasta=os.path.join(config["strainphlan_markers_dir"], "species", "{sp}.fna"),
    resources:
        runtime=3 * 60,
        mem_mb=32 * 1024,
    params:
        sp="{sp}",
        outdir=lambda wc, output: os.path.dirname(output.fasta),
    container:
        config["docker_biobakery"]
    conda:
        "envs/metaphlan.yaml"
    shell:
        """
        extract_markers.py --clade {params.sp} --database {input.db} --output_dir {params.outdir}
        if [ ! -s "{output.fasta}" ]
        then
            echo "output fasta empty; double check path to metaphlan database pkl file"
            exit 1
        fi
    """


rule strainphlan_run:
    input:
        sp_markers=rules.extract_sp_markers.output.fasta,
        sample_pkls=expand(
            os.path.join(
                config["strainphlan_markers_dir"], "samples", "{sample}.json.bz2"
            ),
            sample=SAMPLES.keys(),
        ),
    output:
        "strainphlan/strainphlan_{sp}_output/RAxML_bestTree.{sp}.StrainPhlAn4.tre",
    params:
        metaphlan_pkl=os.path.join(
            config["metaphlan_db"],
            os.path.basename(os.path.dirname(config["metaphlan_db"] + "/")) + ".pkl",
        ),
        marker_in_n_samples=config["marker_in_n_samples"],
        sp="{sp}",
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        references_str=ref_string,
    container:
        config["docker_biobakery"]
    conda:
        "envs/metaphlan.yaml"
    threads: 16
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        strainphlan \
            --database {params.metaphlan_pkl} \
            --marker_in_n_samples {params.marker_in_n_samples} \
            --samples {input.sample_pkls} \
        {params.references_str} \
            --clade_markers {input.sp_markers} \
            --output_dir {params.outdir} \
            --clade {params.sp} \
            --phylophlan_mode fast \
            --nproc {threads}

        """
