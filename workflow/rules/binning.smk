import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


configfile: os.path.join(workflow.basedir, "../../config/config.yaml")


envvars:
    "TMPDIR",


validate(config, os.path.join(workflow.basedir, "../../config/config.schema.yaml"))

SHARDS = make_shard_names(config["nshards"])


onstart:
    with open("config_used.yml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,


BINNING_TOOLS = ["concoct", "metabat2", "maxbin2"]
binstats = expand(
    "metawrap/rawbinning_{sample}/{tool}/{tool}_bins/{tool}.done",
    sample=config["sample"],
    tool=BINNING_TOOLS,
)

refined_binstats_each = expand(
    "metawrap/refined_binning_{sample}/{tool}_bins.stats",
    sample=config["sample"],
    tool=BINNING_TOOLS,
)
refined_stats = f'metawrap/refined_binning_{config["sample"]}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats'
stats_mqc = (
    f'metawrap/refined_binning_{config["sample"]}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats_mqc.tsv',
)
contigs = (
    f'metawrap/refined_binning_{config["sample"]}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.contigs',
)


rule all:
    input:
        binstats,
        refined_stats,
        refined_binstats_each,
        stats_mqc,
        contigs,


rule unzip_rename_fastq_for_metawrap:
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp("tmp_{sample}_1.fastq"),
        R2=temp("tmp_{sample}_2.fastq"),
    shell:
        """
        zcat {input.R1} > {output.R1}
        zcat {input.R2} > {output.R2}
        """


rule metawrap_binning:
    """
    """
    input:
        R1="tmp_{sample}_1.fastq",
        R2="tmp_{sample}_2.fastq",
        assembly=config["assembly"],
    output:
        stats="metawrap/rawbinning_{sample}/{tool}/{tool}_bins/{tool}.done",
    params:
        outdir=lambda wc, output: os.path.dirname(output.stats),
    container:
        config["docker_metawrap"]
    threads: 64
    resources:
        mem_mb=32 * 1024,
        runtime="12:00",
    shell:
        """
        metawrap binning -o {params.outdir} -t {threads} -a {input.assembly} --{wildcards.tool} {input.R1} {input.R2}
        touch {output.stats}
        """


rule metawrap_refine_binning:
    input:
        binputs=binstats,
    output:
        stats=f'metawrap/refined_binning_{{sample}}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
        stats_mqc=f'metawrap/refined_binning_{{sample}}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats_mqc.tsv',
        contigs=f'metawrap/refined_binning_{{sample}}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.contigs',
        bin1=f"metawrap/refined_binning_{{sample}}/{BINNING_TOOLS[0]}_bins.stats",
        bin2=f"metawrap/refined_binning_{{sample}}/{BINNING_TOOLS[1]}_bins.stats",
        bin3=f"metawrap/refined_binning_{{sample}}/{BINNING_TOOLS[2]}_bins.stats",
    params:
        outdir=lambda wc, output: os.path.dirname(output.stats),
        binput_dirs=lambda wc, input: [os.path.dirname(x) for x in input.binputs],
        completeness=config["metawrap_compl_thresh"],
        contamination=config["metawrap_contam_thresh"],
    # this is a different container that has checkm installed
    container:
        config["docker_metawrap_refine"]
    threads: 32
    resources:
        mem_mb=32 * 1024,
        runtime="12:00",
    shell:
        """
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.binput_dirs[0]} -B {params.binput_dirs[1]} -C {params.binput_dirs[2]} -c {params.completeness} -x {params.contamination}
        echo -e "#id: 'metawrap'\n#plot_type: 'table'\n#section_name: 'Bin Refinement'" > {output.stats}_mqc.tsv && cat {output.stats} >> {output.stats}_mqc.tsv
        """
