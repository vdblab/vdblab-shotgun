import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


# the str is needed as when running from github workflow.current_basedir is a Githubfile, not a string or a path so os.path objects
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])


onstart:
    with open("config_used.yaml", "w") as outfile:
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
covermreports = f'coverm/{config["sample"]}_metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.coverage_mqc.tsv'


rule all:
    input:
        binstats,
        refined_stats,
        refined_binstats_each,
        stats_mqc,
        contigs,
        covermreports,


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
        outdir=lambda wc, output: os.path.dirname(os.path.dirname(output.stats)),
    container:
        config["docker_metawrap"]
    threads: 64
    resources:
        mem_mb=32 * 1024,
        runtime=12 * 60,
    shell:
        """
        metawrap binning -o {params.outdir} -t {threads} -a {input.assembly} --{wildcards.tool} {input.R1} {input.R2}
        touch {output.stats}
        """


rule metawrap_refine_binning:
    """The names for the params binput_dirs is due to the crazy naming of the output of metawrap.
    The only consistant file we can use as a trigger is the <tool>.done file, but the refine module needs the dir beneath it under a path like
    metawrap/rawbinning_473/concoct/concoct_bins/concoct_bins
    """
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
        checkm_db=config["checkm_db"],
    # this is a different container that has checkm installed
    container:
        config["docker_metawrap"]
    threads: 32
    # give it 82 gb memory because checkm estimates pplacer will need 40GB per core, and we want this to run reasonably fast
    resources:
        mem_mb=82 * 1024,
        runtime=12 * 60,
    shell:
        """
        export  CHECKM_DATA_PATH={params.checkm_db}
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.binput_dirs[0]} -B {params.binput_dirs[1]} -C {params.binput_dirs[2]} -c {params.completeness} -x {params.contamination}
        echo -e "#id: 'metawrap'\n#plot_type: 'table'\n#section_name: 'Bin Refinement'" > {output.stats}_mqc.tsv && cat {output.stats} >> {output.stats}_mqc.tsv
        """


rule coverm:
    """ This calculates bin coverage.  The bin stats file is used as the
    trigger, but we actually want the directory of the bins

    specifying inputs with --coupled and -1 -2 seem to give equivalent results.
    Note that the output only refers tothe forward file as the sample name.
    minimap2 is supposedly faster and more accurate than BWA, so we use that
    as the mapper.
    we return all the available methods as of the time of writing this rule
    except for coverage_histogram which has to run separately
    ("Cannot specify the coverage_histogram method with any other coverage methods")

    --min-covered-fraction is required to be set to 0 for certain cov metrics
     --genome-fasta-extension is fa since thats how metawrap outputs it
    """
    input:
        R1=config["R1"],
        R2=config["R2"],
        stats=f'metawrap/refined_binning_{{sample}}/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
    output:
        mqc=f'coverm/{{sample}}_metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.coverage_mqc.tsv',
        bams=directory(f'coverm/{{sample}}_metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bams/'),
    params:
        bindir=lambda wc, input: input.stats.replace(".stats", ""),
        fastq_string = lambda wc, input: " ".join([config["R1"][x] + " " + config["R2"][x] for x in range(0, len(config["R1"]))])
    container:
        config["docker_coverm"]
    threads: 16
    shell:
        """
        coverm genome --genome-fasta-directory {params.bindir} \
          --coupled {params.fastq_string} \
          --mapper minimap2-sr \
          --methods mean relative_abundance trimmed_mean \
            covered_bases variance length count reads_per_base rpkm tpm \
          --output-file {output.mqc}.tmp --threads {threads} \
          --bam-file-cache-directory {output.bams} \
          --min-covered-fraction 0 \
          --genome-fasta-extension fa
        # add in the Multiqc header info
        echo -e "# plot_type: 'table'\n# section_name: 'Bin Coverage Statistics'" > {output.mqc}
        cat {output.mqc}.tmp >> {output.mqc}
        rm {output.mqc}.tmp
        """
