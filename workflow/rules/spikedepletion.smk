import os
import json
import yaml
import shutil

from pathlib import Path
from snakemake.utils import validate


include: "common.smk"


# the str is needed as when running from github workflow.current_basedir is a Githubfile, not a string or a path so os.path objects
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")


# validate(
#     config,
#     os.path.join(str(workflow.current_basedir), "../../config/config.schema.yaml"),
# )


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
    ani_filter


covermreports = f'coverm/{config["sample"]}_bins.coverage_mqc.tsv'





rule all:
    input:
        covermreports,
        f"coverm/{config['sample']}_bins.coverage_mqc.tsv",
        f'{config["sample"]}_nospike_R1.fastq.gz',
        f'{config["sample"]}_nospike_R2.fastq.gz',

rule add_spikes_to_bin_dir:
    input:
        bindir=config["bindir"],
    output:
        bindir=temp(directory("{sample}_bins_with_spikes_raw"))
    shell:"""
    mkdir {output.bindir}
    cp {input.bindir}/* {output.bindir}/
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/013/045/GCA_000013045.1_ASM1304v1/GCA_000013045.1_ASM1304v1_genomic.fna.gz | gunzip > {output.bindir}/Salinibacter_ruber.fa
    # getting Trichoderma reesei  QM6a instead of  Trichoderma reesei ATCC 13631 .  Might be fine license-wise but not sure
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/167/675/GCF_000167675.1_v2.0/GCF_000167675.1_v2.0_genomic.fna.gz | gunzip  > {output.bindir}/Trichoderma_reesei.fa
    # getting CBA1121 instead of Haloarcula hispanica ATCC 33960  for the same reason.
    curl -o - ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/729/095/GCF_008729095.1_ASM872909v1/GCF_008729095.1_ASM872909v1_genomic.fna.gz | gunzip  > {output.bindir}/Haloarcula_hispanica.fa
    """

rule fastani:
    input:
        bindir="{sample}_bins_with_spikes_raw"
    output:
        csv="{sample}_ani_results.tsv",
        manifest="{sample}_manifest.txt"
    container:"docker://ghcr.io/vdblab/fastani:1.34"
    threads: 32
    shell: """
    ls {input.bindir}/*fa > {output.manifest}
    fastANI --ql {output.manifest} --rl {output.manifest} --minFraction 0 --threads {threads} -o {output.csv}
    """

rule ani_filter:
    """ 95 based on https://www.nature.com/articles/s41467-018-07641-9
    """
    input:
        csv="{sample}_ani_results.tsv",
        bindir="{sample}_bins_with_spikes_raw",
        manifest="{sample}_manifest.txt",
    output:
        pairs_for_vis="{sample}_close_ani.tsv",
        bindir=temp(directory("{sample}_bins_with_spikes")),
    params:
        ani_thresh = 95
    run:
        import os
        import sys
        import shutil

        contigs_to_ignore = []
        # could use pandas, but this is pretty straightforward
        with open(input.csv, "r") as inf, open(output.pairs_for_vis,  "w") as outf:
            for line in inf:
                (query, reference, ANI, count_mappings, total_frags) = line.strip().split()
                if query == reference:
                    continue
                # we generate vis for all comparisons exceeding 80%, just in case there are
                # any bin-bin comparisons of interest
                if float(ANI) > 80:
                    outf.write(line)
                # ... but only exclude contigs if they exceed the specified
                # threshold species-level
                if float(ANI) < params.ani_thresh:
                    continue
                if "Haloarcula_hispanica" in line or "Salinibacter_ruber" in line or "Trichoderma_reesei" in line:
                    # determine if the bin resembing the cannonical reference is the ani query or the reference
                    if "Haloarcula_hispanica" in query or "Salinibacter_ruber" in query or "Trichoderma_reesei" in query:
                        bin_index  = 1
                    else:
                        bin_index  = 0
                    contigs_to_ignore.append([query, reference][bin_index])
        os.makedirs(output.bindir)
        with open(input.manifest, "r") as inf:
            for line in inf:
                if line.strip() not in contigs_to_ignore:
                    dest =  os.path.join(output.bindir, os.path.basename(line.strip()))
                    if not os.path.exists(dest):
                        shutil.copyfile(line.strip(), dest)
                else:
                    print(f"Excluding {line.strip()}")
        print("done")



rule ani_vis:
    """ visualize the ones we exclude for being too close to the spikes,
    and any other with close ANI (>80%)
    """
    input:
        pairs_for_vis="{sample}_close_ani.tsv",
        bindir="{sample}_bins_with_spikes_raw",
    output:
        figdir=directory("{sample}_bin_ref_homology_figures"),
    container:"docker://ghcr.io/vdblab/fastani:1.34"
    params:
        ani_thresh=.9
    threads: 16
    shell: """
    mkdir {output.figdir}
    cat {input.pairs_for_vis} | cut -f 1,2,3 | while read query ref ANI
    do
        echo "Visualizing  $query vs  $ref (ANI of $ANI)"
        outbase={output.figdir}/$(basename $query)_$(basename $ref)_${{ANI}}
        fastANI -q $query -r $ref --minFraction 0 --threads {threads} --visualize -o ${{outbase}}.out

        Rscript /FastANI/scripts/visualize.R $query $ref  ${{outbase}}.out.visual
    done

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
        bindir=rules.ani_filter.output.bindir,
    output:
        mqc='coverm/{sample}_bins.coverage_mqc.tsv',
        bam=os.path.join('coverm/{sample}_bams/', 'coverm-genome.' + os.path.basename(config["R1"][0]) + ".bam"),
    params:
        fastq_string = lambda wc, input: " ".join([config["R1"][x] + " " + config["R2"][x] for x in range(0, len(config["R1"]))])
    container:
        config["docker_coverm"]
    threads: 32
    shell:
        """
        coverm genome --genome-fasta-directory {input.bindir} \
          --coupled {params.fastq_string} \
          --mapper minimap2-sr \
          --methods mean relative_abundance trimmed_mean \
            covered_bases variance length count reads_per_base rpkm tpm \
          --output-file {output.mqc}.tmp --threads {threads} \
          --bam-file-cache-directory $(dirname  {output.bam}) \
          --min-covered-fraction 0 \
          --genome-fasta-extension fa
        ls $(dirname  {output.bam})
        # add in the Multiqc header info
        echo -e "# plot_type: 'table'\n# section_name: 'Bin Coverage Statistics'" > {output.mqc}
        cat {output.mqc}.tmp >> {output.mqc}
        rm {output.mqc}.tmp
        """


rule get_nonspike_reads:
    input:
        bam=rules.coverm.output.bam,
    output:
        bed="{sample}_spike_regions.bed",
        nonspike=temp("{sample}_nonspike.bam"),
        R1="{sample}_nospike_R1.fastq.gz",
        R2="{sample}_nospike_R2.fastq.gz",
    threads: 4
    params:
        R1tmp=lambda wc, output: output.R1.replace(".gz", ""),
        R2tmp=lambda wc, output: output.R2.replace(".gz", ""),
    container: "docker://ghcr.io/vdblab/bowtie2:2.5.0"
    shell:"""
    # get header, turn into a BED file
    samtools view -H  {input.bam} | grep "Salinibacter\|Trichoderma\|Haloarcula" | cut -f 2,3 | sed "s|SN:||g" | sed "s|LN:|1\t|g" | sort -k1,1 -k2,2n > {output.bed}
    # inspired by https://www.biostars.org/p/473204/
    samtools view -L {output.bed} -U unsorted_{output.nonspike} -o /dev/null -@ {threads} {input.bam}
    # make sure you NAME-sort before running samtools fastq https://www.biostars.org/p/454942/
    samtools sort -n -o {output.nonspike} unsorted_{output.nonspike}
    rm unsorted_{output.nonspike}
    samtools fastq -1 {params.R1tmp} -2 {params.R2tmp} -0 /dev/null -s /dev/null -n -@ {threads} {output.nonspike}
    samtools flagstat {output.nonspike}
    pigz ./*fastq
    """
