import sys
import subprocess
import pandas as pd

main_dir = "/data/peledj/baichom1/Projects/Peled_Analysis/2024/01_Enterococcus_Nutrition_MTX/cazyme_pipeline/"
num_reads_file = "".join([main_dir, "num_reads_per_exp.txt"])
dbcan_dir = "".join([main_dir, "dbcan/"])
align_dir = "".join([main_dir, "bwa_align/"])
save_dir = "".join([main_dir, "detailed_dbcan/"])


def write_dbcan_info_file(
    aligned_file, substrate_file, cgc_file, overview_file, save_file, num_reads
):
    """
    Will handle merging the information from the various files into one summary file.
    """
    # contains an overview from run_dbcan
    overview_df = pd.read_csv(overview_file, sep="\t")

    # counts df currently contains the per-position alignment read count (of the 5' end of reads)
    # we need to summarize over the whole gene region to get our counts per gene.
    counts_df = pd.read_csv(
        aligned_file,
        sep="\t",
        names=["bam_name", "position", "counts", "length", "fraction"],
    ).drop(["position", "fraction"], axis="columns")
    counts_df = (
        counts_df.groupby("bam_name")
        .agg({"counts": "sum", "length": "max"})
        .reset_index()
    )

    # substrate file has the annotated substrates for various CGCs - "k_name" is the key:
    sub_df = pd.read_csv(
        substrate_file,
        header=None,
        skiprows=1,
        names=[
            "k_name_parts",
            "PULID",
            "substrate",
            "substrate_bitscore",
            "signature pairs",
            "dbCAN-sub substrate",
            "dbCAN-sub substrate score",
        ],
        sep="\t",
    )
    sub_df["k_name"] = sub_df["k_name_parts"].map(lambda x: x.split("|")[0])
    sub_df["cgc"] = sub_df["k_name_parts"].map(lambda x: x.split("|")[1])
    sub_df = sub_df[["k_name", "cgc", "substrate", "substrate_bitscore"]]

    # the cgc file acts as a "key file" connecting the "k_names"s used in the substrate file, and the bam_name names;
    cgc_faa_df = pd.read_csv(
        cgc_file,
        sep="\t",
        comment="+",
        names=[
            "_0",
            "caz_type",
            "_2",
            "_3",
            "cgc",
            "k_name",
            "_6",
            "_7",
            "bam_name",
            "_9",
            "_10",
            "_11",
        ],
    )
    cgc_faa_df = cgc_faa_df[["caz_type", "cgc", "k_name", "bam_name"]]

    caz_df = overview_df.merge(
        counts_df, how="left", left_on="Gene ID", right_on="bam_name"
    ).drop("bam_name", axis="columns")
    caz_df = caz_df.merge(
        cgc_faa_df, how="left", left_on="Gene ID", right_on="bam_name"
    )

    # Note that substrates can match with multiple "bam_names".
    # in the cgc file there is a one to many mapping of k_names/cgc's and bam_names.
    caz_df = caz_df.merge(sub_df, how="left", on=["k_name", "cgc"])
    caz_df["RPM"] = ((10**6) * caz_df["counts"]) / num_reads
    caz_df.to_csv(save_file, sep="\t")


def count_num_reads_compressed_file(file_name):
    """
    Small helper function which will take a compressed fastq file and return the number of reads.
        -file_name = str (or path) the location of the file to count number of reads for.
    """
    command = f"echo $(($(zcat {file_name} | wc -l)/4))"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    return float(result.stdout)


if __name__ == "__main__":
    """
    Arguments required (given in the order expected:):
        - aligned_file = the depth file generated from aligning the host-depleted
        reads to the metaerg annotated genes.
        - substrate_file = the run_dbcan substrate file.
        - cgc_file = the run_dbcan cgc.out file.
        - overview_file = the run_dbcan overview file.
        - save_file = name of the file to write the RPM outputs to.
        - r1 = the location of the r1 files for calculating the number of reads in the host-depleted file.
    """
    if "snakemake" not in globals():
        # assume taking in variables from the main argv:
        aligned_file = sys.argv[1]
        substrate_file = sys.argv[2]
        cgc_file = sys.argv[3]
        overview_file = sys.argv[4]
        save_file = sys.argv[5]
        r1 = sys.argv[6]
    else:
        aligned_file = snakemake.input.coverage
        substrate_file = snakemake.input.substrate
        cgc_file = snakemake.input.cgc
        overview_file = snakemake.input.overview
        save_file = snakemake.output.rpm_file
        r1 = snakemake.input.r1
    num_reads = count_num_reads_compressed_file(str(r1))
    print(f"Num reads: {num_reads}")
    write_dbcan_info_file(
        aligned_file, substrate_file, cgc_file, overview_file, save_file, num_reads
    )
