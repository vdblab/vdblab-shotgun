#! /usr/bin/env python3
import os
import pandas as pd
from statistics import  mean, stdev, StatisticsError

def process_sortmerna(paths, outpath,  sample_name):
    res = {}
    for f in paths:
        with open(f, "r")  as inf:
            path = os.path.basename(f)
            for line in inf:
                line = line.strip()
                if "Total reads =" in line:
                    res[path] = {"total": line.split("=")[1].split()[0],
                              "hits": None,
                              "hit_perc": None}
                elif "Total reads passing E-value threshold" in line:
                    res[path]["hits"] = line.split("=")[1].split()[0]
                    res[path]["hit_perc"] = line.split("=")[1].split()[1].replace("(", "").replace(")", "")
                elif "Total reads failing E-value threshold" in line:
                    res[path]["nohits"] = line.split("=")[1].split()[0]
                    res[path]["nohit_perc"] = line.split("=")[1].split()[1].replace("(", "").replace(")", "")
                else:
                    pass

    total_reads = sum([int(v["total"]) for k, v in res.items()])
    total_hits = sum([int(v["hits"]) for k, v in res.items()])
    total_nohits = sum([int(v["nohits"]) for k, v in res.items()])
    mean_perc = round(mean([float(v["hit_perc"]) for k, v in res.items()]), 4)
    mean_noperc = round(mean([float(v["nohit_perc"]) for k, v in res.items()]), 4)
    try:
        stdev_perc = round(stdev([float(v["hit_perc"]) for k, v in res.items()]), 4)
    except StatisticsError as e:
        stdev_perc = 0
    multiqc_header = "# id: 'sortmerna'\n# plot_type: 'table'\n#section_name: 'SortMeRNA'"
    with open(paths[0], "r")  as inf, open(outpath, "w") as outf:
        for line in inf:
            if "Total reads =" in line:
                outf.write(line.split("=")[0] + f"= {total_reads}\n")
            elif "Total reads passing E-value threshold" in line:
                outf.write(line.split("=")[0] + f"= {total_hits} ({mean_perc})\n")
            elif "Total reads failing E-value threshold" in line:
                outf.write(line.split("=")[0] + f"= {total_nohits} ({mean_noperc})\n")
            elif "Reads file: " in line:
                outf.write(f"    Reads file: {sample_name}\n")
            else:
                outf.write(line)

    return([
        multiqc_header,
        "\t".join([str(x) for x in ["sample", "total_reads", "total_hits", "mean_perc", "stdev_perc"]]),
        "\t".join([str(x) for x in [sample_name, total_reads, total_hits, mean_perc, stdev_perc]])
        ])

def process_kneaddata(path, sample_name):
    """ the resutls from teh 4 stage host depletion scheme are processed the same way
    """
    vals = []
    with open(path, "r") as inf:
        for i, line in enumerate(inf):
            if i == 0:
                header = line.strip()
            else:
                sline = line.strip().split("\t")[1:]
                vals.append(sline)
    sums = [sample_name]
    # for each metric, sum the reads across all shards
    for metrici, v  in enumerate(vals[0]):
        sums.append(sum([float(x[metrici]) for x in vals]))
    multiqc_header = "# id: 'kneaddata'\n# plot_type: 'table'\n#section_name: 'Host Depletion'"
    return([
        multiqc_header,
        header,
        "\t".join([str(x) for x in sums])
        ])

def process_bbdup_trimming(paths, sample_name, outpath):
    res = {}
    resdf = []
    for f in paths:
        # parse the header with string manipulation
        with open(f, "r")  as inf:
            path = os.path.basename(f)
            for line in inf:
                line = line.strip()
                if "#Total" in line:
                    res[path] = {"total": int(line.split("\t")[1]),
                              "hits": None,
                              "hit_perc": None}
                elif "#Matched" in line:
                    res[path]["hits"] = int(line.split("\t")[1].split()[0])
                    #res[path]["hit_perc"] = round(res[path]["hits"]/res[path]["total"])
                else:
                    pass
        # parse the bulk of the stuff
        tmp = pd.read_csv(f, sep="\t", comment="#", names=["Name","Reads", "ReadsPct"])
        tmp["file"] = f
        resdf.append(tmp[["file", "Name", "Reads"]])

    total_reads = sum([int(v["total"]) for k, v in res.items()])
    total_hits = sum([int(v["hits"]) for k, v in res.items()])
    hit_perc = str(round(total_hits/total_reads * 100, 8)  ) + "%"
    # combine all them, group_by the hit name, get sum
    bodydf = pd.concat(resdf).groupby(["Name"])["Reads"].sum().reset_index()
    bodydf["Reads"] = bodydf["Reads"].astype(int)
    # recalculate the percentages
    pd.options.display.float_format = '{:.8f}'.format
    bodydf["ReadsPCT"] = bodydf["Reads"] / total_reads * 100
    bodydf["ReadsPCT"] = bodydf.ReadsPCT.map(lambda x: '{:.8f}'.format(x)) + "%"

    # write out the header, then add the merged table
    with open(outpath, "w") as outf:
        outf.write(f"#File\t{sample_name}\n")
        outf.write(f"#Total\t{total_reads}\n")
        outf.write(f"#Matched\t{total_hits}\t{hit_perc}\n")
        outf.write(f"#Name\tReads\tReadsPct\n")
    bodydf.to_csv(outpath, mode='a', index=False,header=False, sep="\t")




def write_table(lines, path):
    with open(path, "w") as outf:
        for lin in lines:
            outf.write(lin + "\n")

if __name__ == "__main__":
    if "snakemake" not in globals():
        # replace this with some stable test data
        inlist =["/data/brinkvd/watersn/tmpout/sortmerna/Sample_FMT_147P_IGO_10170_1/Sample_FMT_147P_IGO_10170_1_sortmerna.log",
                 "/data/brinkvd/watersn/tmpout/sortmerna/Sample_739J_IGO_10170_4/Sample_739J_IGO_10170_4_sortmerna.log"]
        print(process_sortmerna(paths=inlist, outpath="sortmerna.log",  sample_name="test"))
        print(process_kneaddata(path="tmp.log", sample_name="tmp"))

        bbduk =["tmpout/trimmed/473_shard001_trimmingAQ.txt", "tmpout/trimmed/473_shard003_trimmingAQ.txt"]
        print(process_bbdup_trimming(paths=bbduk, outpath="trim.log",  sample_name="test"))
    else:
        with open(snakemake.log.e, "w") as ef, open(snakemake.log.o, "w") as of:
            sys.stderr,  sys.stdout = (ef, of)
            outlines = process_sortmerna(
                paths=snakemake.input.sortmernas,
                outpath=snakemake.output.sortmerna_log,
                sample_name=snakemake.params.sample_name)
            write_table(outlines, snakemake.output.sortmerna_report)

            outlines = process_kneaddata(path=snakemake.input.knead, sample_name=snakemake.params.sample_name)
            write_table(outlines, snakemake.output.knead)

            process_bbdup_trimming(
                paths=snakemake.input.bbtrim,
                outpath=snakemake.output.bbtrim,
                sample_name=snakemake.params.sample_name)
