import gzip
import os
import sys

from Bio import SeqIO
from functools import partial

clusters = []

infile = snakemake.input.gbk #infile = sys.argv[1]

open_fun = partial(gzip.open, mode='rt') if  infile.endswith(".gz") else open

# if infile.endswith(".gz"):
#     open_fun = gzip.opn
# else:
#     open_fun = open

with open_fun(infile) as inf:
    for seqRecord in SeqIO.parse(inf, "genbank"):
        for seqFeature in seqRecord.features:
            if seqFeature.type == "cluster" or True:
                clusters.append(
                    {'contigID': seqRecord.id,
                     'start': seqFeature.location.start,
                     'end': seqFeature.location.end,
                     'gctype': "Unknown"})
                if "product" in seqFeature.qualifiers:
                    clusters[-1]['gctype'] = seqFeature.qualifiers['product'][0]

with open (snakemake.output.tab, "w") as outf:
    for clus in clusters:
        outf.write("\t".join([str(x) for x in [
            os.path.basename(infile),
            clus["contigID"],
            clus["start"],
            clus["end"],
            clus["gctype"]]]) + "\n")
