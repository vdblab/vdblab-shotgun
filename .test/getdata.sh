set -eux
# Get raw SRR:
for srr in SRR18369973	SRR21986403 SRR19479275
do

if [ ! -f "$PWD/$srr/${srr}_1.fastq.gz" ]
then

    singularity run -B "$PWD" docker://hnh0303/sratools:2.10.0 fasterq-dump $srr -O $PWD/$srr/
    gzip $PWD/$srr/*.fastq
fi
done

# Get assembly for annotation tests:
# (We will fetch E coli K-12 MG16555 and B sub 168) - then create a modified version of the E coli genome with a mannose carbohydrate utilization gene added.
# Picking mannose as a carbohydrate utilizaiton gene to include based on this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3558604/#:~:text=sugars%20(galactose%2C%20hexuronates%2C-,mannose,-and%20ribose)%20that
# GCF_000005845.2 = E coli MG1655
# GCF_000009045.1 = B sub 168 
for gcf in GCF_000005845.2 GCF_000009045.1; do
if [ ! -d "$PWD/$gcf/" ]; then
    singularity run -B "$PWD" docker://biocontainers/ncbi-datasets-cli:15.12.0_cv23.1.0-4 datasets download genome accession $gcf --include genome,gff3
    unzip ncbi_dataset.zip
    mv ncbi_dataset/data/$gcf/ $PWD/$gcf/
    for tmpfile in ncbi_dataset.zip ncbi_dataset README.md; do
        rm -rf $tmpfile
    done
fi
done

# We will pull out the B sub manA gene:
