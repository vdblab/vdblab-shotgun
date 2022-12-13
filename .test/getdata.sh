

for srr in SRR18369973	SRR21986403 SRR19479275
do

if [ ! -f "$PWD/$srr/${srr}_1.fastq.gz" ]
then

    singularity run -B "$PWD" docker://hnh0303/sratools:2.10.0 fasterq-dump $srr -O $PWD/$srr/
    gzip $PWD/$srr/*.fastq
fi
done
