
singularity run -B "$PWD" docker://hnh0303/sratools:2.10.0 fasterq-dump SRR21986403 -O $PWD/SRR21986403/
gzip $PWD/SRR21986403/*.fastq
