#! /bin/bash
set -eux

mode=$1

case $mode in

  full | all)
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --snakefile workflow/Snakefile \
	  --directory tmpall/ \
	  --config \
	  sample=473 \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  multiqc_config=${PWD}/multiqc_config.yaml nshards=1 \
	  dedup_reads=False kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/ \
	  stage=all
      ;;
  preprocess )
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --snakefile workflow/Snakefile \
	  --directory tmppre/   \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  multiqc_config=${PWD}/multiqc_config.yaml \
	  nshards=2 \
	  dedup_reads=False kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/ \
	  stage=preprocess
      ;;
  biobakery | bb)
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --snakefile workflow/Snakefile \
	  --directory tmpbio/ \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz]  \
	  stage=biobakery
      ;;

  mtx )
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --snakefile workflow/Snakefile_mtx \
          --directory tmpmtx/ \
          --config \
          sample=473  \
          R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
          R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  mpa_profile=/data/brinkvd/data/shotgun/test/C011815_metaphlan3_profile.txt
      ;;

  kraken | kraken2 )
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --snakefile workflow/Snakefile \
	  --directory tmpkraken/ \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/ \
	  stage=kraken
      ;;
  assembly)
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile workflow/Snakefile \
	  --directory tmpassembly/ \
	  --config sample=473 \
	  R1=[${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz] \
	  R2=[${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz] \
	  stage=assembly
      ;;
  bin)
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile workflow/Snakefile \
	  --directory tmpbinning/ \
	  --config sample=473 \
	  assembly=${PWD}/tmpassembly/473.assembly.fasta \
	  R1=[${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz] \
	  R2=[${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz] \
	  stage=binning
      ;;
  annotate)
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile workflow/Snakefile \
	  --directory tmpannotate/ \
	  --config sample=473 \
	  assembly=${PWD}/tmpassembly/473.assembly.fasta  \
	  stage=annotate
      ;;

  rgi )
      snakemake \
	  --profile ${SNAKEPROFILE} \
	  --singularity-args "--bind $PWD" \
	  --snakefile workflow/Snakefile \
	  --directory tmprgi/ \
	  --config sample=473 \
	  R1=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz \
	  R2=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz \
	  stage=rgi
      ;;
  *)
    echo -e "unknown mode; please chose from all, preprocess, biobakery, bin, kraken2, assembly, annotate, rgi. Exiting\n"
    ;;
esac
