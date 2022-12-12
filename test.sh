#! /bin/bash
set -eux

mode=$1

case $mode in

  full | all)
      snakemake \
	  --profile ${PWD}/msk-lsf/ \
	  --snakefile vdb_shotgun/Snakefile \
	  --directory tmpall/ \
	  --config \
	  sample=473 \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  multiqc_config=${PWD}/vdb_shotgun/multiqc_config.yaml nshards=1 \
	  dedup_reads=False kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/
      ;;
  preprocess )
      snakemake \
	  --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/rules/preprocess.smk \
	  --directory tmppre/   \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  multiqc_config=${PWD}/vdb_shotgun/multiqc_config.yaml \
	  nshards=2 \
	  dedup_reads=False kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/
      ;;
  biobakery | bb)
      snakemake \
	  --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/rules/biobakery.smk \
	  --directory tmpbio/ \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
      ;;

  mtx )
      snakemake \
          --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/Snakefile_mtx \
          --directory tmpmtx/ \
          --config \
          sample=473  \
          R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
          R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  mpa_profile=/data/brinkvd/data/shotgun/test/C011815_metaphlan3_profile.txt
      ;;

  kraken | kraken2 )
      snakemake \
	  --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/rules/kraken.smk \
	  --directory tmpkraken/ \
	  --config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	  kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/
      ;;
  assembly)
      snakemake \
	  --profile ${PWD}/msk-lsf/ --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile vdb_shotgun/rules/assembly.smk  --directory tmpassembly/ \
	  --config sample=473 \
	  R1=[${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz] \
	  R2=[${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz]
      ;;
  bin)
      snakemake \
	  --profile ${PWD}/msk-lsf/ --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile vdb_shotgun/rules/binning.smk  --directory tmpbinning/ \
	  --config sample=473 \
	  assembly=${PWD}/tmpassembly/473.assembly.fasta \
	  R1=[${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz] \
	  R2=[${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz]
      ;;
  annotate)
      snakemake \
	  --profile ${PWD}/msk-lsf/ --cores 32 --singularity-args "-B ${PWD}" \
	  --snakefile vdb_shotgun/rules/annotate.smk  --directory tmpannotate/ \
	  --config sample=473 \
	  assembly=${PWD}/tmpassembly/473.assembly.fasta \
      ;;

  rgi )
      snakemake \
      --profile ${PWD}/msk-lsf/ --singularity-args "--bind $PWD" \
      --snakefile vdb_shotgun/rules/RGI.smk  --directory tmprgi/ \
      --config sample=473 \
      R1=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz \
      R2=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz
      ;;
  *)
    echo -e "unknown mode; please chose from all, preprocess, biobakery, bin, kraken2, assembly, annotate, rgi. Exiting\n"
    ;;
esac

