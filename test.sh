#! /bin/bash
set -ux
set -o pipefail
set -o errexit

mode=$1
rawdataset=$2
case $rawdataset in
    tiny)
	echo " WARNING: this dataset will raise errors during binning/annotation and will have empty metaphlan results due to its size"
	nshards=1
	R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz]
	R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz]
	addnconf="dedup_platform=HiSeq"
	;;
    small)
	nshards=2
	R1=[${PWD}/.test/SRR18369973/SRR18369973_1.fastq.gz]
	R2=[${PWD}/.test/SRR18369973/SRR18369973_2.fastq.gz]
	addnconf="dedup_platform=SRA"
	;;
    medium)
	nshards=2
	R1=[${PWD}/.test/SRR21986403/SRR21986403_1.fastq.gz]
	R2=[${PWD}/.test/SRR21986403/SRR21986403_2.fastq.gz]
	addnconf="dedup_platform=SRA"
	;;
    *)
	echo -e "unknown dataset; please chose from tiny. Exiting\n"
	exit 1
	;;
esac



common_args="--snakefile workflow/Snakefile --profile ${SNAKEPROFILE} --rerun-incomplete --restart-times 0 --cores 32"
case $mode in

    full | all)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/,/scratch/" \
            --directory tmpall_${rawdataset}/ \
	    --config \
	    sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    multiqc_config=${PWD}/multiqc_config.yaml nshards=$nshards \
	    dedup_reads=False \
	    stage=all
	;;
    preprocess )
	# the --notemp is here so we can do the unittests afterward
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
            --directory tmppre_${rawdataset}/ \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    multiqc_config=${PWD}/multiqc_config.yaml \
	    nshards=$nshards \
	    stage=preprocess
	;;
    testpreprocess )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmppre_testing/   \
	    --config \
	    sample=473  \
	    R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	    R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	    $addnconf \
	    nshards=$nshards \
	    stage=preprocess \
	    --generate-unit-tests
	;;
    testpreprocess_build )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmppre_testing/   \
	    --notemp \
	    --config \
	    sample=473  \
	    R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	    R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	    $addnconf \
	  nshards=$nshards \
	  stage=preprocess
	;;
    biobakery | bb)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmpbio/ \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=biobakery
	;;

    mtx )
	snakemake \
	    --profile ${SNAKEPROFILE} \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --snakefile workflow/Snakefile_mtx \
            --directory tmpmtx/ \
            --config \
            sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    mpa_profile=/data/brinkvd/data/shotgun/test/C011815_metaphlan3_profile.txt
	;;

    kraken | kraken2 )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmpkraken/ \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/ \
	    stage=kraken
	;;
    assembly)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmpassembly/ \
	    --config sample=473 \
	    R1=[${PWD}/.test/SRR21986403/SRR21986403_1.fastq.gz] \
	    R2=[${PWD}/.test/SRR21986403/SRR21986403_2.fastq.gz] \
	    stage=assembly
	;;
    bin)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --config sample=473 \
	    --directory tmpbin_$dataset/ \
	    assembly=${PWD}/tmpassembly/473.assembly.fasta  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=binning
	;;
    annotate)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmpannotate/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    assembly=${PWD}/tmpassembly/473.assembly.fasta  \
	    stage=annotate
	;;
    rgi )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmprgi/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=rgi
	;;

  figs )
        for stage in all preprocess biobakery binning kraken assembly annotate rgi
	do
	    snakemake \
		$common_args \
		--singularity-args "-B ${PWD},/data/brinkvd/" \
		--directory tmprgi/ \
		--config sample=473 \
		R1=$R1 \
		R2=$R2 \
		$addnconf \
		nshards=2 \
		assembly=${PWD}/tmpassembly/473.assembly.fasta  \
		stage=$stage --dag  | sed "s|color.*rounded\"|color = \"grey\", style=\"rounded\"|g" > images/${stage}_dag.dot &&  dot -Tpng images/${stage}_dag.dot -o images/${stage}_dag.png
	done

	;;
    *)
	echo -e "unknown mode; please chose from all, preprocess, biobakery, bin, kraken2, assembly, annotate, rgi, figs. Exiting\n"
	;;
esac
