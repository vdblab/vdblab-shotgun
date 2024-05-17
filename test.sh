#! /bin/bash
set -ux
set -o pipefail
set -o errexit

stage=$1
rawdataset=$2
case $rawdataset in
    tiny)
	echo " WARNING: this dataset will raise errors during binning/annotation and will have empty metaphlan results due to its size"
	nshards=1
	R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz]
	R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz]
	addnconf="dedup_platform=HiSeq"
	;;
    tinybz)
	echo " WARNING: this dataset will raise errors during binning/annotation and will have empty metaphlan results due to its size"
	nshards=1
	R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.bz2]
	R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.bz2]
	addnconf="dedup_platform=HiSeq"
	;;
    equalcov)
	nshards=1
	R1=[$PWD/.test/simulated/1_depth100000_equalcoverage_R1.fastq.gz]
	R2=[$PWD/.test/simulated/1_depth100000_equalcoverage_R2.fastq.gz]
	addnconf="dedup_platform=SRA" # art doesnt give Illumina headers
	;;
    equalreads)
	nshards=1
	if [[ "$stage" == pre* ]] # note double brackets
	   then
	       # simulated a high duplicate library
	       if [ ! -f "$PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz" ];
		  then
		      mkdir -p $PWD/.test/simulated/duplicated/
		      cat $PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz > $PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz
		      cat $PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz >> $PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz
		      cat $PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz > $PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz
		      cat $PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz >> $PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz
	       fi
	       R1=[$PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz]
	       R2=[$PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz]
	else
	    R1=[$PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz]
	    R2=[$PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz]
	fi

	addnconf="dedup_platform=SRA" # art doesnt give Illumina headers
	;;
    equalreads4shards )
	nshards=4
	if [[ "$stage" == pre* ]] # note double brackets
	   then
	       if [ ! -f "$PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz" ];
	       then
		   # simulated a high duplicate library
		   mkdir -p $PWD/.test/simulated/duplicated/
		   cat $PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz > $PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz
		   cat $PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz >> $PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz
		   cat $PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz > $PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz
		   cat $PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz >> $PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz
	       fi
	       R1=[$PWD/.test/simulated/duplicated/1_depth100000_R1.fastq.gz]
	       R2=[$PWD/.test/simulated/duplicated/1_depth100000_R2.fastq.gz]
	else
	    R1=[$PWD/.test/simulated/1_depth100000_equalreads_R1.fastq.gz]
	    R2=[$PWD/.test/simulated/1_depth100000_equalreads_R2.fastq.gz]
	fi

	addnconf="dedup_platform=SRA" # art doesnt give Illumina headers
	;;
    small)
	nshards=2
	if [ ! -d "${PWD}/.test/SRR18369973/" ]
	then
	    echo "${PWD}/.test/SRR18369973/ not found; please run the getdata.sh script found in .test/ to fetch two test datasets"
	    exit 1
	fi
	R1=[${PWD}/.test/SRR18369973/SRR18369973_1.fastq.gz]
	R2=[${PWD}/.test/SRR18369973/SRR18369973_2.fastq.gz]
	addnconf="dedup_platform=SRA"
	;;
    multilib)
	# this is to test handling of multiple fastqs (lanes, typically)
	nshards=2
	R1=[${PWD}/.test/SRR21986403/SRR21986403_1.fastq.gz,${PWD}/.test/SRR18369973/SRR18369973_1.fastq.gz]
	R2=[${PWD}/.test/SRR21986403/SRR21986403_2.fastq.gz,${PWD}/.test/SRR18369973/SRR18369973_2.fastq.gz]
	addnconf="dedup_platform=SRA"
	;;
    medium)
	nshards=2
	R1=[${PWD}/.test/SRR21986403/SRR21986403_1.fastq.gz]
	R2=[${PWD}/.test/SRR21986403/SRR21986403_2.fastq.gz]
	addnconf="dedup_platform=SRA"
	;;
    *)
	echo -e "unknown dataset; please chose from tiny, small, medium, or multilib. Exiting\n"
	exit 1
	;;
esac
echo $R1

common_args="--snakefile workflow/Snakefile  --rerun-incomplete --restart-times 0 --cores 32"
case $stage in

     all)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/,/scratch/" \
	    --directory tmp${stage}_${rawdataset}/ \
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
	    --directory tmp${stage}_${rawdataset}/ \
	    --notemp \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    nshards=$nshards \
	    stage=preprocess
	;;
    preprocess-se )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --notemp \
	    --config \
	    sample=473  \
	    lib_layout=single \
	    R1=$R1 \
	    $addnconf \
	    nshards=$nshards \
	    stage=preprocess
	;;
    preprocess-gha )
	# runs just to dedup for easy execution on github actions
	snakemake \
	    $common_args \
            --cores 1 \
            --jobs 1 \
            --resources mem_mb=5000 \
	    --use-singularity \
            --singularity-prefix /github/workspace/.singularity/ \
            --singularity-args '-B /github/' \
	    --directory tmp${stage}_${rawdataset}/ \
	    --notemp \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    multiqc_config=${PWD}/multiqc_config.yaml \
	    nshards=$nshards \
	    stage=preprocess -f dedup/473_R1.fastq.gz
	;;
    testpreprocess )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmppreprocess_testing/   \
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
	    --directory tmppreprocess_testing/   \
	    --notemp \
	    --config \
	    sample=473  \
	    R1=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	    R2=[$PWD/.test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
	    $addnconf \
	  nshards=$nshards \
	  stage=preprocess
	;;
    biobakery )
	# use the preprocessed, host-depleted input so we can get better estimated reads from metaphlan
	if [ ! -d "tmppreprocess_${rawdataset}/" ]
	then
	    echo "tmppreprocess_${rawdataset}/ not present; please run `bash test.sh preprocess $rawdataset`"
	    exit 1
	fi
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config \
	    sample=473  \
	    R1=[$PWD/tmppreprocess_${rawdataset}/hostdepleted/473_R1.fastq.gz] \
	    R2=[$PWD/tmppreprocess_${rawdataset}/hostdepleted/473_R2.fastq.gz] \
	    $addnconf \
	    stage=biobakery -f metaphlan/473_metaphlan3_profile.txt
	;;

    biobakery-se )
	# use the preprocessed, host-depleted input so we can get better estimated reads from metaphlan
	if [ ! -d "tmppreprocess-se_${rawdataset}/" ]
	then
	    echo "tmppreprocess-se_${rawdataset}/ not present; please run `bash test.sh preprocess-se $rawdataset`"
	    exit 1
	fi

	# just run through metaphlan in the interest of time; humann needs lots of resources
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config \
	    sample=473  \
	    R1=[$PWD/tmppreprocess-se_${rawdataset}/hostdepleted/473_R1.fastq.gz] \
	    lib_layout=single \
	    $addnconf \
	    stage=biobakery -f metaphlan/473_metaphlan3_profile.txt
	;;

    mtx )
	snakemake \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --snakefile workflow/Snakefile_mtx \
	    --directory tmp${stage}_${rawdataset}/ \
            --config \
            sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    mpa_profile=/data/brinkvd/data/shotgun/test/C011815_metaphlan3_profile.txt
	;;

    kraken )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config \
	    sample=473  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    kraken2_db=/data/brinkvd/resources/dbs/kraken/k2_pluspf_08gb_20230314/ \
	    stage=kraken
	;;
    assembly )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=assembly assembler=spades
	;;
    bin)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config sample=473 \
	    assembly=${PWD}/.test/473/473.assembly.fasta  \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=binning
	;;
    annotate)
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    assembly=${PWD}/.test/473/473.assembly.fasta  \
	    stage=annotate
	;;
    rgi )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmprgi_${rawdataset}/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    $addnconf \
	    stage=rgi
	;;
    downsample|ds )
	snakemake \
	    $common_args \
	    --singularity-args "-B ${PWD},/data/brinkvd/" \
	    --directory tmp${stage}_${rawdataset}/ \
	    --config sample=473 \
	    R1=$R1 \
	    R2=$R2 \
	    depths=[1000] \
	    reps=[1,2] \
	    $addnconf \
	    stage=downsample
	;;

    figs )
	for stage in all preprocess biobakery binning kraken assembly annotate rgi
	do
	    snakemake \
		$common_args \
		--singularity-args "-B ${PWD},/data/brinkvd/" \
		--directory tmp${stage}_${rawdataset}/ \
		--config sample=473 \
		R1=$R1 \
		R2=$R2 \
		$addnconf \
		nshards=2 \
		assembly=${PWD}/.test/assembly.fna  \
		stage=$stage --dag  | sed "s|color.*rounded\"|color = \"grey\", style=\"rounded\"|g" > images/${stage}_dag.dot &&  dot -Tpng images/${stage}_dag.dot -o images/${stage}_dag.png
	done

	;;
    *)
	echo -e "unknown stage; please chose from all, preprocess, preprocess-se, preprocess-gha, biobakery, biobakery-se, bin, kraken, assembly, annotate, rgi, figs. Exiting\n"
	;;
esac
