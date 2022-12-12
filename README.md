# Shotgun Pipeline

## Main Pipeline
### Usage

Single Sample:
```sh
snakemake \
  --profile ${PWD}/msk-lsf/ \
  --snakefile vdb_shotgun/Snakefile \
  --directory tmpout/ \
  --config \
    sample=473 \
	R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
    nshards=4 \
  --dry-run
```


### Outputs

* MultiQC-ready reports
* Microbe relative abundances (MetaPhlAn3, Kraken2)
* Metabolic pathway relative abundances (HUMAnN3)
* Metagenome assembled genomes (MetaSPAdes)
* AMR profiles with Abricate and RGI
* MAGs with MetaWRAP (Metabat2, CONCOCT, Metabin2)
*

### Workflow

The rule DAG for a single sample looks like this:

<img src="https://github.com/vdblab/vdblab-pipelines/raw/main/vdb_shotgun/images/main_dag.png" alt="Main Shotgun Pipeline DAG" width="900">

The following tools are used:

* BBTools ([site](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | [paper](https://doi.org/10.1371/journal.pone.0185056))
* SeqKit ([site](https://bioinf.shenwei.me/seqkit/) | [paper](https://doi.org/10.1371/journal.pone.0163962))
* KneadData ([site](https://huttenhower.sph.harvard.edu/kneaddata/))
* SortMeRNA ([site](https://bioinfo.lifl.fr/RNA/sortmerna/) | [paper](https://doi.org/10.1093/bioinformatics/bts611))
* FastQC ([site](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
* MultiQC ([site](https://multiqc.info) | [paper](http://dx.doi.org/10.1093/bioinformatics/btw354))
* Kraken2 ([site](https://ccb.jhu.edu/software/kraken2/) | [paper](https://doi.org/10.1186/s13059-019-1891-0))
* MetaSPAdes ([site](http://cab.spbu.ru/software/spades/) | [paper](https://doi.org/10.1101/gr.213959.116))
* MetaQUAST ([site](http://quast.sourceforge.net/metaquast) | [paper](https://doi.org/10.1093/bioinformatics/btv697))
* MetaPhlAn3 ([site](http://huttenhower.sph.harvard.edu/metaphlan) | [paper](https://doi.org/10.7554/eLife.65088))
* HUMAnN3 ([site](https://huttenhower.sph.harvard.edu/humann) | [paper](https://doi.org/10.7554/eLife.65088))
* RGI ([site](https://github.com/arpcard/rgi) | [paper](https://doi.org/10.1093/nar/gkz935))


## Annotation Pipeline
This pipeline will annotate metagenome assembled contigs.

### Usage

```sh
snakemake \
  --profile msk-lsf/ \
  --snakefile vdb_shotgun/annotate.smk \
  --directory tmpassembly/ \
  --config assembly=/full/path/to/contigs.fasta \
  --dry-run
```

### Outputs
* Gene prediction and annotation (MetaErg)
* Secondary metabolite gene clusters (antiSMASH)
* Antimicrobial resistance and virulence genes (ABRicate, AMRFinderPlus)
* Carbohydrate active enzyme (CAZyme) annotation (dbCAN3)

### Workflow

This is the rule DAG:

<img src="https://github.com/vdblab/vdblab-pipelines/raw/main/vdb_shotgun/images/annotate_dag.png" alt="Annotate Shotgun Pipeline DAG" width="500">

The following tools are used:

* MetaErg ([site](https://github.com/xiaoli-dong/metaerg) | [paper](https://doi.org/10.3389/fgene.2019.00999))
* antiSMASH ([site](https://antismash.secondarymetabolites.org/#!/about) | [paper](https://doi.org/10.1093/nar/gkab335))
* ABRicate ([site](https://github.com/tseemann/abricate))
* AMRFinderPlus ([site](https://github.com/ncbi/amr) | [paper](https://doi.org/10.1038/s41598-021-91456-0))
* dbCAN ([site](https://github.com/linnabrown/run_dbcan) | [paper](https://doi.org/10.1093/nar/gky418))


## Strainphlan Pipeline
This pipeline is simple and runs StrainPhlAn for each specified species.
It also accepts a manifest file containing a list of paths to the `.pkl` files generated by metaphlan as part of the main pipeline.

### Usage

```sh
snakemake \
  --profile msk-lsf/ \
  --snakefile vdb_shotgun/strainphlan.smk \
  --directory tmpstrain/ \
  --config \
    manifest=${PWD}/vdb_shotgun/strainphlan_manifest.txt \
    metaphlan_pkl=/data/brinkvd/resources/biobakery_workflows_dbs/mpa_v30_CHOCOPhlAn_201901.pkl \
    marker_in_n_samples=5 \
  --dry-run
```

### Outputs
For each input species:

* Multiple sequence alignment of strains detected in samples
* Phylogenetic tree of strains detected in samples

### Workflow

The rule DAG for two example input species looks like this:

<img src="https://github.com/vdblab/vdblab-pipelines/raw/main/vdb_shotgun/images/strainphlan_dag.png" alt="StrainPhlAn Shotgun Pipeline DAG" width="500">

## Databases
Details of database creation can be found in `database-README.md`



## Modular execution
### Preprocessing
```sh
snakemake --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/Snakefile \
    --directory tmpall/   \
	--config \
	  sample=473  \
	  R1=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz] \
	  R2=[/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz] \
      multiqc_config=${PWD}/vdb_shotgun/multiqc_config.yaml nshards=1 \
	dedup_reads=False kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/

```

### Biobakery
```sh
snakemake --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/rules/biobakery.smk \
    --directory tmpbio/ \
	--config \
	sample=473  \
	R1=${PWD}/tmppre/kneaddata/473_knead_paired_2.fastq.gz \
	R2=${PWD}/tmppre/kneaddata/473_knead_paired_2.fastq.gz
```


### Kraken2
```sh
snakemake --profile ${PWD}/msk-lsf/  --snakefile vdb_shotgun/rules/kraken.smk \
    --directory tmpkraken/
	--config \
	    sample=473  \
	R1=${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz \
	R2=${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz \
	kraken2_db=/data/brinkvd/watersn/minikraken2_v2_8GB_201904_UPDATE/
```

### Assembly
```sh
snakemake --profile ${PWD}/msk-lsf/ --singularity-args "-B ${PWD}" \
    --snakefile vdb_shotgun/rules/assembly.smk  --directory tmpassembly/ \
	--config sample=473 \
	R1=${PWD}/tmppre/trimmed/473_shard001_trim_R1.fastq.gz \
	R2=${PWD}/tmppre/trimmed/473_shard001_trim_R2.fastq.gz

```

### MultiQC
Just run MultiQC on a directory, no need to Snakemake
```sh

cp -r tmppre/reports tmpreports
ccp tmpassembly/quast/quast_473/report.tsv ./tmpreports/
ver="v1.12"
docker run -V $PWD:$PWD docker://ewels/multiqc:${ver} multiqc \
    --config vdb_shotgun/multiqc_config.yaml --force \
	--title "a multiqc report for some test data" \
	-b "generated by ${ver}" --filename multiqc_report.html \
	reports/ --interactive
```


### RGI
```sh
snakemake --profile ${PWD}/msk-lsf/ --singularity-args "--bind $PWD" \
    --snakefile vdb_shotgun/rules/RGI.smk  --directory tmprgi/ \
	--config sample=473 \
	R1=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz \
	R2=$PWD/.test/shotgun/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz

```

### Host-depletion a la https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000393

```
snakemake --profile ${PWD}/msk-lsf-lite/  --snakefile vdb_shotgun/rules/hostdeplete.smk \
    --directory tmphost/   --config sample=473 \
	R1=/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R1_001.fastq.gz \
	R2=/data/brinkvd/data/shotgun/test/473/473_IGO_12587_1_S132_L003_R2_001.fastq.gz \
	human_index_base=/data/brinkvd/resources/bowtie_indexes/chm13.draft_v1.0_plusY/chm13.draft_v1.0_plusY \
	mouse_index_base=/data/brinkvd/resources/bowtie_indexes/GRCm39/GRCm39
```
