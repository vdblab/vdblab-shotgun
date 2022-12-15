# VDB Shotgun Pipeline Development

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Apptainer/Singularity](https://apptainer.org/): while in many cases we do provide conda envs the only method of execution we support is via containers.
- (optional) [A Snakemake Profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles): this coordinates the execution of jobs on whatever hardware you are using.
- Fetch test data by running `cd .test/ && bash getdata.sh`.  This script will pull in the `small` and `medium` publicly available datasets.

## About the Data

### `tiny`
The `tiny` dataset is a subset of a patient stool sample; its creation is described in `.test/473/README.md`.  It is used for quickly testing processes, but the small number of reads (<1000) means that it fails metaphlan profiling, binning, and annotation.

### `small`
This is a whole-genome sequencing run of an E. coli isolate deposited on SRA.  It has enough reads to be processed by all stages of the pipeline, but because it is a single isolate its output is not particularly insightful.


### `medium`

This is a stool shotgun sequencing run deposited on SRA.  It is small but closer to real-world data while executing in a reasonable amount of time for development.


## The test script
Currently the script `test.sh` contains logic to run different pipeline stages.

#### USAGE:

```sh
bash test.sh <pipeline_stage> <dataset>
# eg
bash test.sh preprocess tiny
```
The ouput will be in a folder called `tmp<pipeline_stage>_<dataset>`.

## Updating dag figures in README.md

The `test.sh` can also be used to update the figures included in the README by running the following:


```
bash test.sh figs tiny
