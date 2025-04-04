# Clinical Applications of and Molecular Insights from RNA-Sequencing in a Rare Disease Cohort

This repository contains scripts and codes that were used in the paper Clinical Applications of and Molecular Insights from RNA-Sequencing in a Rare Disease Cohort (link)

**NOTE that the workflow itself is specifically optimized to run under the internal HPC environment. To be able to run it in a local computer, please remove the ```module load``` part in the Snakefile and make sure you have the required tools installed in your path**

# files in this repo

* hg38config.yaml (The yaml file that contains paths to different annotation files and tool)

* submit_snake.sh (the script to submit the snakemake workflow to the HPC slurm)

* Snakefile (The main snakefile that contains the commands for each steps)

* scripts folder (The folder that contains customized QC scripts and the junction analysis script)

* example.txt (the example format of the sample.txt used in the submit_snake.sh)

# RNA-Seq pipeline

## Basic Workflow

The analysis pipeline consists of these proceeses:

1. Alignment using STAR
2. Expression quantitation using RSEM
3. Splicing analysis
4. Several QC steps to ensure sample quality

## GTEx

GTEX data is used as the control for abnormal junctions detection. The GTEX data are processed in the same way as the samples and stored on a sqlite database

## Reference Genome

This pipeline is using Homo_sapiens.GRCh38.104 from Ensembl as reference genome

## software requirements

HPF module dependencies (if you are running it on your local computer, please ensure you have these installed in your path)

* samtools/1.9
* star/2.7.0f
* java/1.8.0_161
* fastqc/0.11.5
* rsem/1.3.3
* R/3.5.1
* python/3.8.0
* bcftools/1.6
* rnaseqc.v2.3.5

R 

* GenomicFeatures
* DBI
* RSQLite
* tidyverse
* reshape2
* optpars

Python

* pandas
* snakemake
* argparse
* warnings
* os
* subprocess
* re


