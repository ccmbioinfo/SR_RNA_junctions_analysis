#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --tmp=100G
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --nodes=1

module load /hpf/largeprojects/ccmbio/yliang/clinical_pipeline/pipeline3.8.0

source /hpf/largeprojects/ccmbio/yliang/clinical_pipeline/python_venv/bin/activate

configfile=/hpf/largeprojects/tgnode/sandbox/yliang_analysis/ashish_analysis/SR_pipeline/hg38config.yaml
sample_file=/hpf/largeprojects/tgnode/sandbox/yliang_analysis/ashish_analysis/SR_pipeline/test.txt
root_dir=/hpf/largeprojects/tgnode/sandbox/yliang_analysis/ashish_analysis/test_pipeline

snakemake -s Snakefile --jobs 1 --configfile $configfile --config sample_file=$sample_file root_dir=$root_dir
