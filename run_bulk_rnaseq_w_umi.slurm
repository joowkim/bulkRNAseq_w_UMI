#! /usr/bin/env bash

#partition - defq, bigmem and xtreme
#SBATCH --job-name=bulk-rnaseq
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH -o mrnaseq.%A.o
#SBATCH -e mrnaseq.%A.e

module load nextflow/23.0.4.1
module load singularity/3.8.0

nextflow run bulk_rnaseq_w_umi.nf -c ./bulk_rnaseq_w_umi_conf/run.config -resume -profile slurm

module unload nextflow/23.0.4.1
module unload singularity/3.8.0