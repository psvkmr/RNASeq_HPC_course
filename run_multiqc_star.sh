#!/bin/bash

#SBATCH --partition=brc
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --job-name=teaching_pipeline
#SBATCH --verbose
#SBATCH --output=/scratch/groups/wg_translbio/HPC_RNASeq_course/logs/run_multiqc_star.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_multiqc_star.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# variables - do not edit
base_dir=/scratch/groups/wg_translbio/HPC_RNASeq_course
user_dir=/scratch/users/${kcl_id}

# create output dir if necessary, and redirect log and err files there
mkdir -p ${user_dir}/teaching/outputs/multiqc

# load multiqc
module load apps/multiqc

# run multiqc over star log outputs
multiqc ${base_dir}/teaching/outputs/bams/*Log* --ignore ${base_dir}/teaching/outputs/bams/*bam* -n ${user_dir}/teaching/outputs/multiqc/multiqc_star.html




