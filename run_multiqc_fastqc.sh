#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --job-name=teaching_pipeline
#SBATCH --verbose
#SBATCH --output=/scratch/users/%u/teaching/logs/run_multiqc_fastqc.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_multiqc_fastqc.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# variables - do not edit
base_dir=/scratch/users/k2142172/teaching
user_dir=/scratch/users/${kcl_id}

# create output dir if necessary, and redirect log and err files there
mkdir -p ${user_dir}/teaching/outputs/multiqc

# load multiqc module
#module load apps/multiqc
multiqc=${base_dir}/packages/multiqc

# run multiqc over the fastqc files created and produce html output
$multiqc ${base_dir}/outputs/fastqc/*.zip --ignore ${base_dir}/outputs/fastqc/*.html -n ${user_dir}/teaching/outputs/multiqc/multiqc_fastqc.html




