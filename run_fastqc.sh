#!/bin/bash

#SBATCH --partition=brc
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --job-name=teaching_pipeline
#SBATCH --verbose
#SBATCH --output=/scratch/groups/wg_translbio/HPC_RNASeq_course/logs/run_fastqc.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_fastqc.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# variables - do not edit
base_dir=/scratch/groups/wg_translbio/HPC_RNASeq_course
user_dir=/scratch/users/${kcl_id}
fastqc=${base_dir}/packages/FastQC/fastqc

# create output dir if necessary, and redirect log and err files there
mkdir -p ${user_dir}/teaching/outputs/fastqc

# load java module - required to run fastqc
module load apps/openjdk

# run fastqc over each fastq file
# fastqc has no option to provide file with list of fastq files to process
# so you can either run the fastqc example command one by one for each fastq file you have...
# ...or use the 'xargs' command to after 'find'ing all the fastq files...
# ...in the given directory
find ${base_dir}/raw_data -name "*fastq.gz" | xargs $fastqc --outdir=${user_dir}/teaching/outputs/fastqc

