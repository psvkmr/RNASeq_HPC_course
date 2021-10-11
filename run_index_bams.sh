#!/bin/bash

#SBATCH --partition=brc
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --job-name=teaching_pipeline
#SBATCH --output=/scratch/groups/wg_translbio/HPC_RNASeq_course/logs/run_index_bams.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_index_bams.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# variables - do not edit
base_dir=/scratch/groups/wg_translbio/HPC_RNASeq_course
user_dir=/scratch/users/${kcl_id}
samtools=${base_dir}/packages/samtools-1.11/bin/samtools

# create output dir if necessary, and redirect log and err files there
mkdir -p ${user_dir}/teaching/outputs/bams/samtools_stats

# create 'bams' variable with path to all bam files
bams=${base_dir}/outputs/bams/*.bam
bams_array=($bams)

# create 'names' array from 'bams' variable to provide name for each output in loop below
names=$(ls $bams | cut -d"/" -f8 | cut -d"." -f1)
names_array=($names)

# run samtools to index, and provide statistics on, every bam file
for i in ${!bams_array[@]}; do
  bam=${bams_array[i]}
  $samtools index -b $bam
  $samtools idxstats $bam > ${user_dir}/teaching/outputs/bams/samtools_stats/${names_array[i]}_idxstats.txt
  $samtools stats $bam > ${user_dir}/teaching/outputs/bams/samtools_stats/${names_array[i]}_stats.txt
done

