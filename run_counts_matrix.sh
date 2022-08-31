#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --job-name=teaching_pipeline
#SBATCH --output=/scratch/users/%u/teaching/logs/run_counts_matrix.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_counts_matrix.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# path variables - do not edit
base_dir=/scratch/users/k2142172/teaching
user_dir=/scratch/users/${kcl_id}
featurecounts=${base_dir}/packages/subread-2.0.1-Linux-x86_64/bin/featureCounts
#gtf=${base_dir}/resources/*.gtf
gtf=/scratch/users/k2142172/resources/GRCh38/*.gtf

# create output dir if necessary, and redirect log and err files there
mkdir -p ${user_dir}/teaching/outputs/gene_expression

# create 'bams' variable containing paths to each bam file
bams=${base_dir}/outputs/bams/*.bam

# run featurecounts over all bams at once to give one gene counts matrix
# matrix will be table-like, genes x samples where any cell is the read counts...
# ...for that gene in that sample
$featurecounts \
  -a $gtf \
  -F GTF \
  -g gene_id \
  -T 8 \
  --verbose \
  -o ${user_dir}/teaching/outputs/gene_expression/gene_counts.tab \
  $bams
