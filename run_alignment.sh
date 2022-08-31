#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --job-name=teaching_pipeline
#SBATCH --output=/scratch/users/%u/teaching/logs/run_alignment.out

# CHANGE '--output' ABOVE TO /scratch/users/<your k number>/logs/run_alignment.out
# CHANGE 'kcl_id' BELOW TO YOUR K NUMBER
kcl_id=k2142172

# variables - do not edit
base_dir=/scratch/users/k2142172/teaching
user_dir=/scratch/users/${kcl_id}
star=${base_dir}/packages/STAR-2.7.8a/bin/Linux_x86_64_static/STAR
#star_index=${base_dir}/resources/STAR
star_index=/scratch/users/k2142172/resources/GRCh38/STAR

# create output dir if necessary
mkdir -p ${user_dir}/teaching/outputs/bams

# assign the fastq files to the variable 'fastqs', and create an array
fastqs=${base_dir}/raw_data/*fastq.gz
fastqs_array=($fastqs)

# using the created fastq variable, extract the names to a new variable and create another array
names=$(ls $fastqs | cut -d"/" -f7 | cut -d"." -f1)
names_array=($names)


# load star genome
$star --genomeDir $star_index --genomeLoad LoadAndExit

# run star to create bam alignment files from fastq files
# runs in a for loop for each fastq file in the 'fastqs' variable
# uses the 'names' array to assign the appropriate sample name to each bam file
for i in ${!fastqs_array[@]}; do
  $star \
    --runThreadN 8 \
    --readFilesIn ${fastqs_array[i]} \
    --readFilesCommand zcat \
    --genomeLoad LoadAndKeep \
    --genomeDir $star_index \
    --outFileNamePrefix ${user_dir}/teaching/outputs/bams/${names_array[i]}_ \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 30000000000 \
    --outSAMunmapped Within
done
