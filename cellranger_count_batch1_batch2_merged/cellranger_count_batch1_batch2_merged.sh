#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -pe smp 4
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=08:00:00
#$ -j yes

#set the paths
base_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01
container_dir=$base_dir/assets/containers
transcriptome_dir=$base_dir/assets/reference_genomes/analysis_hapoe_chr
export SINGULARITY_BINDPATH="$base_dir"

#change the working directory
cd $base_dir/results/HMGB1_snRNA_mm10/05_cellranger_count_batch1_batch2_merged/

#process command line arguments to the script
sample_id=$1
fastq_dir=$2
fastq_dir="${fastq_dir/;/,}"
sample_name=$3

#run cell ranger count
singularity exec $container_dir/cellranger.sif cellranger count \
--id=$sample_id \
--transcriptome=$transcriptome_dir/adppg-mm10-apoe-chr-mapt-chr \
--fastqs=$fastq_dir \
--sample=$sample_name \
--include-introns=true \


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
