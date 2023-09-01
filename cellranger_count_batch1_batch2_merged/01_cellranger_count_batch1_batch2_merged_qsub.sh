#!/bin/bash

#base project path
base_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01

# path to the cell ranger script
cellranger_script=$base_dir/scripts/YH_NK01/HMGB1_snRNA_mm10/cellranger_count_batch1_batch2_merged/cellranger_count_batch1_batch2_merged.sh

#make the results output directory
mkdir $base_dir/results/HMGB1_snRNA_mm10/05_cellranger_count_batch1_batch2_merged

# qsub cell ranger count script for all independent libraries
while IFS="," read -r rec_column1 rec_column2 rec_column3; 
do 
	smp_id=$rec_column1
	fastq_dir=$rec_column2
	smp=$rec_column3
	qsub $cellranger_script $smp_id $fastq_dir $smp
done < <(tail -n +2 $base_dir/data/HMGB1_snRNA_mm10/samplesheet_batch1_batch2_merged.csv)
