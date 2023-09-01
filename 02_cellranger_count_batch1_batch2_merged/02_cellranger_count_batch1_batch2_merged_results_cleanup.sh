#!/bin/bash

##########################################################################
## This script takes two command line arguments (in the below order)
## 1. CellRanger results folder path to be cleaned up
## 2. Temporary or trash folder path
##
## Example run:
## 02_cellranger_count_batch1_batch2_merged_results_cleanup.sh \
## /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/results/HMGB1_snRNA_mm10/05_cellranger_count_batch1_batch2_merged \
## /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/results/HMGB1_snRNA_mm10
##########################################################################

#check the command line arguments
if [ -z "$1" ]
  then
    echo "No argument supplied. Please provide folder path to be cleaned up as a command line argument."
    exit 1	
elif [ -z "$2" ]
	then
    echo "Second argument not supplied. Please provide folder path to the trash folder as a command line argument."
    exit 1
fi

#change the working directory
cd $1

tmp_folder=$2/$(basename $1)
mkdir -p $tmp_folder

#clean up the results folder
for d in *; do
	mkdir $tmp_folder/$d
	mv $d/* $tmp_folder/$d/
	cp -R $tmp_folder/$d/outs/filtered_feature_bc_matrix $d
	cp -R $tmp_folder/$d/outs/raw_feature_bc_matrix $d
	cp $tmp_folder/$d/outs/web_summary.html $d
done
