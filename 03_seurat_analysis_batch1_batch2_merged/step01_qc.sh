#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -pe smp 1
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes
#$ -P neuroppg

#set paths
data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/
script_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/scripts/YH_NK01/HMGB1_snRNA_mm10/seurat_analysis_batch1_batch2_merged
container_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/assets/containers
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/r_sn_rna_seq.sif Rscript $script_dir/src/01_qc.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
