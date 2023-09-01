#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                  "results/HMGB1_snRNA_mm10/06_seurat_analysis_batch1_batch2_merged/data")
indir <- paste0(basedir,
                "/03_clustering/")
outdir <- paste0(basedir,
                 "/08_de_genes/")

#create the output directory if it doesn't exist
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#1. DE gene and DE pathway analyses of clusters 9, 14, 21,22, 24, 29 for the following pairs:
#a) fE4-HMGB1i vs fE4-saline
#b) fE3-HMGB1i vs fE3-saline
#c) fE3-saline vs fE4-saline
for(cluster_of_interest in c(9,14,21,22,24,29)){
  ## a) fE4-HMGB1i vs fE4-saline
  de_cluster_e4hmgbi_vs_e4saline <- FindMarkers(object = dat,
                                                ident.1 = "PS19-fE4_HMGB1 Inhibitors",
                                                ident.2 = "PS19-fE4_Saline",
                                                verbose = TRUE,
                                                group.by="genotype_trt",
                                                subset.ident = cluster_of_interest,
                                                assay = "SCT",
                                                slot = "data",
                                                test.use = "wilcox",
                                                logfc.threshold = 0.1)
  #add gene symbols as a column
  de_cluster_e4hmgbi_vs_e4saline <- cbind(gene=rownames(de_cluster_e4hmgbi_vs_e4saline),
                                          de_cluster_e4hmgbi_vs_e4saline)
  #write out the results to a csv file
  write.csv(de_cluster_e4hmgbi_vs_e4saline,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_fE4-HMGB1i_vs_fE4-saline_",
                          "post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  
  
  ## b) fE3-HMGB1i vs fE3-saline
  de_cluster_e3hmgbi_vs_e3saline <- FindMarkers(object = dat,
                                                ident.1 = "PS19-fE3_HMGB1 Inhibitors",
                                                ident.2 = "PS19-fE3_Saline",
                                                verbose = TRUE,
                                                group.by="genotype_trt",
                                                subset.ident = cluster_of_interest,
                                                assay = "SCT",
                                                slot = "data",
                                                test.use = "wilcox",
                                                logfc.threshold = 0.1)
  #add gene symbols as a column
  de_cluster_e3hmgbi_vs_e3saline <- cbind(gene=rownames(de_cluster_e3hmgbi_vs_e3saline),
                                          de_cluster_e3hmgbi_vs_e3saline)
  #write out the results to a csv file
  write.csv(de_cluster_e3hmgbi_vs_e3saline,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_fE3-HMGB1i_vs_fE3-saline_",
                          "post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  
  
  ## c) fE3-saline vs fE4-saline
  de_cluster_e3saline_vs_e4saline <- FindMarkers(object = dat,
                                                 ident.1 = "PS19-fE3_Saline",
                                                 ident.2 = "PS19-fE4_Saline",
                                                 verbose = TRUE,
                                                 group.by="genotype_trt",
                                                 subset.ident = cluster_of_interest,
                                                 assay = "SCT",
                                                 slot = "data",
                                                 test.use = "wilcox",
                                                 logfc.threshold = 0.1)
  #add gene symbols as a column
  de_cluster_e3saline_vs_e4saline <- cbind(gene=rownames(de_cluster_e3saline_vs_e4saline),
                                           de_cluster_e3saline_vs_e4saline)
  #write out the results to a csv file
  write.csv(de_cluster_e3saline_vs_e4saline,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_fE3-saline_vs_fE4-saline_",
                          "post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  
}


#2. find the list of background genes
all_data <- GetAssayData(dat, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                        pca_dim,
                        "_res_",
                        cluster_res,
                        ".csv"),
          row.names = FALSE)

#3. Cluster 14,21 and 22 vs other 
#Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 24, 27, 29, 31, 33, 34)
for(cluster_of_interest in c(14,21,22)){
  de_cluster_vs_other_ex_neurons <- FindMarkers(object = dat,
                                                ident.1 = cluster_of_interest,
                                                ident.2 = c(1, 5, 6, 9, 11, 17,
                                                            19, 20, 24, 27, 29,
                                                            31, 33, 34),
                                                assay = "SCT",
                                                slot = "data",
                                                test.use = "wilcox",
                                                logfc.threshold = 0.1)
  #add gene symbols as a column
  de_cluster_vs_other_ex_neurons <- cbind(gene=rownames(de_cluster_vs_other_ex_neurons),
                                          de_cluster_vs_other_ex_neurons)
  #write out the results to a csv file
  write.csv(de_cluster_vs_other_ex_neurons,
            file = paste0("de_genes_cluster_",
                          cluster_of_interest,
                          "_vs_otherexcitatoryneurons_",
                          "post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                          pca_dim,
                          "_res_",
                          cluster_res,
                          ".csv"),
            row.names = FALSE)
  
}


#save the session info
writeLines(capture.output(sessionInfo()), 
           "../../sessionInfo_step08_de_genes.txt")

#record logs
print("********** Script completed! **********")


################## END ################## 
