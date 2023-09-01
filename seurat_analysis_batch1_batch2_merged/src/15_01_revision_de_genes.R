#!/usr/bin/env Rscript

#load required packages
library(Seurat)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                  "results/HMGB1_snRNA_mm10/06_seurat_analysis_batch1_batch2_merged/data")
indir <- paste0(basedir, "/03_clustering/")
outdir <- paste0(basedir, "/15_revision/")

#create the output directory if it doesn't exist
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}
setwd(outdir)

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir, "sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim, "_res_", cluster_res, ".rds"))


###############################################################################
#1. DE pathway (KEGG) enrichment of the general cell cluster 9 versus 
# other Ex Neurons (clusters 1, 5, 6, 11, 17, 19, 20, 22, 24, 27, 29, 31, 33, 34).
###############################################################################
de_cluster_vs_other_ex_neurons <- FindMarkers(object = dat,
                                              ident.1 = 9,
                                              ident.2 = c(1, 5, 6, 11, 17, 19, 20, 22, 
                                                          24, 27, 29, 31, 33, 34),
                                              assay = "SCT",
                                              slot = "data",
                                              test.use = "wilcox",
                                              logfc.threshold = 0.1)
#add gene symbols as a column
de_cluster_vs_other_ex_neurons <- cbind(gene=rownames(de_cluster_vs_other_ex_neurons),
                                        de_cluster_vs_other_ex_neurons)
#write out the results to a csv file
write.csv(de_cluster_vs_other_ex_neurons,
          file = paste0("de_genes_cluster_9_vs_otherExNeurons_post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                        pca_dim, "_res_", cluster_res,".csv"),
          row.names = FALSE)


###############################################################################
#2. DE pathway (KEGG) enrichment of the general cell cluster 24 versus 
# other Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 22, 27, 29, 31, 33, 34).
###############################################################################
de_cluster_vs_other_ex_neurons <- FindMarkers(object = dat,
                                              ident.1 = 24,
                                              ident.2 = c(1, 5, 6, 9, 11, 17, 19, 20, 
                                                          22, 27, 29, 31, 33, 34),
                                              assay = "SCT",
                                              slot = "data",
                                              test.use = "wilcox",
                                              logfc.threshold = 0.1)
#add gene symbols as a column
de_cluster_vs_other_ex_neurons <- cbind(gene=rownames(de_cluster_vs_other_ex_neurons),
                                        de_cluster_vs_other_ex_neurons)
#write out the results to a csv file
write.csv(de_cluster_vs_other_ex_neurons,
          file = paste0("de_genes_cluster_24_vs_otherExNeurons_post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                        pca_dim, "_res_", cluster_res,".csv"),
          row.names = FALSE)


###############################################################################
#3. Volcano plot of DEGs of the general cell cluster 21 versus 
# other Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 22, 24, 27, 29, 31, 33, 34).
###############################################################################
de_cluster_vs_other_ex_neurons <- FindMarkers(object = dat,
                                              ident.1 = 21,
                                              ident.2 = c(1, 5, 6, 9, 11, 17, 19, 20, 22, 
                                                          24, 27, 29, 31, 33, 34),
                                              assay = "SCT",
                                              slot = "data",
                                              test.use = "wilcox",
                                              logfc.threshold = 0.1)
#add gene symbols as a column
de_cluster_vs_other_ex_neurons <- cbind(gene=rownames(de_cluster_vs_other_ex_neurons),
                                        de_cluster_vs_other_ex_neurons)
#write out the results to a csv file
write.csv(de_cluster_vs_other_ex_neurons,
          file = paste0("de_genes_cluster_21_vs_otherExNeurons_post_clustering_sct_HMGBi_batch1_batch2_pcadims_",
                        pca_dim, "_res_", cluster_res,".csv"),
          row.names = FALSE)



#record logs
print("********** Script completed! **********")


################## END ################## 

