#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                  "results/HMGB1_snRNA_mm10/06_seurat_analysis_batch1_batch2_merged/data")
indir <- paste0(basedir,
                "/06_subcluster_astrocyte_clusters_13_25/")
outdir <- paste0(basedir,
                 "/12_subclusters_de_genes/astrocyte/")

#create the output directory if it doesn't exist
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}
setwd(outdir)

#load the Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astro <- readRDS(paste0(indir,
                            "astrocyte_data_post_subcluster_",
                            "HMGBi_batch1_batch2_pcadims_",
                            pca_dim_astro,
                            "_res_",
                            cluster_res_astro,
                            ".rds"))

# 1. DE gene for:
## 1) Astrocyte subcluster 5 vs other astrocyte subclusters excluding 6, 13, 17.
## 2) For astrocyte subcluster 5, PS19-fE4-HMGB1-i vs PS19-fE4-saline.
## 3) For astrocyte subcluster 5, PS19-fE3-HMGB1-i vs PS19-fE3-saline.
## 4) Astrocyte subcluster 6 vs other astrocyte subclusters excluding 5, 13, 17.
## 5) Astrocyte subcluster 13 vs other astrocyte subclusters excluding 5, 6, 17.
## 6) For astrocyte subcluster 13, PS19-fE4-HMGB1-i vs PS19-fE4-saline.
## 7) For astrocyte subcluster 13, PS19-fE3-HMGB1-i vs PS19-fE3-saline.
## 8) Astrocyte subcluster 17 vs other astrocyte subclusters excluding 5, 6, 13.
## 9) For astrocyte subcluster 17, PS19-fE4-HMGB1-i vs PS19-fE4-saline.
## 10) For astrocyte subcluster 17, PS19-fE3-HMGB1-i vs PS19-fE3-saline.

clusters_for_de <- c(5,6,13,17)
other_clusters <- seq(1, max(as.numeric(dat_astro$seurat_clusters)))
other_clusters <- other_clusters[!(other_clusters %in% clusters_for_de)]
print(other_clusters)

for(de_clus in clusters_for_de){
  de_list <- FindMarkers(object = dat_astro, 
                         ident.1 = de_clus,
                         ident.2 = other_clusters,
                         assay = "SCT", 
                         slot = "data", 
                         test.use = "wilcox",
                         logfc.threshold = 0.1)
  #add gene symbols as a column
  de_list <- cbind(gene=rownames(de_list), 
                   de_list)
  
  #write out the results to a csv file
  write.csv(de_list,
            file = paste0("de_genes_astrocyte_subcluster_",
                          de_clus,
                          "_vs_other_astrocytes_excluding_subclusters_",
                          "5_6_13_17_",
                          "HMGBi_batch1_batch2_pcadims_",
                          pca_dim_astro,
                          "_res_",
                          cluster_res_astro,
                          ".csv"),
            row.names = FALSE)

  
  if(de_clus %in% c(5,13,17)){
    ## fE4-HMGB1i vs fE4-saline
    de_cluster_e4hmgbi_vs_e4saline <- FindMarkers(object = dat_astro,
                                                  ident.1 = "PS19-fE4_HMGB1 Inhibitors",
                                                  ident.2 = "PS19-fE4_Saline",
                                                  verbose = TRUE,
                                                  group.by="genotype_trt",
                                                  subset.ident = de_clus,
                                                  assay = "SCT",
                                                  slot = "data",
                                                  test.use = "wilcox",
                                                  logfc.threshold = 0.1)
    #add gene symbols as a column
    de_cluster_e4hmgbi_vs_e4saline <- cbind(gene=rownames(de_cluster_e4hmgbi_vs_e4saline),
                                            de_cluster_e4hmgbi_vs_e4saline)
    #write out the results to a csv file
    write.csv(de_cluster_e4hmgbi_vs_e4saline,
              file = paste0("de_genes_astrocyte_subcluster_",
                            de_clus,
                            "_fE4-HMGB1i_vs_fE4-saline_",
                            "HMGBi_batch1_batch2_pcadims_",
                            pca_dim_astro,
                            "_res_",
                            cluster_res_astro,
                            ".csv"),
              row.names = FALSE)
  }
  
  if(de_clus %in% c(5,17)){
    ## fE3-HMGB1i vs fE3-saline
    de_cluster_e3hmgbi_vs_e3saline <- FindMarkers(object = dat_astro,
                                                  ident.1 = "PS19-fE3_HMGB1 Inhibitors",
                                                  ident.2 = "PS19-fE3_Saline",
                                                  verbose = TRUE,
                                                  group.by="genotype_trt",
                                                  subset.ident = de_clus,
                                                  assay = "SCT",
                                                  slot = "data",
                                                  test.use = "wilcox",
                                                  logfc.threshold = 0.1)
    #add gene symbols as a column
    de_cluster_e3hmgbi_vs_e3saline <- cbind(gene=rownames(de_cluster_e3hmgbi_vs_e3saline),
                                            de_cluster_e3hmgbi_vs_e3saline)
    #write out the results to a csv file
    write.csv(de_cluster_e3hmgbi_vs_e3saline,
              file = paste0("de_genes_astrocyte_subcluster_",
                            de_clus,
                            "_fE3-HMGB1i_vs_fE3-saline_",
                            "HMGBi_batch1_batch2_pcadims_",
                            pca_dim_astro,
                            "_res_",
                            cluster_res_astro,
                            ".csv"),
              row.names = FALSE)
  }
}



#2. find the list of background genes
all_data <- GetAssayData(dat_astro, assay = "SCT")
background_genes <- rownames(all_data[rowSums(all_data[])>0,])
write.csv(as.data.frame(background_genes),
          file = paste0("nonzero_bakcground_gene_list_",
                        "astrocyte_data_HMGBi_batch1_batch2_pcadims_",
                        pca_dim_astro,
                        "_res_",
                        cluster_res_astro,
                        ".csv"),
          row.names = FALSE)



#save the session info
writeLines(capture.output(sessionInfo()), 
           "../../../sessionInfo_step12_01_subclusters_astrocyte_de_genes.txt")

#record logs
print("********** Script completed! **********")


################## END ################## 

