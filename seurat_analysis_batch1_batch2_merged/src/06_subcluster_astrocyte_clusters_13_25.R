#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

#set all directory paths
outdir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                 "results/HMGB1_snRNA_mm10/06_seurat_analysis_batch1_batch2_merged")
indir <- paste0(outdir,
                "/data/03_clustering/")
#create the output directory if it doesn't exist
if(!dir.exists(file.path(outdir, "plot/06_subcluster_astrocyte_clusters_13_25"))){
  dir.create(file.path(outdir, "plot/06_subcluster_astrocyte_clusters_13_25"), recursive = T)
  dir.create(file.path(outdir, "data/06_subcluster_astrocyte_clusters_13_25"), recursive = T)
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

#sub cluster the astrocytes clusters 13 and 25
dat_astrocytes <- subset(x = dat, idents = c(13,25))

#normalization
dat_astrocytes <- SCTransform(dat_astrocytes, 
                              method="glmGamPoi")
#PCA
dat_astrocytes <- RunPCA(dat_astrocytes, 
                         assay = "SCT", 
                         verbose = FALSE, 
                         approx=FALSE)
pdf(paste0(outdir,
           "/plot/06_subcluster_astrocyte_clusters_13_25/",
           "astrocyte_data_HMGBi_batch1_batch2_pcaplot.pdf"))
print(ElbowPlot(dat_astrocytes))
print(ElbowPlot(dat_astrocytes, ndims = 50))
dev.off()

#identify the list of marker genes for astrocytes 
marker_genes_astrocytes <- read.csv(paste0("/gladstone/bioinformatics/",
                                           "adnetworksppg/Project_1_Huang/NK01/",
                                           "data/snRNA_mm10/",
                                           "Hippocampus_cellcluster_MarkerGenes.csv"))
marker_genes_astrocytes <- 
  marker_genes_astrocytes$marker_gene_symbol[
    marker_genes_astrocytes$hippocampal_region == "Astrocytes"]
marker_genes_astrocytes <- c(marker_genes_astrocytes,
                             "Apoe","hapoE-transgene",
                             "Mapt","Human-MAPT")

#set the number of PCs and resolution to be used for clustering as recommended by Yadong
astro_pca_dim <- 15
astro_cluster_res <- 0.9

#perform clustering and UMAP 
dat_astrocytes_processed <- FindNeighbors(dat_astrocytes, 
                                          dims = 1:astro_pca_dim)
dat_astrocytes_processed <- RunUMAP(dat_astrocytes_processed, 
                                    dims = 1:astro_pca_dim)
dat_astrocytes_processed <- FindClusters(dat_astrocytes_processed, 
                                         resolution =  astro_cluster_res)

#rename the cluster ids to start from 1 instead of 0
dat_astrocytes_processed$orig_seurat_clusters <- dat_astrocytes_processed$seurat_clusters
dat_astrocytes_processed$seurat_clusters <- factor(as.numeric(
  as.character(dat_astrocytes_processed$orig_seurat_clusters))+1)
Idents(dat_astrocytes_processed) <- dat_astrocytes_processed$seurat_clusters

#visualizations
#i. #Generate UMAP without cluster labels
DimPlot(dat_astrocytes_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                 "astrocyte_data_HMGBi_batch1_batch2", 
                                 "_pcadims_",
                                 astro_pca_dim,
                                 "_res_",
                                 astro_cluster_res,
                                 "_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#ii. Generate UMAP with cluster labels
DimPlot(dat_astrocytes_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                 "astrocyte_data_HMGBi_batch1_batch2", 
                                 "_pcadims_",
                                 astro_pca_dim,
                                 "_res_",
                                 astro_cluster_res,
                                 "_labeled_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#iii. Generate UMAP split by genotypes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_genotype_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_astrocytes_processed, 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype", 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#iv. Generate UMAP split by treatment
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_treatment_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_astrocytes_processed, 
              group.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "trt", 
              group.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#v. Generate UMAP split by genotype and treatment
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_genotype_and_trt_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_astrocytes_processed, 
              group.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype_trt", 
              group.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_astrocytes_processed, 
              split.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#vi. Generate featureplot of all astrocyte marker genes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features = marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_labeled_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features = marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE,
                    label = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_split_genotype_and_trt_featureplot.pdf"),
    width = 25,
    height = 7)
for(f in 1:length(marker_genes_astrocytes)){
  print(FeaturePlot(dat_astrocytes_processed, 
                    features= marker_genes_astrocytes[f], 
                    pt.size = 1.5,
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE))
  print(FeaturePlot(dat_astrocytes_processed, 
                    features= marker_genes_astrocytes[f],
                    pt.size = 1.5,
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE))
}
dev.off()

#vii. Generate dotplot of all astrocyte marker genes
pdf(file = paste0(outdir,
                  "/plot/06_subcluster_astrocyte_clusters_13_25/",
                  "astrocyte_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  astro_pca_dim,
                  "_res_",
                  astro_cluster_res,
                  "_dotplot.pdf"),
    width = 14,
    height = 7)
print(DotPlot(dat_astrocytes_processed, 
              features = marker_genes_astrocytes
))
print(DotPlot(dat_astrocytes_processed, 
              features = marker_genes_astrocytes, 
              group.by = "genotype_trt"
))
dev.off()

#viii. #Generate visualizations using the various metadata
meta_features <- c("sample_number", "sex", "age_dose1","age_perfused", 
                   "date_of_nuclear_isolation", "sequencing_batch")
for(mf in meta_features){
  DimPlot(dat_astrocytes_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE,
          group.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                   "astrocyte_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   astro_pca_dim,
                                   "_res_",
                                   astro_cluster_res, "_",
                                   mf, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  DimPlot(dat_astrocytes_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by = mf,
          split.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                   "astrocyte_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   astro_pca_dim,
                                   "_res_",
                                   astro_cluster_res, "_",
                                   mf, "_split_umap.pdf")),
           plot = .,
           width = 30,
           height = 10)
}
DimPlot(dat_astrocytes_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        group.by = "sample_number",
        split.by = "sample_number", 
        reduction = "umap",
        ncol = 7) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                 "astrocyte_data_HMGBi_batch1_batch2", 
                                 "_pcadims_",
                                 astro_pca_dim,
                                 "_res_",
                                 astro_cluster_res, "_sample_number_split_umap_v2.pdf")),
         plot = .,
         width = 30,
         height = 20)

#ix. # Generate feature plot for QC data
metadata_cols_fp <- c("nCount_RNA","nFeature_RNA", "percent.mt")
for(metacol in metadata_cols_fp){
  FeaturePlot(dat_astrocytes_processed, 
              pt.size = 1.5,
              raster=FALSE,
              order = TRUE, 
              label = FALSE, 
              features=metacol, 
              reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/06_subcluster_astrocyte_clusters_13_25/",
                                   "astrocyte_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   astro_pca_dim,
                                   "_res_",
                                   astro_cluster_res, "_", 
                                   metacol, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}


#save the astrocyte processed sub-clustered data
saveRDS(dat_astrocytes_processed, 
        file = paste0(outdir,
                      "/data/06_subcluster_astrocyte_clusters_13_25/",
                      "astrocyte_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                      astro_pca_dim,
                      "_res_",
                      astro_cluster_res,
                      ".rds"))

#Find differentially expressed genes for ID'ing the clusters
MarkersRes <- FindAllMarkers(dat_astrocytes_processed, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("/data/06_subcluster_astrocyte_clusters_13_25/",
                                  "marker_genes_astrocyte_data_post_subcluster",
                                  "_HMGBi_batch1_batch2_pcadims_",
                                  astro_pca_dim,
                                  "_res_",
                                  astro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step06_subcluster_astrocyte_clusters_13_25.txt"))

#record logs
print("*********** Script completed! ***********")

############### END ###############

