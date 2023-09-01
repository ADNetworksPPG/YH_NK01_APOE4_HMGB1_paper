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
if(!dir.exists(file.path(outdir, "plot/07_subcluster_microglia_clusters_7_26_30"))){
  dir.create(file.path(outdir, "plot/07_subcluster_microglia_clusters_7_26_30"), recursive = T)
  dir.create(file.path(outdir, "data/07_subcluster_microglia_clusters_7_26_30"), recursive = T)
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

#sub cluster the microglia clusters 7, 26 and 30
dat_microglia <- subset(x = dat, idents = c(7,26,30))

#normalization
dat_microglia <- SCTransform(dat_microglia, 
                             method="glmGamPoi")
#PCA
dat_microglia <- RunPCA(dat_microglia, 
                        assay = "SCT", 
                        verbose = FALSE, 
                        approx=FALSE)
pdf(paste0(outdir,
           "/plot/07_subcluster_microglia_clusters_7_26_30/",
           "microglia_data_HMGBi_batch1_batch2_pcaplot.pdf"))
print(ElbowPlot(dat_microglia))
print(ElbowPlot(dat_microglia, ndims = 50))
dev.off()

#identify the list of marker genes for microglia 
marker_genes_microglia <- read.csv(paste0("/gladstone/bioinformatics/",
                                          "adnetworksppg/Project_1_Huang/NK01/",
                                          "data/snRNA_mm10/",
                                          "Hippocampus_cellcluster_MarkerGenes.csv"))
marker_genes_microglia <- 
  marker_genes_microglia$marker_gene_symbol[
    marker_genes_microglia$hippocampal_region == "Microglia"]
marker_genes_microglia <- c(marker_genes_microglia,
                            "Apoe","hapoE-transgene",
                            "Mapt","Human-MAPT")
marker_genes_homeostatic_microglia <- c("P2ry12","Csf1r","Hexb","Cst3",
                                        "Cx3cr1","Siglech","Tgfbr1","Selplg",
                                        "Mef2a","Serinc3")
marker_genes_dam <- c("Cd9","Fth1","Plp1")
marker_genes_erbb_pathway <- c("Nrg1","Nrg2","Nrg3","Camk2a",
                               "Akt3","Ptk2","Mapk8","Erbb4")
marker_genes_alzheimers_disease <- c("Grin1","Grin2a","Grin2b","Itpr1",
                                     "Itpr2","Ryr2","Ryr3","Plcb1",
                                     "Plcb4","Cacna1c","Ppp3ca")
marker_genes_microglia <- unique(c(marker_genes_microglia,
                                   marker_genes_homeostatic_microglia,
                                   marker_genes_dam,
                                   marker_genes_erbb_pathway,
                                   marker_genes_alzheimers_disease))

#set the number of PCs and resolution to be used for clustering as recommended by Yadong
micro_pca_dim <- 15
micro_cluster_res <- 0.9

#perform clustering and UMAP 
dat_microglia_processed <- FindNeighbors(dat_microglia, 
                                         dims = 1:micro_pca_dim)
dat_microglia_processed <- RunUMAP(dat_microglia_processed, 
                                   dims = 1:micro_pca_dim)
dat_microglia_processed <- FindClusters(dat_microglia_processed, 
                                        resolution =  micro_cluster_res)

#rename the cluster ids to start from 1 instead of 0
dat_microglia_processed$orig_seurat_clusters <- dat_microglia_processed$seurat_clusters
dat_microglia_processed$seurat_clusters <- factor(as.numeric(
  as.character(dat_microglia_processed$orig_seurat_clusters))+1)
Idents(dat_microglia_processed) <- dat_microglia_processed$seurat_clusters

#visualizations
#i. #Generate UMAP without cluster labels
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/07_subcluster_microglia_clusters_7_26_30/",
                                 "microglia_data_HMGBi_batch1_batch2", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#ii. Generate UMAP with cluster labels
DimPlot(dat_microglia_processed, 
        pt.size = 1,
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap") %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/07_subcluster_microglia_clusters_7_26_30/",
                                 "microglia_data_HMGBi_batch1_batch2", 
                                 "_pcadims_",
                                 micro_pca_dim,
                                 "_res_",
                                 micro_cluster_res,
                                 "_labeled_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#iii. Generate UMAP split by genotypes
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_genotype_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_microglia_processed, 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype", 
              group.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#iv. Generate UMAP split by treatment
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_treatment_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_microglia_processed, 
              group.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "trt", 
              group.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#v. Generate UMAP split by genotype and treatment
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_genotype_and_trt_umap.pdf"),
    width = 15,
    height = 7)
print(DimPlot(dat_microglia_processed, 
              group.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype_trt", 
              group.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE))
print(DimPlot(dat_microglia_processed, 
              split.by = "genotype_trt",
              pt.size = 1,
              raster=FALSE,
              order = TRUE,
              label = TRUE))
dev.off()

#vi. Generate featureplot of all microglia marker genes
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features = marker_genes_microglia[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_labeled_featureplot.pdf"),
    width = 14,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features = marker_genes_microglia[f],
                    pt.size = 1.5,
                    raster=FALSE,
                    order = TRUE,
                    label = TRUE))
}
dev.off()
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_split_genotype_and_trt_featureplot.pdf"),
    width = 25,
    height = 7)
for(f in 1:length(marker_genes_microglia)){
  print(FeaturePlot(dat_microglia_processed, 
                    features= marker_genes_microglia[f], 
                    pt.size = 1.5,
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE))
  print(FeaturePlot(dat_microglia_processed, 
                    features= marker_genes_microglia[f],
                    pt.size = 1.5,
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE))
}
dev.off()

#vii. Generate dotplot of all microglia marker genes
pdf(file = paste0(outdir,
                  "/plot/07_subcluster_microglia_clusters_7_26_30/",
                  "microglia_data_HMGBi_batch1_batch2", 
                  "_pcadims_",
                  micro_pca_dim,
                  "_res_",
                  micro_cluster_res,
                  "_dotplot.pdf"),
    width = 14,
    height = 7)
print(DotPlot(dat_microglia_processed, 
              features = marker_genes_microglia
))
print(DotPlot(dat_microglia_processed, 
              features = marker_genes_microglia, 
              group.by = "genotype_trt"
))
dev.off()

#viii. #Generate visualizations using the various metadata
meta_features <- c("sample_number", "sex", "age_dose1","age_perfused", 
                   "date_of_nuclear_isolation", "sequencing_batch")
for(mf in meta_features){
  DimPlot(dat_microglia_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE,
          group.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/07_subcluster_microglia_clusters_7_26_30/",
                                   "microglia_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   micro_pca_dim,
                                   "_res_",
                                   micro_cluster_res, "_",
                                   mf, "_umap.pdf")),
           plot = .,
           width = 12,
           height = 7)
  DimPlot(dat_microglia_processed, 
          pt.size = 1,
          raster = FALSE, 
          order = TRUE, 
          label = FALSE, 
          group.by = mf,
          split.by = mf, 
          reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/07_subcluster_microglia_clusters_7_26_30/",
                                   "microglia_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   micro_pca_dim,
                                   "_res_",
                                   micro_cluster_res, "_",
                                   mf, "_split_umap.pdf")),
           plot = .,
           width = 30,
           height = 10)
}

#ix. # Generate feature plot for QC data
metadata_cols_fp <- c("nCount_RNA","nFeature_RNA", "percent.mt")
for(metacol in metadata_cols_fp){
  FeaturePlot(dat_microglia_processed, 
              pt.size = 1.5,
              raster=FALSE,
              order = TRUE, 
              label = FALSE, 
              features=metacol, 
              reduction = "umap") %>%
    ggsave(file = file.path(outdir, 
                            paste0("plot/07_subcluster_microglia_clusters_7_26_30/",
                                   "microglia_data_HMGBi_batch1_batch2", 
                                   "_pcadims_",
                                   micro_pca_dim,
                                   "_res_",
                                   micro_cluster_res, "_", 
                                   metacol, "_featureplot.pdf")),
           plot = .,
           width = 12,
           height = 7)
}


#save the microglia processed sub-clustered data
saveRDS(dat_microglia_processed, 
        file = paste0(outdir,
                      "/data/07_subcluster_microglia_clusters_7_26_30/",
                      "microglia_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                      micro_pca_dim,
                      "_res_",
                      micro_cluster_res,
                      ".rds"))

#Find differentially expressed genes for ID'ing the clusters
MarkersRes <- FindAllMarkers(dat_microglia_processed, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("/data/07_subcluster_microglia_clusters_7_26_30/",
                                  "marker_genes_microglia_data_post_subcluster",
                                  "_HMGBi_batch1_batch2_pcadims_",
                                  micro_pca_dim,
                                  "_res_",
                                  micro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


#Find differentially expressed genes for ID'ing the clusters
#no marker genes detected for subclusters 2 and 5
#re-run FindAllMarkers with default only.pos=FALSE and logfc threshold = 0.25
MarkersRes <- FindAllMarkers(dat_microglia_processed, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             #only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.25)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("/data/07_subcluster_microglia_clusters_7_26_30/",
                                  "marker_genes_only_pos_false_logfc_thresh_0.25_",
                                  "microglia_data_post_subcluster",
                                  "_HMGBi_batch1_batch2_pcadims_",
                                  micro_pca_dim,
                                  "_res_",
                                  micro_cluster_res,
                                  ".csv")),
          row.names = FALSE)


#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step07_subcluster_microglia_clusters_7_26_30.txt"))

#record logs
print("*********** Script completed! ***********")

############### END ###############

