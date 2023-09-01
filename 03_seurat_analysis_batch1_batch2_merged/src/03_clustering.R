#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

#increase RAM usage so the futurized Seurat functions... 
#can access larger global variables
options(future.globals.maxSize = 4000 * 1024^2)

#set all directory paths
outdir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                 "results/HMGB1_snRNA_mm10/06_seurat_analysis_batch1_batch2_merged")

#create the output directory if it doesn't exist
if(!dir.exists(file.path(outdir, "plot/03_clustering"))){
  dir.create(file.path(outdir, "plot/03_clustering"), recursive = T)
  dir.create(file.path(outdir, "data/03_clustering"), recursive = T)
}

#read in the SCTransform normalized dataset
data_sct <- readRDS(paste0(outdir,
                           "/data/02_merge_and_normalize/",
                           "merged_data_batch1_batch2_post_sct.rds"))

# run dimensionality reduction using PCA
number_of_pcs <- length(data_sct$SCT@var.features)
data_sct_pca <- RunPCA(data_sct, 
                       assay = "SCT", 
                       npcs = number_of_pcs, 
                       verbose = FALSE, 
                       approx=FALSE)

pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_PCA_plots.pdf")),
    width = 12,
    height = 7)
print(ElbowPlot(data_sct_pca, ndims = number_of_pcs))
print(ElbowPlot(data_sct_pca, ndims = 50))
print(ElbowPlot(data_sct_pca, ndims = 20))
print(VizDimLoadings(data_sct_pca, dims = 1:2, reduction = "pca"))
print(DimPlot(data_sct_pca, reduction = "pca", raster = FALSE))
print(DimHeatmap(data_sct_pca, dims = 1:2, cells = 500, balanced = TRUE))
print(DimHeatmap(data_sct_pca, dims = 3:4, cells = 500, balanced = TRUE))
dev.off()

print("***** PCA completed! *****")

##################
#Based on email from Yadong...
#use 15 PCs and a resolution of 0.7 for clusteting.
##################
#get the list of marker genes
marker_genes_matrix <- read.csv(paste0("/gladstone/bioinformatics/",
                                       "adnetworksppg/Project_1_Huang/NK01/",
                                       "data/snRNA_mm10/",
                                       "Hippocampus_cellcluster_MarkerGenes.csv"))
features_to_plot <- marker_genes_matrix$marker_gene_symbol
#remove the features that give an error
features_to_plot <- features_to_plot[!(features_to_plot %in% 
                                         c("Rfx3.2","Man1"))]
#additional features
features_to_plot <- c(features_to_plot,"Apoe","hapoE-transgene","Mapt",
                      "Human-MAPT", "Gad1", "Gad2", "Acsbg1", 
                      "S100b", "Cldn10", "Trem2", "Cx3cr1", 
                      "Slc17a7", "Hmgb1", "Nrg1", "Nrg2", 
                      "Nrg3","Erbb4", "Hspa8", "Hsp90ab1", "Hsp90aa1")

pca_dim <- 15
cluster_res <- 0.7

data_sct_cluster <- RunUMAP(data_sct_pca, 
                            dims = 1:pca_dim, 
                            verbose = FALSE)
data_sct_cluster <- FindNeighbors(data_sct_cluster, 
                                  dims = 1:pca_dim, 
                                  verbose = FALSE)
data_sct_cluster <- FindClusters(data_sct_cluster, 
                                 resolution = cluster_res, 
                                 verbose = TRUE)

#rename the cluster ids to start from 1 instead of 0
data_sct_cluster$orig_seurat_clusters <- data_sct_cluster$seurat_clusters
data_sct_cluster$seurat_clusters <- factor(as.numeric(
  as.character(data_sct_cluster$orig_seurat_clusters))+1)
Idents(data_sct_cluster) <- data_sct_cluster$seurat_clusters

#add another metadata column for genotype and treatment
data_sct_cluster$genotype_trt <- paste0(data_sct_cluster$genotype,
                                        "_",
                                        data_sct_cluster$trt)
data_sct_cluster$genotype_trt <- factor(data_sct_cluster$genotype_trt, 
                                        levels = c("PS19-fE3_Saline",
                                                   "PS19-fE3_HMGB1 Inhibitors",
                                                   "PS19-fE4_Saline",
                                                   "PS19-fE4_HMGB1 Inhibitors"))

saveRDS(data_sct_cluster, 
        file = paste0(outdir,
                      "/data/03_clustering/",
                      "sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#visualize the data
#Generate UMAP without cluster labels
DimPlot(data_sct_cluster, 
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap", 
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_clustering/",
                                 "cluster_sct_data_batch1_batch2", 
                                 "_pcadims_",
                                 pca_dim,
                                 "_res_",
                                 cluster_res,
                                 "_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#Generate UMAP with cluster labels
DimPlot(data_sct_cluster, 
        raster = FALSE, 
        order = TRUE, 
        label = TRUE, 
        reduction = "umap", 
        ncol = 2) %>%
  ggsave(file = file.path(outdir, 
                          paste0("plot/03_clustering/",
                                 "cluster_sct_data_batch1_batch2", 
                                 "_pcadims_",
                                 pca_dim,
                                 "_res_",
                                 cluster_res,
                                 "_labeled_umap.pdf")),
         plot = .,
         width = 12,
         height = 7)

#Generate UMAP split by genotypes
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_split_genotype_umap.pdf")),
    width = 25,
    height = 10)
DimPlot(data_sct_cluster, 
        group.by = "genotype",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype", 
        group.by = "genotype",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype",
        raster=FALSE,
        order = TRUE,
        label = TRUE, 
        reduction = "umap")
dev.off()

#Generate UMAP split by treatment
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_split_trt_umap.pdf")),
    width = 25,
    height = 10)
DimPlot(data_sct_cluster, 
        group.by = "trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "trt", 
        group.by = "trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "trt",
        raster=FALSE,
        order = TRUE,
        label = TRUE, 
        reduction = "umap")
dev.off()

#Generate UMAP split by genotype and treatment
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_split_genotype_and_trt_umap.pdf")),
    width = 25,
    height = 10)
DimPlot(data_sct_cluster, 
        group.by = "genotype_trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype_trt", 
        group.by = "genotype_trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype_trt",
        raster = FALSE, 
        order = TRUE, 
        label = FALSE, 
        reduction = "umap")
DimPlot(data_sct_cluster, 
        split.by = "genotype_trt",
        raster=FALSE,
        order = TRUE,
        label = TRUE, 
        reduction = "umap")
dev.off()

#Generate featureplots for all the marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_featureplot.pdf")),
    width = 25,
    height = 15)
for(i in 1:(ifelse(length(features_to_plot)%%2,
                   length(features_to_plot)+1,
                   length(features_to_plot))/8)){
  if((8*i) > length(features_to_plot)){
    print(FeaturePlot(data_sct_cluster, 
                      features = 
                        features_to_plot[((8*(i-1))+1):
                                           (length(features_to_plot))], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = FALSE, 
                      reduction = "umap", 
                      ncol = 4))
  }else{
    print(FeaturePlot(data_sct_cluster, 
                      features= features_to_plot[((8*(i-1))+1):(8*i)], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = FALSE, 
                      reduction = "umap", 
                      ncol = 4))
  }
}
dev.off()

#Generate featureplots with cluster labels for all the marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_labeled_featureplot.pdf")),
    width = 25,
    height = 15)
for(i in 1:(ifelse(length(features_to_plot)%%2,
                   length(features_to_plot)+1,
                   length(features_to_plot))/8)){
  if((8*i) > length(features_to_plot)){
    print(FeaturePlot(data_sct_cluster, 
                      features= 
                        features_to_plot[((8*(i-1))+1):
                                           (length(features_to_plot))], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = TRUE, 
                      reduction = "umap", 
                      ncol = 4))
  }else{
    print(FeaturePlot(data_sct_cluster, 
                      features= features_to_plot[((8*(i-1))+1):(8*i)], 
                      raster = FALSE, 
                      order = TRUE, 
                      label = TRUE, 
                      reduction = "umap", 
                      ncol = 4))
  }
}
dev.off()

#Generate featureplots split by genotype and treatment for all the marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_split_genotype_and_trt_featureplot.pdf")),
    width = 25,
    height = 7)
for(i in 1:length(features_to_plot)){
  print(FeaturePlot(data_sct_cluster, 
                    features= features_to_plot[i], 
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = FALSE, 
                    reduction = "umap"))
  print(FeaturePlot(data_sct_cluster, 
                    features= features_to_plot[i],
                    split.by = "genotype_trt",
                    raster = FALSE, 
                    order = TRUE, 
                    label = TRUE, 
                    reduction = "umap"))
}
dev.off()

#Generate violin plots for all the marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_violinplot.pdf")),
    width = 25,
    height = 15)
for(i in 1:(ifelse(length(features_to_plot)%%2,
                   length(features_to_plot)+1,
                   length(features_to_plot))/8)){
  if((8*i) > length(features_to_plot)){
    print(VlnPlot(data_sct_cluster, 
                  features= 
                    features_to_plot[((8*(i-1))+1):
                                       (length(features_to_plot))],
                  ncol = 4))
  }else{
    print(VlnPlot(data_sct_cluster, 
                  features= features_to_plot[((8*(i-1))+1):(8*i)], 
                  ncol = 4))
  }
}
dev.off()

#Generate violin plots split by genotype and treatment for all the marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_split_genotype_and_trt_violinplot.pdf")),
    width = 25,
    height = 7)
for(i in 1:length(features_to_plot)){
  print(VlnPlot(data_sct_cluster, 
                features= features_to_plot[i], 
                split.by = "genotype_trt"))
  print(VlnPlot(data_sct_cluster, 
                features= features_to_plot[i],
                group.by = "genotype_trt"))
}
dev.off()

#Generate a dot plot for all marker genes of interest
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_dotplot.pdf")),
    width = 20,
    height = 24)
print(DotPlot(data_sct_cluster, 
              features = features_to_plot, 
              cols = c("blue", "red")) + 
        RotatedAxis())   
dev.off()

#Generate a dot plot for all marker genes of interest group by genotype and treatment
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_group_genotype_and_trt_dotplot.pdf")),
    width = 20,
    height = 24)
print(DotPlot(data_sct_cluster, 
              features = features_to_plot, 
              cols = c("blue", "red"),
              group.by = "genotype_trt"
) + 
  RotatedAxis())   
dev.off()

#Make apoE and mapt dotplot for all genotype and treatment groups
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_apoe_mapt_group_genotype_and_trt_dotplot.pdf")),
    width = 7,
    height = 5)
print(DotPlot(data_sct_cluster, 
              features = c("Apoe","hapoE-transgene"), 
              group.by = "genotype_trt"))
print(DotPlot(data_sct_cluster, 
              features = c("Mapt","Human-MAPT"), 
              group.by = "genotype_trt"))
dev.off()

#Histograms for Apoe expression
dat <- data.frame(genotype_trt = factor(data_sct_cluster$genotype_trt),
                  Mouse_apoE = GetAssayData(data_sct_cluster)["Apoe",],
                  Human_apoE = GetAssayData(data_sct_cluster)["hapoE-transgene",])
p1 <- ggplot(dat, 
             aes(x=Mouse_apoE, fill = genotype_trt)) + 
  geom_histogram(position = "identity") +
  scale_x_continuous(breaks = c(0,1,2)) +
  facet_grid(.~genotype_trt) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none")
p2 <- ggplot(dat, 
             aes(x=Human_apoE, fill = genotype_trt)) + 
  geom_histogram(position = "identity") +
  scale_x_continuous(breaks = c(0,1,2)) +
  facet_grid(.~genotype_trt) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none")
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_mouse_and_huma_apoe_expression_",
                            "genotype_and_trt_histogram.pdf")),
    width = 10,
    height = 7)
grid.arrange(p1,p2)
dev.off()

p1 <- ggplot(dat[dat$Mouse_apoE > 0,], 
             aes(x=Mouse_apoE, fill = genotype_trt)) + 
  geom_histogram(position = "identity") +
  scale_x_continuous(breaks = c(0,1,2)) +
  facet_grid(.~genotype_trt) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none")
p2 <- ggplot(dat[dat$Human_apoE > 0,], 
             aes(x=Human_apoE, fill = genotype_trt)) + 
  geom_histogram(position = "identity") +
  scale_x_continuous(breaks = c(0,1,2)) +
  facet_grid(.~genotype_trt) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = "none")
pdf(file = file.path(outdir, 
                     paste0("plot/03_clustering/",
                            "cluster_sct_data_batch1_batch2", 
                            "_pcadims_",
                            pca_dim,
                            "_res_",
                            cluster_res,
                            "_mouse_and_huma_apoe_greater_than_zero_",
                            "expression_genotype_and_trt_histogram.pdf")),
    width = 10,
    height = 7)
grid.arrange(p1,p2)
dev.off()


#Find differentially expressed genes for ID'ing the clusters
MarkersRes <- FindAllMarkers(data_sct_cluster, 
                             assay = "SCT", 
                             slot = "data", 
                             test.use = "wilcox",
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.5)
write.csv(MarkersRes,
          file = file.path(outdir, 
                           paste0("data/03_clustering/",
                                  "marker_genes_post_clustering_sct_",
                                  "batch1_batch2_pcadims_",
                                  pca_dim,
                                  "_res_",
                                  cluster_res,
                                  ".csv")),
          row.names = FALSE)

#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step03_clustering.txt"))

#record logs
print("***** Script completed! *****")


########################## END ########################## 
