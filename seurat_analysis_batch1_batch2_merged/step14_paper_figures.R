#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(scales)
library(tools)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "10_HMGB1i_snRNAseq_results/06_seurat_analysis_batch1_batch2_merged/")
outdir <- paste0(basedir, "14_paper_figures/")
if(!dir.exists(outdir)){
  dir.create(outdir)
}
setwd(outdir)


#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(basedir,
                      "03_clustering/",
                      "sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

#load the astrocytes Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astrocytes <- readRDS(paste0(basedir,
                                 "06_subcluster_astrocyte_clusters_13_25/",
                                 "astrocyte_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                                 pca_dim_astro,
                                 "_res_",
                                 cluster_res_astro,
                                 ".rds"))

#load the microglia Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_microglia <- readRDS(paste0(basedir,
                                "07_subcluster_microglia_clusters_7_26_30/",
                                "microglia_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                                pca_dim_micro,
                                "_res_",
                                cluster_res_micro,
                                ".rds"))


##############################################################
##1. Make a "Genotype-split UMAP with Clusters 7, 9, 21, 24, 
#    26, and 30 Highlighted", following the order and 
#    dimension in the attached Figure 7B
##############################################################
dat$clusters_of_interest <- ifelse(dat$seurat_clusters %in% c(7,9,21,24,26,30),
                                   dat$seurat_clusters,
                                   0)
(DimPlot(dat, 
         split.by = "genotype_trt",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(6)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Cluster 7", 
                                  "Cluster 9", "Cluster 21",
                                  "Cluster 24","Cluster 26", "Cluster 30"),
                       values = c("grey",hue_pal()(6)))) %>%
  ggsave(file = "umap_genotype_trt_split_clusters_of_interest.pdf",
         plot = .,
         width = 12,
         height = 7)


##############################################################
##2. Make a "Genotype-split UMAP with Microglia Subclusters 
#    1, 2, 7, and 13 Highlighted", following the order and 
#    dimension in the attached Figure 7F.
##############################################################
dat_microglia$clusters_of_interest <- ifelse(dat_microglia$seurat_clusters %in% c(1,2,7,13),
                                             dat_microglia$seurat_clusters,
                                             0)
(DimPlot(dat_microglia, 
         split.by = "genotype_trt",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(4)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Subcluster 1", "Subcluster 2", 
                                  "Subcluster 7","Subcluster 13"),
                       values = c("grey",hue_pal()(4)))) %>%
  ggsave(file = "microglia_umap_genotype_trt_split_clusters_of_interest.pdf",
         plot = .,
         width = 12,
         height = 7)


##############################################################
##3. Make a "Genotype-split UMA with Astrocyte Subclusters 
#    5, 6, 13, and 17 Highlighted", following the order and 
#    dimension in the attached Figure 7K.
##############################################################
dat_astrocytes$clusters_of_interest <- ifelse(dat_astrocytes$seurat_clusters %in% c(5,6,13,17),
                                              dat_astrocytes$seurat_clusters,
                                             0)
(DimPlot(dat_astrocytes, 
         split.by = "genotype_trt",
         group.by = "clusters_of_interest",
         raster = FALSE, 
         order = TRUE, 
         cols = c("grey",hue_pal()(4)),
         reduction = "umap") + 
    scale_color_manual(labels = c("Other", "Subcluster 5", "Subcluster 6", 
                                  "Subcluster 13","Subcluster 17"),
                       values = c("grey",hue_pal()(4)))) %>%
  ggsave(file = "astrocyte_umap_genotype_trt_split_clusters_of_interest.pdf",
         plot = .,
         width = 10,
         height = 7)



##############################################################
##4. Make a "Genotype-split APOE Feature Plot" for all cell 
#    clusters (label clusters), following the same order and 
#    dimension as the "Genotype-split UMAP".
##############################################################
#option #1
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2.5) + 
   theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_trt_split_hAPOE_maxcutoff_2.5.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #2
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_trt_split_hAPOE_maxcutoff_2.0.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #3
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) +  
   theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_trt_split_hAPOE_maxcutoff_1.5.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #4
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_trt_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #5
(FeaturePlot(dat, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "featureplot_genotype_trt_split_hAPOE_default.pdf",
         plot = .,
         width = 12,
         height = 7)


##############################################################
##5. Make a Genotype-split APOE Feature Plot" for all microglia 
#    subclusters (label clusters), following the same order 
#    and dimension as the "Genotype-split UMAP".
##############################################################
#option #1
(FeaturePlot(dat_microglia, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "microglia_featureplot_genotype_trt_split_hAPOE_maxcutoff_1.5.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #2
(FeaturePlot(dat_microglia, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "microglia_featureplot_genotype_trt_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 12,
         height = 7)
#option #3
(FeaturePlot(dat_microglia, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "microglia_featureplot_genotype_trt_split_hAPOE_default.pdf",
         plot = .,
         width = 12,
         height = 7)


##############################################################
##6. Make a "Genotype-split APOE Feature Plot" for all 
#    astrocyte subclusters (label clusters), following the 
#    same order and dimension as the "Genotype-split UMAP".
##############################################################
#option #1
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2.5) + 
   theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_trt_split_hAPOE_maxcutoff_2.5.pdf",
         plot = .,
         width = 10,
         height = 7)
#option #2
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 2) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_trt_split_hAPOE_maxcutoff_2.0.pdf",
         plot = .,
         width = 10,
         height = 7)
#option #3
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1.5) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_trt_split_hAPOE_maxcutoff_1.5.pdf",
         plot = .,
         width = 10,
         height = 7)
#option #4
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE,
             max.cutoff = 1) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_trt_split_hAPOE_maxcutoff_1.0.pdf",
         plot = .,
         width = 10,
         height = 7)
#option #5
(FeaturePlot(dat_astrocytes, 
             features = "hapoE-transgene",
             split.by = "genotype_trt",
             raster = FALSE, 
             order = TRUE,
             label = TRUE) + 
    theme(legend.position = "right")) %>%
  ggsave(file = "astrocyte_featureplot_genotype_trt_split_hAPOE_default.pdf",
         plot = .,
         width = 10,
         height = 7)



##############################################################
##7. Volcano plots:
#    1. For cluster 21, make volcano plots of (a) cluster 21 
#       vs other ex neurons and (b) PS19-fE4-HMGB1-i vs PS19-fE4-Saline.
#
#    2. For microglia subclusters 1, 2, 7, and 13, make the 
#    following volcano plots, with the cutoff of 0.4 log2FC and p<0.05 
#    and highlighting top 20 downregulated and top 20 upregulated DE genes (both directions).
#       a) subcluster 1 vs other microglia subclusters and PS19-fE4-HMGB1-i vs PS19-fE4-Saline of subcluster 1.
#       b) subcluster 2 vs other microglia subclusters and PS19-fE4-HMGB1-i vs PS19-fE4-Saline of subcluster 2.
#       c) subcluster 7 vs other microglia subclusters and PS19-fE4-HMGB1-i vs PS19-fE4-Saline of subcluster 7.
#       d) subcluster 13 vs other microglia subclusters and PS19-fE4-HMGB1-i vs PS19-fE4-Saline of subcluster 13.
#
#    3. For astrocyte subclusters 6 and 13, make the 
#    following volcano plots, with the cutoff of 0.4 log2FC and p<0.05 
#    and highlighting top 20 downregulated and top 20 upregulated DE genes 
#    (if there are, for both directions).
#       a) subcluster 6 vs other astrocyte subclusters.
#       b) subcluster 13 vs other astrocyte subclusters.
##############################################################
setwd(basedir)

#get the list of all the de gene files
de_files <- c("08_de_genes/de_genes_cluster_21_vs_otherexcitatoryneurons_post_clustering_sct_HMGBi_batch1_batch2_pcadims_15_res_0.7.csv",
              "08_de_genes/de_genes_cluster_21_fE4-HMGB1i_vs_fE4-saline_post_clustering_sct_HMGBi_batch1_batch2_pcadims_15_res_0.7.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_1_vs_other_microglia_excluding_subclusters_1_2_7_13_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_1_fE4-HMGB1i_vs_fE4-saline_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_2_vs_other_microglia_excluding_subclusters_1_2_7_13_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_2_fE4-HMGB1i_vs_fE4-saline_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_7_vs_other_microglia_excluding_subclusters_1_2_7_13_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_7_fE4-HMGB1i_vs_fE4-saline_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_13_vs_other_microglia_excluding_subclusters_1_2_7_13_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/microglia/de_genes_microglia_subcluster_13_fE4-HMGB1i_vs_fE4-saline_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/astrocyte/de_genes_astrocyte_subcluster_6_vs_other_astrocytes_excluding_subclusters_5_6_13_17_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv",
              "12_subclusters_de_genes/astrocyte/de_genes_astrocyte_subcluster_13_vs_other_astrocytes_excluding_subclusters_5_6_13_17_HMGBi_batch1_batch2_pcadims_15_res_0.9.csv")


############
#volcano plots
############
#generate volcano plot for each of the de gene lists
for(i in 1:length(de_files)){
  #these files have DE genes that are significant
  signif_res <- read.table(de_files[i], 
                           sep = ",",header = TRUE)
  
  #set the height based on the input
  height_selected <- 12
  
  #set the limits for x-axis
  xlim_selected <- c(-3,3)
  if(max(abs(signif_res$avg_log2FC)) < 2){
    xlim_selected <- c(-2,2)
  }
  
  if(grepl("08_de_genes", de_files[i], fixed = TRUE)){
    #generate the volcano plot with defaults
    pdf(paste0(outdir,
               "volcano_plot_",
               tools::file_path_sans_ext(basename(de_files[i])),
               ".pdf"),
        height = height_selected,
        width = 20)
    print(EnhancedVolcano(signif_res,
                          lab = signif_res$gene,
                          title = tools::file_path_sans_ext(basename(de_files[i])),
                          x = 'avg_log2FC',
                          y = 'p_val',
                          xlim = xlim_selected,
                          drawConnectors = TRUE,
                          legendLabSize = 10,
                          ylab = bquote(~-Log[10] ~ italic(p-value)),
                          legendLabels = c("NS", expression(Log[2] ~ FC), 
                                           "p-value", 
                                           expression(p-value ~ and
                                                      ~ log[2] ~ FC))
    ))
    dev.off()
    
  } else{
    #With the cutoff of 0.4 log2FC and p<0.05, label the top 20 downregulated 
    #and top 20 upregulated DE genes (both directions).
    signif_res_to_label_up <- signif_res %>% filter(p_val < 0.05 & avg_log2FC > 0.4) %>% arrange(desc(abs(avg_log2FC))) %>% head(20)
    signif_res_to_label_down <- signif_res %>% filter(p_val < 0.05 & avg_log2FC < -0.4) %>% arrange(desc(abs(avg_log2FC))) %>% head(20)
    signif_res_to_label <- rbind.data.frame(signif_res_to_label_up, signif_res_to_label_down)
    write.csv(signif_res_to_label, 
              file = file.path(outdir, paste0("top_20_upANDdown_genes_",basename(de_files[i]))))
    genes_to_label <- signif_res_to_label$gene
    
    #generate the volcano plot
    pdf(paste0(outdir,
               "volcano_plot_",
               tools::file_path_sans_ext(basename(de_files[i])),
               ".pdf"),
        height = height_selected,
        width = 20)
    print(EnhancedVolcano(signif_res,
                          lab = signif_res$gene,
                          title = tools::file_path_sans_ext(basename(de_files[i])),
                          selectLab = genes_to_label,
                          x = 'avg_log2FC',
                          y = 'p_val',
                          FCcutoff = 0.4,
                          pCutoff = 0.05,
                          max.overlaps = Inf,
                          xlim = xlim_selected,
                          drawConnectors = TRUE,
                          legendLabSize = 10,
                          ylab = bquote(~-Log[10] ~ italic(p-value)),
                          legendLabels = c("NS", expression(Log[2] ~ FC), 
                                           "p-value", 
                                           expression(p-value ~ and
                                                      ~ log[2] ~ FC))
    ))
    dev.off()
  }
}


##############################################################
##8. Dotplot of the select genes for all astrocyte subclusters.
##############################################################
genes_to_plot <- "Luzp2, Slc1a2, Slc1a3, Mfge8, Wdr17, Gpc5, Grm3, Rorb, Rgs20, Tspan7, Ptprt, Nrxn1, Gfap, Id3, Aqp4, Myoc, Id1, Fabp7, Ctsb, Vim, Osmr, Serpina3n, Gsn, Ggta1, Trpm3, Csmd1, C1qa, Apod, Cryab, Ptgds, Plp1, Fth1, C4b, Cd9, Sparcl1, Plce1, Sgcd, Fos, S100a6, hapoE-transgene, Clu, Mertk, Msi2, Npas3, Trpm3, Prex2, Kcnip4, Dpp10, Meg3, Ptprd, Tenm2, Celf2, Snhg11, Nrxn3, Nkain2, Gria1, Syt1, Ttr, Dlg2, Cntnap2, Rbfox1, Opcml, Nrg3, Ubb, Tmsb4x, Calm1, Calm2, Cst3, Nrgn, Camk1d, Hspa8, Hsp90aa1, Hsp90ab1"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_astrocytes, 
         features = genes_to_plot, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes_all_astrocyte_subclusters.pdf",
         plot = .,
         width = 20,
         height = 9)


##############################################################
##9. Dotplot of the select genes for all microglia subclusters.
##############################################################
genes_to_plot <- "Hexb, Cst3, Cx3cr1, Ctsd, Csf1r, Ctss, Sparc, Tmsb4x, P2ry12, C1qa, C1qb, Tmem119, Apod, Ttr, Ptgds, Cryab, Plp1, Enpp2, Ctsd, Tyrobp, Ctsb, hapoE-transgene, B2m, Fth1, Lyz2, Trem2, Axl, Cst3, Cst7, Ctsl, Lpl, Cd9, Csf1, Ccl6, Itgax, Clec7a, Lilrb4a, Timp2, Pde4b, Nkain2, St18, Prr5l, Pcdh9, Dpp10, Kcnip4, Meg3, Csmd1, Csmd3, Nrg1, Nrg3, Nrxn1, Nrxn3, Celf2, Tenm2, Snhg11, Rbfox1, Ptprd, Fgf14, Cntnap2, Dlg2, Cadm2, Nav2, Nav3, Srgap2, Neat1, Slc24a2, Magi2, Plp1, Tmeff2, Frmd5, Ank2, Grik2, Malat1, Acaca, Cmss1, Smyd3, Arsb, Pak7, Tgfbr1, Calm1, Ubb, Tmsb4x, Calm2, Nrgn, Camk1d, Hspa8, Hsp90aa1, Hsp90ab1, Ptprd"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_microglia, 
         features = genes_to_plot, 
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_with_selected_genes_all_microglia_subclusters.pdf",
         plot = .,
         width = 20,
         height = 9)


#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path("sessionInfo_step14_paper_figures.txt"))

print("********** Script completed! **********")

################## END ################## 
