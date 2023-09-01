#KEGG Enrichment Analysis of the enriched genes
#this script uses clusterProfiler::enrichKEGG() function

#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(scales)
library(tools)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/10_HMGB1i_snRNAseq_results/",
                  "06_seurat_analysis_batch1_batch2_merged/")
indir <- paste0(basedir, "15_revision/")
outdir <- indir
setwd(indir)


###############################################################################
#1. DE pathway (KEGG) enrichment of the general cell cluster 9 versus 
# other Ex Neurons (clusters 1, 5, 6, 11, 17, 19, 20, 22, 24, 27, 29, 31, 33, 34).

#2. DE pathway (KEGG) enrichment of the general cell cluster 24 versus 
# other Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 22, 27, 29, 31, 33, 34).

#3. DE pathway (KEGG) enrichment of the general cell cluster 21 versus 
# other Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 22, 24, 27, 29, 31, 33, 34).
###############################################################################
de_files <- list.files(pattern = "de_genes")
background_genes <- read.table(
  paste0("../08_de_genes/nonzero_bakcground_gene_list_",
         "post_clustering_sct_HMGBi_batch1_batch2_pcadims_15_res_0.7.csv"),
  sep = ",",
  header = TRUE)

for(i in 1:length(de_files)){
  #these files have DE genes that are significant
  signif_res <- read.table(de_files[i], 
                           sep = ",",
                           header = TRUE)
  colnames(signif_res)[1] <- "Geneid"
  
  #get the gene entrez IDs
  entrezIDs <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keys = background_genes$background_genes%>% 
                                       as.character(),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  entrezIDs_signif <- AnnotationDbi::select(org.Mm.eg.db, 
                                            keys = signif_res$Geneid %>% 
                                              as.character(),
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")
  
  #-----------------------------------
  #KEGG analysis
  #------------------------------
  ekeggbp <- enrichKEGG(
    gene     = entrezIDs_signif$ENTREZID %>% subset(., !is.na(.)),
    universe = entrezIDs$ENTREZID %>% subset(., !is.na(.)),
    organism    = "mmu",
    minGSSize = 10,
    pvalueCutoff = 0.8,
    keyType = "ncbi-geneid"
  )
  
  #translating gene IDs to human readable symbols
  ekeggbp <- setReadable(ekeggbp, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  
  #Visualize
  ## save images
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-dot.pdf"),
      height = 8)
  print(dotplot(ekeggbp, showCategory = 20, orderBy="GeneRatio") )
  dev.off()
  
  x2 <- enrichplot::pairwise_termsim(ekeggbp)
  pdf(paste0(outdir,
             tools::file_path_sans_ext(de_files[i]),
             "_ekegg-emap.pdf"),
      height = 8)
  print(emapplot(x2, showCategory = 40) )
  dev.off()
  
  #save the list of enriched pathways
  write.csv(ekeggbp,file = paste0(outdir,
                                  tools::file_path_sans_ext(de_files[i]),
                                  "_ekegg.csv"))
  
}



###############################################################################
#3. Volcano plot of the general cell cluster 21 versus 
# other Ex Neurons (clusters 1, 5, 6, 9, 11, 17, 19, 20, 22, 24, 27, 29, 31, 33, 34).

#4.Volcano plot of DEGs of general cluster 21 of PS19-E4-HMGB1-inhibitor versus 
#  PS19-E4-Saline (you already generated a DEG list, just need to make a volcano plot).
###############################################################################
#get the list of all the de gene files
de_files <- c("de_genes_cluster_21_vs_otherExNeurons_post_clustering_sct_HMGBi_batch1_batch2_pcadims_15_res_0.7.csv",
              "../08_de_genes/de_genes_cluster_21_fE4-HMGB1i_vs_fE4-saline_post_clustering_sct_HMGBi_batch1_batch2_pcadims_15_res_0.7.csv")


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
  
  #generate the volcano plot with defaults
  #version 1
  pdf(paste0(outdir,
             "volcano_plot_",
             tools::file_path_sans_ext(basename(de_files[i])),
             "_p_cutoff_1e-5_fc_cutoff_0.4.pdf"),
      height = height_selected,
      width = 20)
  print(EnhancedVolcano(signif_res,
                        lab = signif_res$gene,
                        title = tools::file_path_sans_ext(basename(de_files[i])),
                        x = 'avg_log2FC',
                        y = 'p_val',
                        FCcutoff = 0.4,
                        drawConnectors = TRUE,
                        legendLabSize = 10,
                        ylab = bquote(~-Log[10] ~ italic(p-value)),
                        legendLabels = c("NS", expression(Log[2] ~ FC), 
                                         "p-value", 
                                         expression(p-value ~ and
                                                    ~ log[2] ~ FC))
  ))
  dev.off()
  
  #version 2
  pdf(paste0(outdir,
             "volcano_plot_",
             tools::file_path_sans_ext(basename(de_files[i])),
             "_p_cutoff_10e-5_fc_cutoff_0.4.pdf"),
      height = height_selected,
      width = 20)
  print(EnhancedVolcano(signif_res,
                        lab = signif_res$gene,
                        title = tools::file_path_sans_ext(basename(de_files[i])),
                        x = 'avg_log2FC',
                        y = 'p_val',
                        FCcutoff = 0.4,
                        pCutoff = 10e-5,
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



###############################################################################
#Dotplots
###############################################################################
#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(basedir, "03_clustering/sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim, "_res_", cluster_res, ".rds"))

#load the astrocytes Seurat object
pca_dim_astro <- 15
cluster_res_astro <- 0.9
dat_astrocytes <- readRDS(paste0(basedir, "06_subcluster_astrocyte_clusters_13_25/",
                                 "astrocyte_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                                 pca_dim_astro, "_res_", cluster_res_astro, ".rds"))

#load the microglia Seurat object
pca_dim_micro <- 15
cluster_res_micro <- 0.9
dat_microglia <- readRDS(paste0(basedir, "07_subcluster_microglia_clusters_7_26_30/",
                                "microglia_data_post_subcluster_HMGBi_batch1_batch2_pcadims_",
                                pca_dim_micro, "_res_", cluster_res_micro, ".rds"))

#5. Dot plot of the following gene expression of PS19-E4-saline and PS19-E4-HMGB1-inhibitor 
# in general cluster 21. Prkcq, Vmp1, Pik3r1, Dapk1, Deptor, Bcl2, Rras2, Prkacb, Bcl2l1, 
#  Prkaa2, Sh3glb1, Pten, Ctsd, Mapk1, Lamp1, Pik3r3, Mapk8, Mapk10, hapoE, Sirt1.
genes_to_plot <- "Prkcq, Vmp1, Pik3r1, Dapk1, Deptor, Bcl2, Rras2, Prkacb, Bcl2l1, Prkaa2, Sh3glb1, Pten, Ctsd, Mapk1, Lamp1, Pik3r3, Mapk8, Mapk10, hapoE-transgene, Sirt1"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
dat$seurat_cluster_genotype_trt <- paste0("cluster_", dat$seurat_clusters, "_", dat$genotype_trt)
dat_subset <- subset(dat, subset = genotype_trt %in% c("PS19-fE4_Saline", "PS19-fE4_HMGB1 Inhibitors"))
(DotPlot(dat_subset, 
         features = genes_to_plot, 
         #idents = 21,
         group.by = "seurat_cluster_genotype_trt",
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_selected_genes_cluster_21_genotype_trt.pdf",
         plot = .,
         width = 16,
         height = 18)


#6. Dot plot of the following gene expression in astrocyte subclusters 5, 6, 13, 17.
#  Tlr2, Tlr3, Tlr4, Tlr 7, Tlr9, Md2, Rage, Tim3, Il1rl1, Il1rl2 Cxcr4
genes_to_plot <- "Tlr2, Tlr3, Tlr4, Tlr7, rna_Tlr9, Md2, Ly96, Rage, Mok, Ager, Tim3, Havcr2, rna_Il1rl1, Il1rl2, Cxcr4"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_astrocytes, 
         features = genes_to_plot, 
         #idents = c(5, 6, 13, 17),
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_1_selected_genes_astrocyte_subclusters_5_6_13_17.pdf",
         plot = .)


#7. Dot plot of the following gene expression in microglia subclusters 1, 2, 7, 13.
#  Tlr2, Tlr3, Tlr4, Tlr 7, Tlr9, Md2, Rage, Tim3, Il1rl1, Il1rl2, Cxcr4
genes_to_plot <- "Tlr2, Tlr3, Tlr4, Tlr7, Tlr9, Md2, Ly96, Rage, Mok, Ager, Tim3, Havcr2, Il1rl1, Il1rl2, Cxcr4"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_microglia, 
         features = genes_to_plot, 
         #idents = c(1, 2, 7, 13),
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_1_selected_genes_microglia_subclusters_1_2_7_13.pdf",
         plot = .)


#8. Dot plot of the following gene expression in microglia subclusters 1, 2, 7, 13.
#  Ddx58, Trim25, Lyn, Ccl4, Birc3, Plcg2, Icam1, Nfkb1, Cflar, Bcl2l1, Cxcl10, Stat1, 
#  Irf7, Tlr3, Tlr7, Cd86, Mapk10, Pik3r1, Ccl12, Cybb, Plcb1, Prkcd, Stat3, Tgfbr2
genes_to_plot <- "Ddx58, Trim25, Lyn, Ccl4, Birc3, Plcg2, Icam1, Nfkb1, Cflar, Bcl2l1, Cxcl10, Stat1, Irf7, Tlr3, Tlr7, Cd86, Mapk10, Pik3r1, Ccl12, Cybb, Plcb1, Prkcd, Stat3, Tgfbr2"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_microglia, 
         features = genes_to_plot, 
         #idents = c(1, 2, 7, 13),
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_2_selected_genes_microglia_subclusters_1_2_7_13.pdf",
         plot = .,
         width = 20)


#9. Dot plot of the following gene expression in astrocyte subclusters 5, 6, 13, 17.
#  Blnk, Tnfrsf11a, Plcg2, Syk, Lyn, Btk, Cflar, Nfkb1, Bcl2l1, Ddx58, Irak1, Ikbkg, Bcl2, Prkcb, 
#  Xiap, Erc1, Cyld, Cd86, Tlr7, Casp8, Ifnar2, Stat1, Pik3cb, Mapk14, Pik3ca, Tbk1, Mapk10, Map2k1, 
#  Map2k4, Akt3, Tgfbr1, Tgfbr2, Prkcd, Cdc42, Smad3, Mapk14, Smad4, Kras, Stat3, Smad2, Plcb1, Bcl2, 
#  Jak2, Plce1, Prkce, Vegfa, F3
genes_to_plot <- "Blnk, Tnfrsf11a, Plcg2, Syk, Lyn, Btk, Cflar, Nfkb1, Bcl2l1, Ddx58, Irak1, Ikbkg, Bcl2, Prkcb, Xiap, Erc1, Cyld, Cd86, Tlr7, Casp8, Ifnar2, Stat1, Pik3cb, Mapk14, Pik3ca, Tbk1, Mapk10, Map2k1, Map2k4, Akt3, Tgfbr1, Tgfbr2, Prkcd, Cdc42, Smad3, Mapk14, Smad4, Kras, Stat3, Smad2, Plcb1, Bcl2, Jak2, Plce1, Prkce, Vegfa, F3"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
(DotPlot(dat_astrocytes, 
         features = genes_to_plot, 
         #idents = c(5, 6, 13, 17),
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_2_selected_genes_astrocyte_subclusters_5_6_13_17.pdf",
         plot = .,
         width = 20)


#10. Dot plot of the following gene expression of PS19-E4-saline and PS19-E4-HMGB1-inhibitor 
#  in microglia subcluster 1, 2, 7, 13.
#  Tyrobp, Ctsb, hapoE, B2m, Fth1, Lyz2, Trem2, Axl, Cst7, Cst1, Cstl, Lpl, Cd9, Csf1, Ccl6, Itgax, Lilrb4a, Pimp2.
genes_to_plot <- "Tyrobp, Ctsb, hapoE-transgene, B2m, Fth1, Lyz2, Trem2, Axl, Cst7, Cst1, Cstl, Lpl, Cd9, Csf1, Ccl6, Itgax, Lilrb4a, Pimp2"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
dat_microglia$seurat_cluster_genotype_trt <- paste0("microglia_subcluster_", dat_microglia$seurat_clusters, 
                                                    "_", dat_microglia$genotype_trt)
dat_microglia_subset <- subset(dat_microglia, subset = genotype_trt %in% c("PS19-fE4_Saline", "PS19-fE4_HMGB1 Inhibitors"))
(DotPlot(dat_microglia_subset, 
         features = genes_to_plot, 
         #idents = c(1, 2, 7, 13),
         group.by = "seurat_cluster_genotype_trt",
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_3_selected_genes_microglia_subclusters_1_2_7_13_genotype_trt.pdf",
         plot = .,
         width = 20,
         height = 12)


#11. Dot plot of the following gene expression of PS19-E4-Saline and PS19-E4-HMGB1-inhibitor 
#  in microglia subclusters 1, 2, 7, 13.
#  Ddx58, Trim25, Lyn, Ccl4, Birc3, Plcg2, Icam1, Nfkb1, Cflar, Bcl2l1, Cxcl10, Stat1, Irf7, Tlr3, 
#  Tlr7, Cd86, Mapk10, Pik3r1, Ccl12, Cybb, Plcb1, Prkcd, Stat3, Tgfbr2
genes_to_plot <- "Ddx58, Trim25, Lyn, Ccl4, Birc3, Plcg2, Icam1, Nfkb1, Cflar, Bcl2l1, Cxcl10, Stat1, Irf7, Tlr3, Tlr7, Cd86, Mapk10, Pik3r1, Ccl12, Cybb, Plcb1, Prkcd, Stat3, Tgfbr2"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
dat_microglia$seurat_cluster_genotype_trt <- paste0("microglia_subcluster_", dat_microglia$seurat_clusters, 
                                                    "_", dat_microglia$genotype_trt)
dat_microglia_subset <- subset(dat_microglia, subset = genotype_trt %in% c("PS19-fE4_Saline", "PS19-fE4_HMGB1 Inhibitors"))
(DotPlot(dat_microglia_subset, 
         features = genes_to_plot, 
         #idents = c(1, 2, 7, 13),
         group.by = "seurat_cluster_genotype_trt",
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_4_selected_genes_microglia_subclusters_1_2_7_13_genotype_trt.pdf",
         plot = .,
         width = 20)


#12. Dot plot of the following gene expression of PS19-E4-Saline and PS19-E4-HMGB1-inhibitor 
# in astrocyte subclusters 5, 6, 13, 17.
#  Blnk, Tnfrsf11a, Plcg2, Syk, Lyn, Btk, Cflar, Nfkb1, Bcl2l1, Ddx58, Irak1, Ikbkg, Bcl2, Prkcb, 
#  Xiap, Erc1, Cyld, Cd86, Tlr7, Casp8, Ifnar2, Stat1, Pik3cb, Mapk14, Pik3ca, Tbk1, Mapk10, Map2k1, 
#  Map2k4, Akt3, Tgfbr1, Tgfbr2, Prkcd, Cdc42, Smad3, Mapk14, Smad4, Kras, Stat3, Smad2, Plcb1, Bcl2, 
#  Jak2, Plce1, Prkce, Vegfa, F3
genes_to_plot <- "Blnk, Tnfrsf11a, Plcg2, Syk, Lyn, Btk, Cflar, Nfkb1, Bcl2l1, Ddx58, Irak1, Ikbkg, Bcl2, Prkcb, Xiap, Erc1, Cyld, Cd86, Tlr7, Casp8, Ifnar2, Stat1, Pik3cb, Mapk14, Pik3ca, Tbk1, Mapk10, Map2k1, Map2k4, Akt3, Tgfbr1, Tgfbr2, Prkcd, Cdc42, Smad3, Mapk14, Smad4, Kras, Stat3, Smad2, Plcb1, Bcl2, Jak2, Plce1, Prkce, Vegfa, F3"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
dat_astrocytes$seurat_cluster_genotype_trt <- paste0("astrocyte_subcluster_", dat_astrocytes$seurat_clusters, 
                                                     "_", dat_astrocytes$genotype_trt)
dat_astrocytes_subset <- subset(dat_astrocytes, subset = genotype_trt %in% c("PS19-fE4_Saline", "PS19-fE4_HMGB1 Inhibitors"))
(DotPlot(dat_astrocytes_subset, 
         features = genes_to_plot, 
         #idents = c(5, 6, 13, 17),
         group.by = "seurat_cluster_genotype_trt",
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_3_selected_genes_astrocyte_subclusters_5_6_13_17_genotype_trt.pdf",
         plot = .,
         width = 20)

#13.Dotplot of the following gene expression of PS19-E4-saline and PS19-E4-HMGB1-inhibitor in 
# astrocyte subcluster 5, 6, 13, 17.
# Gfap, Id3, Aqp4, Id1, Fabp7, Ctsb, Vim, Osmr, Sepina3n, Ggta1, C1qa, Apod, Cryab, Ptgds, Plp1, 
# Fth1, C4b, Cd9, Sparcl1, Plce1, Sgcd, Fos, S100a6, hAPOE, Clu, Ttr.
genes_to_plot <- "Gfap, Id3, Aqp4, Id1, Fabp7, Ctsb, Vim, Osmr, Sepina3n, Ggta1, C1qa, Apod, Cryab, Ptgds, Plp1, Fth1, C4b, Cd9, Sparcl1, Plce1, Sgcd, Fos, S100a6, hapoE-transgene, Clu, Ttr"
genes_to_plot <- unique(as.character(unlist(strsplit(as.character(genes_to_plot), ", "))))
dat_astrocytes$seurat_cluster_genotype_trt <- paste0("astrocyte_subcluster_", dat_astrocytes$seurat_clusters, 
                                                     "_", dat_astrocytes$genotype_trt)
dat_astrocytes_subset <- subset(dat_astrocytes, subset = genotype_trt %in% c("PS19-fE4_Saline", "PS19-fE4_HMGB1 Inhibitors"))
(DotPlot(dat_astrocytes_subset, 
         features = genes_to_plot, 
         #idents = c(5, 6, 13, 17),
         group.by = "seurat_cluster_genotype_trt",
         cols = c("blue", "red")) + 
    RotatedAxis())%>%
  ggsave(file = "dotplot_4_selected_genes_astrocyte_subclusters_5_6_13_17_genotype_trt.pdf",
         plot = .,
         width = 20)


#save the session info
writeLines(capture.output(sessionInfo()), 
           "../sessionInfo_step15_revision.txt")

#record logs
print("********** Script completed! **********")


################## END ################## 

