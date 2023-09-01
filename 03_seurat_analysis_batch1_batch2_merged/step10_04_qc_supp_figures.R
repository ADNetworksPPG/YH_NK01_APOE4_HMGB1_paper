#Quality control metrics per cluster
#supplementary figure

#run locally on Ayushi's laptop

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

#set all directory paths
basedir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                  "10_HMGB1i_snRNAseq_results/06_seurat_analysis_batch1_batch2_merged/")
indir <- paste0(basedir, 
                "03_clustering/")
outdir <- indir

#load the Seurat object
pca_dim <- 15
cluster_res <- 0.7
dat <- readRDS(paste0(indir,
                      "sct_data_batch1_batch2_post_cluster_pcadims_",
                      pca_dim,
                      "_res_",
                      cluster_res,
                      ".rds"))

##the number of cells in each mouse in each cluster
all_metadata <- dat[[]]

#make a table of the required data
avg_qc_features_per_cluster <- all_metadata %>% group_by(seurat_clusters) %>% 
  summarise(total_cells = n(),
            avg_ngene = mean(nFeature_RNA),
            sd_ngene = sd(nFeature_RNA),
            se_ngene = sd(nFeature_RNA) / sqrt(length(nFeature_RNA)),
            avg_nUMI = mean(nCount_RNA),
            sd_nUMI = sd(nCount_RNA),
            se_nUMI = sd(nCount_RNA) / sqrt(length(nCount_RNA)),
            avg_mt =  mean(percent.mt),
            sd_mt =  sd(percent.mt),
            se_mt = sd(percent.mt) / sqrt(length(percent.mt))) %>% 
  as.data.frame()

##########
#fig a
##########
# Create plot with legend
ggp1_legend <- ggplot(avg_qc_features_per_cluster, 
                      aes(x=seurat_clusters, y=total_cells, 
                          fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  guides(fill=guide_legend(ncol=4)) + 
  scale_fill_discrete(name = "Cluster Identity", 
                      labels = c("1 - Ex Neuron (GC)",
                                 "2 - Subiculum Neuron",
                                 "3 - Oligodendrocyte",
                                 "4 - Oligodendrocyte",
                                 "5 - Ex Neuron (GC)",
                                 "6 - Ex Neuron (CA1)",
                                 "7 - Microglia",
                                 "8 - In Neuron",
                                 "9 - Ex Neuron (CA2/3)",
                                 "10 - In Neuron",
                                 "11 - Ex Neuron (CA1)",
                                 "12 - In Neuron",
                                 "13 - Astrocyte",
                                 "14 - Subiculum Neuron",
                                 "15 - In Neuron",
                                 "16 - OPC",
                                 "17 - Ex Neuron",
                                 "18 - Oligodendrocyte",
                                 "19 - Subiculum Neuron",
                                 "20 - Subiculum Neuron",
                                 "21 - Ex Neuron",
                                 "22 - Ex Neuron (CA2/3)",
                                 "23 - Unknown",
                                 "24 - Ex Neuron (GC)",
                                 "25 - Astrocyte",
                                 "26 - Microglia",
                                 "27 - Ex Neuron",
                                 "28 - OPC",
                                 "29 - Ex Neuron (CA1)",
                                 "30 - Microglia",
                                 "31 - Ex Neuron (CA1)",
                                 "32 - Choroid Plexus",
                                 "33 - Ex Neuron",
                                 "34 - Ex Neuron (CA2/3)",
                                 "35 - Unknown",
                                 "36 - Unknown",
                                 "37 - Unknown")) +
  theme(text = element_text(size = 26))
# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(ggp1_legend)

##########
#fig b
##########
p1 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=total_cells, fill=seurat_clusters)) +
  geom_bar(position=position_dodge(), stat="identity") + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Number of Cells per Cluster") +
  xlab("Cell Cluster") + 
  ylab("Number of Cells")

##########
#fig c
##########
p2 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=avg_ngene, fill=seurat_clusters)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg_ngene-se_ngene, ymax=avg_ngene+se_ngene), 
                width=.2,
                position=position_dodge(.9)) + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Average Genes per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Average Number of Genes") 

##########
#fig d
##########
p3 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=avg_nUMI, fill=seurat_clusters)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg_nUMI-se_nUMI, ymax=avg_nUMI+se_nUMI), 
                width=.2,
                position=position_dodge(.9)) + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Average nUMI per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Average nUMI")

##########
#fig e
##########
p4 <- ggplot(avg_qc_features_per_cluster, 
             aes(x=seurat_clusters, y=avg_mt, fill=seurat_clusters)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=avg_mt-se_mt, ymax=avg_mt+se_mt), 
                width=.2,
                position=position_dodge(.9)) + 
  scale_x_discrete(expand = c(0.04,0.04))+
  theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    text = element_text(size = 20),
    legend.position = "none") + 
  ggtitle("Average % Mitochondrial Genes per Cells by Cluster") +
  xlab("Cell Cluster") + 
  ylab("Average % Mito Genes")



# Draw plots with shared legend
#with standard error
pdf(paste0(outdir, "quality_control_supplemental_figure_with_SE_bars.pdf"),
    width = 25, height = 26)
grid.arrange(shared_legend,
             arrangeGrob(p1, p2, p3,p4, ncol = 2),
             heights=c(2, 10)) 
dev.off()

print("********** Script completed! **********")

################## END ################## 
