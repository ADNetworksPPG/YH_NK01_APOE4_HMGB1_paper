#!/usr/bin/env Rscript

#load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

#set all directory paths
basedir <- paste0("/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/",
                  "results/HMGB1_snRNA_mm10/")

indir <- paste0(basedir,"05_cellranger_count_batch1_batch2_merged/")
outdir <- paste0(basedir,"06_seurat_analysis_batch1_batch2_merged")
path_postfix_filtered_matrix <- "/filtered_feature_bc_matrix"

#create the output directory if it doesn't exist
if(!dir.exists(file.path(outdir, "plot/01_qc"))){
  dir.create(file.path(outdir, "plot/01_qc"), recursive = T)
  dir.create(file.path(outdir, "data/01_qc"), recursive = T)
}
setwd(outdir)

#get the list of all samples
samples_list <- list.dirs(indir,full.names=FALSE,recursive=FALSE)
#read in the sample metadata
meta.data <- read.csv(paste0("../../../data/HMGB1_snRNA_mm10/",
                             "analysis_summary_batch1_batch2_merged.csv"))

#create a data frame to record cell count before and after filtering per sample
cells_per_sample <- data.frame(sample=character(), 
                               filter_stage=character(), 
                               number_of_cells=numeric(),
                               nFeature_RNA_cutoff=numeric(),
                               percent.mt_cutoff=numeric(),
                               stringsAsFactors = FALSE)

#read in data for each sample and perform QC
for(i in 1:length(samples_list)){
  split.obj.name <- strsplit(samples_list[i],"and")[[1]][1]
  datadir <- paste0(indir,samples_list[i],path_postfix_filtered_matrix)
  obj.data <- Read10X(data.dir = datadir)
  this.obj <- CreateSeuratObject(counts = obj.data)
  cells_pre_filter <- ncol(this.obj)
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(samples_list[i],
                                                   "pre_filter",
                                                   cells_pre_filter,
                                                   0,
                                                   0)
  this.obj[["percent.mt"]] <- PercentageFeatureSet(this.obj, pattern = "^mt-")
  
  #plot the distribution of QC metrics
  pdf(paste0(outdir,"/plot/01_qc/",samples_list[i],"_pre_qc_plot.pdf"))
  # Visualize QC metrics as a violin plot
  print(VlnPlot(this.obj, 
                features = c("nFeature_RNA", 
                             "nCount_RNA",
                             "percent.mt"),
                ncol = 3))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "nFeature_RNA"))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "percent.mt"))
  dev.off()
  
  #filter top 1% (99th quantile)
  percent.mt.cutoff <- 0.25
  nfeature.cutoff <- quantile(this.obj$nFeature_RNA, 0.99)
  this.obj <- subset(this.obj,
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < nfeature.cutoff & 
                       percent.mt < percent.mt.cutoff )
  
  #add metadata
  this.obj$sample_number <- samples_list[i]
  this_obj_metadata <- meta.data[grep(split.obj.name,meta.data$Sample..),]
  this.obj$genotype <- this_obj_metadata$Genotype
  this.obj$trt <- this_obj_metadata$Treatment
  this.obj$mouse_number <- this_obj_metadata$Mouse..
  this.obj$sex <- this_obj_metadata$Sex
  this.obj$age_dose1 <- this_obj_metadata$Age.at.Dose.1
  this.obj$age_perfused <- this_obj_metadata$Age.Perfused
  this.obj$date_of_nuclear_isolation <- this_obj_metadata$Date.of.sn.Isolation
  this.obj$sequencing_batch <- this_obj_metadata$Batch
  
  #get number of cells post filtering
  cells_post_filter <- ncol(this.obj)
  
  print(samples_list[i])
  print(head(this.obj[[]]))
  
  cells_per_sample[nrow(cells_per_sample)+1,] <- c(samples_list[i],
                                                   "post_filter",
                                                   cells_post_filter,
                                                   nfeature.cutoff,
                                                   percent.mt.cutoff)
  
  #plot the distribution of QC metrics
  pdf(paste0(outdir,"/plot/01_qc/",samples_list[i],"_post_qc_plot.pdf"))
  # Visualize QC metrics as a violin plot
  print(VlnPlot(this.obj, 
                features = c("nFeature_RNA", 
                             "nCount_RNA",
                             "percent.mt"),
                ncol = 3))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "nFeature_RNA"))
  print(FeatureScatter(this.obj, 
                       feature1 = "nCount_RNA", 
                       feature2 = "percent.mt"))
  dev.off()
  
  #save the filtered Seurat object
  saveRDS(this.obj, file = paste0(outdir,"/data/01_qc/",samples_list[i],".rds"))
  
}

#save the qc metadata
write.csv(cells_per_sample, 
          file = paste0(outdir,"/data/01_qc/qc_per_sample_metadata.csv"), 
          row.names = FALSE)

#create a plot for cell counts before and after filtering for each sample
pdf(paste0(outdir,
           "/plot/01_qc/",
           "qc_cells_per_sample_plot.pdf"),
    width = 14)
print(ggplot(cells_per_sample, 
             aes(factor(sample), 
                 as.numeric(number_of_cells), 
                 fill = rev(filter_stage))) + 
        geom_bar(stat="identity", 
                 position = "dodge") + 
        scale_fill_brewer(palette = "Set1",
                          labels = c("Pre-QC", "Post-QC"))+
        guides(fill=guide_legend(title="QC stage")) +
        xlab("Sample id") +
        ylab("Number of cells") +
        ggtitle("Cells per sample - Pre and Post QC") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
)
dev.off()


#save the session info
writeLines(capture.output(sessionInfo()), 
           file.path(outdir, "sessionInfo_step01_qc.txt"))

#record logs
print("***** Script completed! *****")


########################## END ##########################                      

