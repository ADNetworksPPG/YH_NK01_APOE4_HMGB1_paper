#KEGG Enrichment Analysis of the enriched genes
#this script uses clusterProfiler::enrichKEGG() function

#run locally on Ayushi's laptop

#load the required packages
library(org.Mm.eg.db)
library(clusterProfiler)
library(tidyr)

#set the working directories
indir <- paste0("~/Dropbox (Gladstone)/YH_NK01-Results/",
                "10_HMGB1i_snRNAseq_results/",
                "06_seurat_analysis_batch1_batch2_merged/08_de_genes/")
setwd(indir)
outdir <- "../09_KEGG_enriched_pathways/"
dir.create(outdir, showWarnings = FALSE)


de_files <- list.files(pattern = "de_genes")
background_genes <- read.table(
  paste0("nonzero_bakcground_gene_list_",
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

#save the session info
writeLines(capture.output(sessionInfo()), 
           "../sessionInfo_step09_network_analysis_kegg.txt")

print("********** Script completed! **********")

################## END ################## 
