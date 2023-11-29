
suppressPackageStartupMessages({
  library(rstudioapi)
  library(DESeq2)
  library(clusterProfiler, quietly = TRUE)
  library(msigdbr, quietly = T)
  library(UpSetR, quietly = TRUE)
  library(enrichplot, quietly = TRUE)
  library(ggplot2)
  library(data.table)
})
#Paths
workingD <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(workingD))

#Parameters
#Which database inside msigdbr?
category <- 'C5'
subcategory <- 'HPO'
#msigdbr_collections()

#Input
experiment <- 'Salmon_polyA_filteredClda'
input <- paste0('../original_files/DeSeq_results/', experiment, '.txt')

#Outputs
resD0 <- paste0(experiment, '/')

resD <- gsub(':','_',paste0(resD0,category,'_', subcategory, '/'))
#resD <- paste0(dirname(workingD), resD)
if (!file.exists(resD)){
  dir.create(file.path(resD))
}



resTSV <- paste0(resD,'GSEA_results',gsub(':','_',paste0(category,'_', subcategory)), '.txt')

#1) Load data
data <- fread(input)
#duplicates_count <- data %>% group_by(ENSID) %>% summarise(count = n())

#I am going to work over the stat parameter. 
#To calculate this i divide the log2FoldChange by the lfcSE.
#The other option, using log2FoldChange, would require shrinkage
#As the data is not a DESEqDataSet object, I cannot do the shrinkage. 
#Stat not always available, check input
original_gene_list <- data$stat
names(original_gene_list)<- data$GENE_ID
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing=T)

#2) Calculate GSEA and write tables of results
#Get genes and categories
db_sets <- msigdbr(species = 'Mus musculus', category = category, 
                   subcategory = subcategory)%>% 
  dplyr::select(gs_name, ensembl_gene)
head(db_sets) #each gene associated with each msig group

#3) Perform GSEA
set.seed(1)
egs <- GSEA(geneList = gene_list, pvalueCutoff = 0.1, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
#head(egs@result)
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = resTSV, sep= "\t", quote = F, row.names = T)

