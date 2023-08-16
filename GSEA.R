## Based on the scripts by Eva Sacristán and Sandra
#González (GENGS CBMSO)

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
category <- 'C2'
subcategory <- 'REACTOME'
#msigdbr_collections()

topCat <- 15


#Input
experiment <- 'Salmon_polyA'
input <- paste0('DeSeq_results/',experiment, '.txt')

#Outputs
resD0 <- paste0('GSEA_results/', experiment, '/')

resD <- gsub(':','_',paste0(resD0,category,'_', subcategory, '/'))
#resD <- paste0(dirname(workingD), resD)
if (!file.exists(resD)){
  dir.create(file.path(resD))
}

resTSV <- paste0(resD,'GSEA_results.txt')
dotplotF <- paste0(resD, "dotplot.jpeg")
geneconceptF <- paste0(resD,"gene_concept_net.jpeg")
ridgeF <- paste0(resD,"GSEA_ridge.jpeg")
upsetF <- paste0(resD,"upset_plot.jpeg")
gseaplotsF <- paste0(resD,'all_gseaplots.jpeg')


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

# Remove items with name '0'
#dat <- dat[names(dat) != 0,] #only human, discuss how to handle duplicates
#dat <- sort(dat$Stat, decreasing=TRUE)

#2) Calculate GSEA and write tables of results
#Get genes and categories
db_sets <- msigdbr(species = 'Mus musculus', category = category, 
                   subcategory = subcategory)%>% 
  dplyr::select(gs_name, ensembl_gene)
head(db_sets) #each gene associated with each msig group

#3) Perform GSEA
set.seed(1)
egs <- GSEA(geneList = gene_list, pvalueCutoff = 0.05, eps = 0, pAdjustMethod = "BH", 
            seed = T, TERM2GENE = db_sets) #for more accurate p value set eps to 0
#head(egs@result)
egs_df <- data.frame(egs@result)
egs_df <- egs_df[, -c(1,2)]

write.table(egs_df, file = resTSV, sep= "\t", quote = F, row.names = T)

#Reconsider the top category number if there are less terms than especified
if (dim(egs_df)[1] < topCat){
  topCat <- dim(egs_df)[1]
}

#4) Plot the results
##Dotplot
jpeg(file = dotplotF, units = 'in', width = 15, height = 10, 
     res = 300)
par(mar = c(2, 2, 2, 5)) 
title <- 'Dot plot with GSEA categories'
dotplot(egs, x = "GeneRatio", color = "p.adjust", showCategory = 15, 
            font.size = 10, title = title)
invisible(dev.off())

##Gene-concept network
jpeg(file = geneconceptF, units = 'in', width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
cnetplot(egs, categorySize="p.adjust", font.size = 15, colorEdge = T)
invisible(dev.off())

##Ridge line plot
jpeg(file = ridgeF, units = 'in', width = 15, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
ridgeplot(egs, fill="p.adjust", orderBy= 'NES', core_enrichment = T, 
          showCategory = topCat)
invisible(dev.off())

##Upset plot (of the 10 first terms)
#Save the genes in each category in a list
genes_top <- as.data.frame(as.factor(head(egs@result$core_enrichment, topCat)))
list_top <- list()
for (i in 1:topCat) {
    list_top[[i]] <- unlist(strsplit(as.character(genes_top[i,1]),split="/"))   
}

#Store all unique gene IDs
uniq <- as.character(unique(names(dat)))
#Get top functions
func_top <- egs$Description[1:topCat]
#Make sparse matrix with 1 for every gene in each category
mat <- matrix(0L, nrow = length(uniq), ncol = length(func_top)) 
for (gene in 1:length(uniq)) {
  for (func in 1:length(func_top)) {
    gen <- uniq[gene]
    if (gen %in% list_top[[func]]) {
      mat[gene,func] =  1
    }}} 

#For all top categories at once:
gseap <- gseaplot2(egs, geneSetID = 1:topCat, pvalue_table = F)
ggsave(gseap, file=gseaplotsF, device = "jpeg", units= "in", 
       height = 15, width = 20)

