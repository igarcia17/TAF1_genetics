setwd("C:/Users/CBM/Desktop/RNAseq_bipolarDisorder/results_GSEA")
library(DESeq2)
library(GOplot)
input <- '../results_DGE/deseq_objects.RData'
load(input)
colores <- c('red', 'green', 'blue', 'orange', 'yellow', 'pink', 'black')

res <- results(dds, contrast=c("condition", 'Unaffected', 'Affected'))
res$FoldChange <- 2^res$log2FoldChange  #have actual fold change
res <- res[colnames(res)[c(1,7,2:6)]] # order columns

allgenes <- as.data.frame(res$log2FoldChange)
rownames(allgenes) <- rownames(res)
colnames(allgenes) <- 'logFC'
allgenes$gene <- rownames(res)

#Bubble plot of 68 categories validated from GSEA
circ0 <- read.csv('bubble_validated.txt', header = T, sep="\t")

#To compute the Z-score
cats <- as.data.frame(circ0$term)
cats$gene <- circ0$gene
cats$set_size <- circ0$count
colnames(cats) <- c('term', 'gene', 'set_size')

# Assuming you have loaded your data frames 'cats' and 'allgenes'

# Create an empty vector to store Z-scores
z_scores <- numeric(nrow(cats))

# Iterate through each row of 'cats' and calculate the Z-score
for (i in 1:nrow(cats)) {
  # Extract the list of genes in the current term and split it into a vector
  genes_in_term <- unlist(strsplit(cats$gene[i], ","))
  
  # Filter 'allgenes' to get log2 fold changes for genes in the current term
  log2foldchanges_in_term <- allgenes$logFC[allgenes$gene %in% genes_in_term]

  # Calculate the Z-score for the term
  upregulated_count <- sum(log2foldchanges_in_term > 0)
  downregulated_count <- sum(log2foldchanges_in_term < 0)
  set_size <- length(log2foldchanges_in_term)
  
  z_score <- (upregulated_count - downregulated_count) / sqrt(set_size)
  
  # Store the calculated Z-score in the vector
  z_scores[i] <- z_score
}

# Add the Z-scores to the 'cats' data frame
cats$ZScore <- z_scores

#Relate to the circ data frame
circ0$zscore <- cats$ZScore[cats$term == circ0$term]
circ <- circ0

#Z-score P value GSEA

jpeg(file = 'zscore_vs_pGSEA2.jpeg', units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ, labels = 2.5, colour = colores, ID =T, table.col=F, table.legend = F)

invisible(dev.off())

#Zscore P value MAGMA
circ <- circ0
circ$adj_pval <- circ$P.MAGMA
jpeg(file = 'zscore_vs_pMAGMA2.jpeg', units = 'in', width = 12, height = 10, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ, labels = 2.3, colour = colores, ID =T, table.col=F, table.legend = F)
invisible(dev.off())

