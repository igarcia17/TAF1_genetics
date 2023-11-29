setwd("C:/Users/CBM/Desktop/TAF1_genetics/1_GSEA_results")
library(readxl)
library(DESeq2)
library(GOplot)

#Initial color vector and input files
colores <- c( 'grey','red')
categoriesF <- "Salmon_polyA_independent_filteredClda/GSEA_results_withMAGMA_filteredClda_independent.xlsx"
input_DESEQ <- '../original_files/DeSeq_results/Salmon_polyA_independent_filteredClda.txt'

#Function to select categories
remove_cats <- function(df, dbs){
  result_df <- df[!(df$category %in% dbs), ]
  return(result_df)
}


#load DeSeq data
res <- read.csv(input_DESEQ, sep = "\t")
rownames(res) <- res$GENE_ID
res$FoldChange <- 2^res$log2FoldChange  #have actual fold change

allgenes <- as.data.frame(res$log2FoldChange)
rownames(allgenes) <- rownames(res)
colnames(allgenes) <- 'logFC'
allgenes$gene <- rownames(res)

#load GSEA categories
circ0 <- read_excel(categoriesF, sheet = "bubble_input")

#To compute the Z-score
cats <- as.data.frame(circ0$term)
cats$gene <- circ0$gene
cats$set_size <- circ0$count
colnames(cats) <- c('term', 'gene', 'set_size')

vec <- runif(57, min=0.000001, max= 0.000006)

# Assuming you have loaded your data frames 'cats' and 'allgenes'

# Create an empty vector to store Z-scores
z_scores <- numeric(nrow(cats))

# Iterate through each row of 'cats' and calculate the Z-score
for (i in 1:nrow(cats)) {
  # Extract the list of genes in the current term and split it into a vector
  genes_in_term <- unlist(strsplit(cats$gene[i], "/"))
  
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
circ$adj_pval <- circ$adj_pval + vec

filename <- 'prueba12zscorespotenciados.jpeg'
jpeg(file = filename, units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ, labels = 0, colour = colores, ID =F, table.col=F, table.legend = F)
invisible(dev.off())


###############################Regarding P value MAGMA risk
#Z-score P value MAGMA risk
circ1 <- circ
#circ1$zscore <- circ1$adj_pval
circ1$adj_pval <- circ$MAGMA_P_MS_risk
circ1 <- circ1[circ1$adj_pval < 0.05, ]

filename <- 'prueba8.jpeg'
jpeg(file = filename, units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ1, labels = 4, colour = colores, ID =T, table.col=F, table.legend = F)
invisible(dev.off())

#Removing unwanted cats
dbs <- c('REACTOME', 'KEGG', 'HALLMARK', 'GOMF', 'GOCC', 'GOBP')
circ1_1 <- remove_cats(circ1, dbs=NULL)

filename <- 'zscore_vs_pMAGMA_risk_GO_HPO.jpeg'
jpeg(file = filename, units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ1_1, labels = 4, colour = colores, ID =T, table.col=F, table.legend = F)
invisible(dev.off())


###############################Regarding P value MAGMA severity
#Zscore P value MAGMA
circ2 <- circ
#circ2$zscore <- circ2$adj_pval
circ2$adj_pval <- circ$MAGMA_P_MS_severity
circ2 <- circ2[circ2$adj_pval < 0.05, ]

filename <- 'zscore_vs_pMAGMA_severity_allcats.jpeg'
jpeg(file = filename, units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ2, labels = 0.5, colour = colores, ID =T, table.col=F, table.legend = F)
invisible(dev.off())

dbs <- c('REACTOME', 'KEGG', 'HALLMARK')
circ2_1 <- remove_cats(circ2, dbs=dbs)

filename <- 'zscore_vs_pMAGMA_severity_GO_HPO.jpeg'
jpeg(file = filename, units = 'in', width = 15, height = 20, res = 300)
par(mar = c(2, 2, 2, 5)) 
GOBubble(circ2_1, labels = 1.23, colour = colores, ID =T, table.col=F, table.legend = F)
invisible(dev.off())
