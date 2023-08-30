library(stats)
seed(13)

oligoden_markersF <- 'general_OL_colM.txt'
active_genesF <- 'top_deficit_elon_FC_SER2P_GENE_menorque_0_7.txt'
output <- paste0('ORA_res_',active_genesF, oligoden_markersF)

experimental_genes <- readLines(active_genesF)
oligodendrocyte_markers <- readLines(oligoden_markersF)

# Set parameters
total_genes <- 20000  # Total number of genes
total_oligodendrocyte_markers <- length(oligodendrocyte_markers)  # Total number of oligodendrocyte markers in your reference list
total_experimental_genes <- length(experimental_genes)  # Total number of genes in your experimental list

# Calculate the p-value using the hypergeometric test
p_value <- 1 - phyper(total_experimental_genes - 1, total_genes, total_oligodendrocyte_markers, total_experimental_genes)

#p-value < 0.05 suggests that the overrepresentation is statistically significant.
#Save
write(p_value, output)
