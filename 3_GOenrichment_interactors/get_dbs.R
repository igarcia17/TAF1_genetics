library('dplyr')
library('msigdbr')

inputF <- 'not_diff_not_sig_GOid.txt'
resF <- 'not_diff_not_sig_final.txt'

#Enriched categories to validate
filter <- as.data.frame(read.csv(inputF, sep = '\t', header = F)[,1])
colnames(filter) <- 'gs_exact_source'

#GO and HPO
GO_set <- as.data.frame(msigdbr(species = 'Homo sapiens', category = 'C5', 
                   subcategory = NULL)%>% 
  dplyr::select(gs_exact_source, ensembl_gene))

#Bind all three data bases
target_sets <- merge(GO_set, filter, by='gs_exact_source')

#Check merge
uni <- unique(target_sets[c('gs_exact_source')])
filt <- unique(filter[c('gs_exact_source')])

#Save
write.table(target_sets, file = resF, sep="\t",  col.names=F, quote = F, row.names = F)
