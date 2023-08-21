library('dplyr')
library('msigdbr')

inputF <- 'GSEAcats_polyA.txt'
resF <- 'validate_sets_polyA.txt'

#Enriched categories to validate
filter <- as.data.frame(read.csv(inputF, sep = '\t', header = F)[,1])
colnames(filter) <- 'gs_name'

#KEGG and REACTOME
C2_set <- as.data.frame(msigdbr(species = 'Homo sapiens', category = 'C2', 
                   subcategory = NULL)%>% 
  dplyr::select(gs_name, ensembl_gene))

#GO and HPO
GO_set <- as.data.frame(msigdbr(species = 'Homo sapiens', category = 'C5', 
                   subcategory = NULL)%>% 
  dplyr::select(gs_name, ensembl_gene))

#Hallmark
H_set <- as.data.frame(msigdbr(species = 'Homo sapiens', category = 'H', 
                                subcategory = NULL)%>% 
                          dplyr::select(gs_name, ensembl_gene))
#Bind all three data bases
sets <- bind_rows(C2_set, GO_set, H_set)
target_sets <- merge(sets, filter, by='gs_name')

#Check merge
uni <- unique(target_sets[c('gs_name')])
dim(filter) == dim(uni)

#Save
write.table(target_sets, file = resF, sep="\t",  col.names=F, quote = F, row.names = F)
