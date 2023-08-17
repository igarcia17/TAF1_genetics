# TAF1_genetics
Genetics analysis of TAF1 and their interactors in relation to multiple sclerosis (MS).
1. [DONE] GSEA RNAseq data
   From two DeSeq2 results, a) Ribominus and b) PolyA, both aligned using Salmon.
   Datasets: GO, HPO, Hallmark, KEGG, REACTOME. Used the identifiers on Mus musculus, NOT the equivalents in human, because of duplicates and lack of correspondance. 
   Very interesting and straightforward results for KEGG of polyA.
   In comparison of both PolyA and Ribominus, only selected categories with p value adjusted < 0.01. 55 in common.
   
2. Gene-set analysis of MAGMA for in silico validation of GSEA results. Check for MS, esophageal-gastric cancer and ALS.
3. GO categories of interactors in column O, sheet DATA of file 'interactores_and_other_data.xlsx' and gene-set analysis of emerged categories.
4. [DONE] Gene-set analysis of MAGMA for those interactors, 4 sets: only WT, only mut, not differential, all together. Found in column O.
5. ORA with oligodendrocytes markers over the genes on the CHIP-seq pause peak.
