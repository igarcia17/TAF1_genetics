# TAF1_genetics
Genetics analysis of TAF1 and their interactors in relation to multiple sclerosis (MS).
1. [DONE] GSEA RNAseq data

   From two DeSeq2 results, a) Ribominus and b) PolyA, both aligned using Salmon.

   Datasets: GO, HPO, Hallmark, KEGG, REACTOME. Used the identifiers on Mus musculus, NOT the equivalents in human, because of duplicates and lack of correspondance. 

   Very interesting and straightforward results for KEGG of polyA.

   In comparison of both PolyA and Ribominus, only selected categories with p value adjusted < 0.01. 55 in common.

   Summary of all results, including those on analyses 2, in this folder.

2. [DONE] Gene-set analysis of MAGMA for in silico validation of GSEA results.

   Done for both ribominus and polyA, checking MS, esophageal cancer and ALS.

3. [DONE] GO categories of interactors in column O, sheet DATA of file 'interactores_and_other_data.xlsx' and gene-set analysis of emerged categories.

   GO categories for the groups 'higher in WT', 'higher in mut' and 'not differential'. Discarded no interactors.

   GO enrichment analysis (powered by PANTHER) at http://geneontology.org/ 

   To get the genes in each GO term, I used the R package msigdb. There were some categories that were not included in this software, therefore this results are PROVISIONAL.
   
   Gene set in MAGMA in multiple sclerosis of the unique categories.
   
4. [DONE] Gene-set analysis of MAGMA for interactors in proteomics studies. 
   - Against the Multiple Sclerosis (MS) gene-based analysis (raw file) constructed over the list of 60k genes (hg38 liftover)
   - Against the MS gene-based analysis (raw file) constructed over the list of 19k genes (in hg19, mainly protein coding)

   There were two lists of interactors. i) Original/list1, 4 sets: only WT, only mut, not differential, all together; ii) New list of interactors based on oligodendrocites, 2 sets: WT>D38 and D38>WT.

5. [DONE] ORA with oligodendrocytes markers over the genes on the CHIP-seq pause peak. Actually, it is an hypergeometric test. All significant.
