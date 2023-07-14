#!/bin/bash

if [ $# -eq 0 ]; then
echo "Error : No argument provided."
exit 1
fi

dataset=$1
size=$2
if [ -z "$3" ]; then
  # Set the default option
annot_file="/home/proyectos/genpsych/PGC_datasets/glist_hg19_from_hg38_v2.txt"
else
  # Use the provided argument
annot_file=$3
fi

snp_loc="./snps_loc/${dataset}_snps.loc"
out_annot="${dataset}_5kb_annot"
annot="./${out_annot}.genes.annot"
pvalfile="./pval_files/${dataset}_pval"
out_genebased="${dataset}_genebased"

#Load modules
module load plink/1.07
module load magma/1.09

#Annotation step
magma --annotate window=5 --snp-loc $snp_loc --gene-loc $annot_file --out $out_annot
echo 'Finished annotation'

#Gene based analysis
magma --bfile ./g1000_eur --gene-annot $annot --pval $pvalfile N=$size --gene-model multi=snp-wise --out $out_genebased
echo 'Finished gene based'


