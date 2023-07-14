#!/bin/bash
module load plink/1.07
module load magma/1.09

#Argument 1 is or 'binary' or 'property', to get what analysis is being done
#Argument 2 stands for the set file or covariates file where the list of genes and sets is provided
#Argument 3 stands for the genebased analysis raw that use as template
#Argument 4 is for the output name

# Check if an argument is provided
if [ -z "$1" ]; then
  echo "No argument provided. Please specify 'binary' or 'property'."
  exit 1
fi

# Check the value of the argument
if [ "$1" = "binary" ]; then
  #Gene-set analysis
magma --gene-results $3 --set-annot $2 col=2,1  --out $4

elif [ "$1" = "property" ]; then
magma --gene-results $3 --gene-covar $2 --out $4
else
  echo "Invalid argument. Please specify 'binary' or 'property'."
fi

