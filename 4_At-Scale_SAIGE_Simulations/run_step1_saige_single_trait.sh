#!/bin/bash
#Extract command line parameters
plink_file=$1
pheno_file=$2
phenotype=$3
out_dir=$4

#Run step 1
Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --plinkFile=$plink_file --phenoFile=$pheno_file --invNormalize=FALSE --phenoCol=$phenotype --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=binary --outputPrefix=./$phenotype.step1_output --minMAFforGRM 0.01 --LOCO T  --IsOverwriteVarianceRatioFile=TRUE --nThreads=1