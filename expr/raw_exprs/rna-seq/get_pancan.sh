#!/bin/bash
# download rna-seq data and annotation

set -euo pipefail

# download exprs
gdc-client download 3586c0da-64d0-4b74-a449-5ff4d9136611

mv ./3586c0da-64d0-4b74-a449-5ff4d9136611/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv ./pancan_exprs.tsv

rm -rf 3586c0da-64d0-4b74-a449-5ff4d9136611


# download annotation
#gdc-client download 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81

#mv ./1b5f413e-a8d1-4d10-92eb-7c4ae739ed81/TCGA-CDR-SupplementalTableS1.xlsx ./pancan_pheno.xlsx

# remove the directory
#rm -rf 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81 

# preprocess the data
#Rscript preprocess.R
