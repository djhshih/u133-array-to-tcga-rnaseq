#!/bin/bash
# download rna-seq data and annotation

set -euo pipefail

# download exprs
gdc-client download 3586c0da-64d0-4b74-a449-5ff4d9136611

mv ./3586c0da-64d0-4b74-a449-5ff4d9136611/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv ./pancan_exprs.tsv

rm -rf 3586c0da-64d0-4b74-a449-5ff4d9136611
