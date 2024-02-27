#!/bin/bash
# code modified from:
# https://github.com/djhshih/analysis-xpa-geo/blob/master/GSE70818/get.sh
# download raw .CEL files from GEO, accession number is gse68793

set -euo pipefail
IFS=$'\n\t'

accession=GSE68793
outdir_name="./raw/"

url=ftp://ftp.ncbi.nlm.nih.gov/geo/series
outdir=raw

mkdir -p $outdir
cd $outdir

# download and extract main raw archive file
group=${accession%???}nnn
file=${accession}_RAW.tar
curl -o ${file} ${url}/${group}/${accession}/suppl/$file
tar -xf ${file} && rm ${file}
gunzip *.gz

# remove all the .txt files
rm -rf *.txt

# create tsv table with all directory
ls -I "*.tsv" | awk -v path="$outdir_name" '{print path $0}' > ../raw_lst.tsv

cd -

sed -i '1i files' raw_lst.tsv

# get pheno data
Rscript ./download_pheno.R

# Run frma
Rscript run_frma.R

