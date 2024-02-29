#!/bin/bash
# download all rna-seq and u133 data

set -euo pipefail
IFS=$'\n\t'

# rna-seq
cd ./rna-seq/
bash get_pancan.sh
cd ..

# u133
cd ./u133/

cd gse68793
bash get.sh
cd ..

cd gse68833
bash get.sh
cd ..

cd gse68850
bash get.sh
cd ..

cd gse82191
bash get.sh
cd ..

cd ..
