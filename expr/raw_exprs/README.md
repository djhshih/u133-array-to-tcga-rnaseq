# Bulk expression data

We are using the The Cancer Genome Atlas (TCGA) dataset.
- the RNA-seq dataset is downloaded from the PanCanAtlas
- the microarray dataset is downloaded from the GEO

If you have installed the `gdc-client`, simply run `get.sh` to download both the RNA-seq and microarray datasets. If you want to download a specific dataset, please refer to the sections below.

## Download and preprocess RNA-seq bulk expression data
Ensure `gdc-client` is installed, otherwise, please refer to [gdc-client github](https://github.com/NCI-GDC/gdc-client) and [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) for instruction on installation.

To download and preprocess the raw TCGA dataset from PanCanAtlas:

- script: `./rna-seq/get_pancan.sh`
- outputs: 
    - `./rna-seq/pancan_exprs.tsv`
    - `./rna-seq/patient_barcode_ref.tsv`

## Download and preprocess u133 microarray bulk expression data
Run the `get.sh` script in each directory in `./u133/`.
