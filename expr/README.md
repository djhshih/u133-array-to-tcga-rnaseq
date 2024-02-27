# Bulk expression data

`./raw_exprs/` contains scripts to download expression from PANCAN (for rna-seq) and from the GEO (for microarray). Expression have been preprocessed (i.e., frma, log, ..etc) and that the column names are the bcr patient barcodes rather than the GEO sample accession ids. Note that phenotype annotation data is also downloaded in order to map to the bcr-patient-barcode.

`/matched_exprs/` contains matched expression; that is, we matched the genes and samples of the rna-seq to the microarray data.
