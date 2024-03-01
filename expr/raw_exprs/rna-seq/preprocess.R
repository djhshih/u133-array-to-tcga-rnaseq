# preprocess and match the exprs to the phenotype annotation

library(io)

exprs_df <- qread("./pancan_exprs.tsv");

dim(exprs_df) # 20531 11070

# keep the entrez ids
gene_ids <- exprs_df$gene_id;
entrez_ids <- unlist(lapply(gene_ids, function(x){
  strsplit(x, "\\|")[[1]][2]
}))

rownames(exprs_df) <- entrez_ids;
exprs_df <- exprs_df[, !colnames(exprs_df) == "gene_id"]

# match samples
exprs_sample_ids <- colnames(exprs_df);
exprs_sample_ids <- substr(exprs_sample_ids, 1, 12); # to match the format

# create a reference table
exprs_sample_ids_df <- data.frame(
  "full_id" = colnames(exprs_df),
  "short_id" = exprs_sample_ids
)
exprs_sample_ids_df[1:5,]
qwrite(exprs_sample_ids_df, "patient_barcode_ref.tsv");

# replace the full column names with the trimmed column names
colnames(exprs_df) <- exprs_sample_ids;

exprs_df[1:5, 1:5]

qwrite(exprs_df, "rna_exprs.rds");