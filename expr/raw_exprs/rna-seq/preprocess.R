# preprocess and match the exprs to the phenotype annotation

library(io)
library(openxlsx)

exprs_df <- qread("./pancan_exprs.tsv");
pheno_df <- read.xlsx("./pancan_pheno.xlsx", rowNames = TRUE);

dim(exprs_df) # 20531 11070
dim(pheno_df) # 11160    33

exprs_df[1:5, 1:5]
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

intersect_samples <- intersect(exprs_sample_ids, pheno_df$bcr_patient_barcode);
length(intersect_samples) # 10230

intersect_samples[1:5]

# match
pheno_df <- pheno_df[match(intersect_samples, pheno_df$bcr_patient_barcode), ]
dim(pheno_df) # 10230    33

exprs_df <- exprs_df[, match(intersect_samples, colnames(exprs_df))]
dim(exprs_df) # 20531 10230

stopifnot(colnames(exprs_df) == pheno_df$bcr_patient_barcode)

exprs_df[1:5, 1:5]
pheno_df[1:5, 1:5]

qwrite(exprs_df, "exprs_df-matched.rds");
qwrite(pheno_df, "pheno_df-matched.tsv");