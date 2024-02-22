# download and clean annotation info for dataset GSE68833

library(io)
library(GEOquery)
Sys.setenv(VROOM_CONNECTION_SIZE = "500000");

gse_data <- getGEO("GSE68833")[[1]];
pheno_df <- pData(phenoData(gse_data));

pheno_df <- pheno_df[, c("title", "geo_accession", "source_name_ch1", 
  "age:ch1", "cytogenetic:ch1", "french-american-british classification:ch1", 
  "organism part:ch1", "Sex:ch1")];

# Clean
pheno_df$title <- sub("EXP-608: Blood_", "", pheno_df$title)

pheno_df$source_name_ch1 <- tolower(gsub(" ", "_", pheno_df$source_name_ch1));

pheno_df$`cytogenetic:ch1` <- tolower(pheno_df$`cytogenetic:ch1`);
pheno_df$`organism part:ch1` <- tolower(gsub("Blood;", "", pheno_df$`organism part:ch1`));
pheno_df$`organism part:ch1` <- gsub(" ", "_", pheno_df$`organism part:ch1`);

pheno_df$`Sex:ch1` <- tolower(pheno_df$`Sex:ch1`);

pheno_df[1:2,]

colnames(pheno_df) <- c("bcr_patient_barcode", "geo_accession", "cancer_type",
  "age", "cytogenetic", "french-american-british_classification", 
  "sample_location", "sex")

rownames(pheno_df) <- pheno_df$bcr_patient_barcode;

pheno_df[1:5, ]

# convert to all character for later merging
pheno_df[] <- lapply(pheno_df, as.character);

write.table(pheno_df, "gse68833_pheno.tsv", row.names = TRUE, 
  sep = "\t", quote = FALSE);