# download and clean annotation info for dataset GSE82191

library(io)
library(GEOquery)
Sys.setenv(VROOM_CONNECTION_SIZE = "500000");

gse_data <- getGEO("GSE82191")[[1]];
pheno_df <- pData(phenoData(gse_data));

pheno_df <- pheno_df[, c("external id:ch1", "geo_accession", "source_name_ch1", "age:ch1", 
  "cancersite:ch1", "clinicaldiagnosis:ch1", "diseasestaging:ch1", 
  "Sex:ch1", "tissue:ch1", "tumorgrading:ch1")];

# Clean

pheno_df$source_name_ch1 <- tolower(gsub(" ", "_", pheno_df$source_name_ch1));
pheno_df$`age:ch1` <- tolower(pheno_df$`age:ch1`);
pheno_df$`cancersite:ch1` <- tolower(gsub(" ", "_", pheno_df$`cancersite:ch1`));
pheno_df$`clinicaldiagnosis:ch1` <- tolower(gsub(" ", "_", pheno_df$`clinicaldiagnosis:ch1`));
pheno_df$`diseasestaging:ch1` <- tolower(gsub(" ", "_", pheno_df$`diseasestaging:ch1`));
pheno_df$`Sex:ch1` <- tolower(gsub(" ", "_", pheno_df$`Sex:ch1`));
pheno_df$`tissue:ch1` <- tolower(gsub(" ", "_", pheno_df$`tissue:ch1`));
pheno_df$`tumorgrading:ch1` <- tolower(gsub(" ", "_", pheno_df$`tumorgrading:ch1`));

pheno_df[pheno_df == "[Not Available]"] = NA;

table(pheno_df$source_name_ch1)
pheno_df[1:2,]

colnames(pheno_df) <- c("bcr_patient_barcode", "geo_accession", "cancer_type",
  "age", "cancer_site", "clinical_diagnosis", "disease_stage",
  "sex", "tissue", "tumor_grade");

# use the TCGA "trimmed" name as the rownames
tcga_ids <- pheno_df$bcr_patient_barcode;
tcga_ids <- substr(tcga_ids, 1, 15);
tcga_ids[duplicated(tcga_ids)]; # none

rownames(pheno_df) <- tcga_ids;

pheno_df[1:5, ]

# convert to all character for later merging
pheno_df[] <- lapply(pheno_df, as.character);

write.table(pheno_df, "gse82191_pheno.tsv", row.names = TRUE, 
  sep = "\t", quote = FALSE);