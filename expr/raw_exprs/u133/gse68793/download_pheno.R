# download and clean annotation info for dataset GSE68793

library(io)
library(GEOquery)
Sys.setenv(VROOM_CONNECTION_SIZE = "500000");

gse_data <- getGEO("GSE68793")[[1]];
pheno_df <- pData(phenoData(gse_data));

pheno_df <- pheno_df[, c("title", "geo_accession", "source_name_ch1", "age:ch1", 
  "clinical m staging:ch1", "disease staging:ch1", "histologic type:ch1", 
  "number pack years smoked:ch1", "pathologic n staging:ch1", 
  "pathologic t staging:ch1",  "Sex:ch1", "smoking history:ch1")];

# Clean
pheno_df$title <- sub("EXP-592: ", "", pheno_df$title)

pheno_df$source_name_ch1 <- tolower(gsub(" ", "_", pheno_df$source_name_ch1));

pheno_df$`disease staging:ch1` <- gsub("Stage ", "", pheno_df$`disease staging:ch1`);
pheno_df$`histologic type:ch1` <- tolower(gsub(" ", "_", pheno_df$`histologic type:ch1`));
pheno_df$`histologic type:ch1` <- sub("lung_squamous_cell_carcinoma-_not_otherwise_specified_\\(nos\\)",
  "lung_squamous_cell_carcinoma", pheno_df$`histologic type:ch1`)

pheno_df$`Sex:ch1` <- tolower(pheno_df$`Sex:ch1`);

pheno_df$`smoking history:ch1` <- tolower(gsub(" ", "_", pheno_df$`smoking history:ch1`));

pheno_df[pheno_df == "[Not Available]"] = NA;

table(pheno_df$source_name_ch1)
pheno_df[1:2,]

colnames(pheno_df) <- c("bcr_patient_barcode", "geo_accession", "cancer_type",
  "age", "clinical_m_stage", "disease_stage", "histologic_type", "cigarette_pack_per_year",
  "pathologic_n_stage", "pathologic_t_stage", "sex", "smoking_history")

# use the TCGA "trimmed" name as the rownames
tcga_ids <- pheno_df$bcr_patient_barcode;
tcga_ids <- substr(tcga_ids, 1, 15);
tcga_ids[duplicated(tcga_ids)]; # TCGA-21-1076-01

pheno_df$bcr_patient_barcode[grepl("TCGA-21-1076-01", tcga_ids)];
# "TCGA-21-1076-01A-01R-0690-01" "TCGA-21-1076-01A-02R-0690-01"
# different portion, pick 01 because all the others are from 01 portion

remove_idx <- which(pheno_df$bcr_patient_barcode == "TCGA-21-1076-01A-02R-0690-01");
pheno_df <- pheno_df[-remove_idx, ]

rownames(pheno_df) <- substr(pheno_df$bcr_patient_barcode, 1, 15);

pheno_df[1:5, ]

# convert to all character for later merging
pheno_df[] <- lapply(pheno_df, as.character);

write.table(pheno_df, "gse68793_pheno.tsv", row.names = TRUE, 
  sep = "\t", quote = FALSE);