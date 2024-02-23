# download and clean annotation info for dataset GSE68850

library(io)
library(GEOquery)
Sys.setenv(VROOM_CONNECTION_SIZE = "500000");

gse_data <- getGEO("GSE68850")[[1]];
pheno_df <- pData(phenoData(gse_data));

pheno_df <- pheno_df[, c("title", "geo_accession", "source_name_ch1", 
  "clinical diagnosis:ch1", "disease state:ch1", "histologic diagnosis:ch1", 
  "organism part:ch1", "pathologic status:ch1", "tissue anatomic site:ch1")];

# Clean
pheno_df$title <- sub("gabri-00003: ", "", pheno_df$title)

pheno_df$source_name_ch1 <- tolower(gsub(" ", "_", pheno_df$source_name_ch1));
pheno_df$`disease state:ch1` <- tolower(gsub(" ", "_", pheno_df$`disease state:ch1`));
pheno_df$`histologic diagnosis:ch1` <- tolower(gsub(" ", "_", pheno_df$`histologic diagnosis:ch1`));
pheno_df$`organism part:ch1` <- tolower(pheno_df$`organism part:ch1`);
pheno_df$`pathologic status:ch1` <- tolower(gsub(" ", "_", pheno_df$`pathologic status:ch1`));
pheno_df$`tissue anatomic site:ch1` <- tolower(gsub(" ", "_", pheno_df$`tissue anatomic site:ch1`));

pheno_df[pheno_df == "[Not Available]"] = NA;

pheno_df[1:2,]

colnames(pheno_df) <- c("bcr_patient_barcode", "geo_accession", "cancer_type",
  "clinical_diagnosis", "disease", "histologic_diagnosis", "tissue_location",
  "pathologic_status", "tissue_anatomic_site")

# use the TCGA "trimmed" name as the rownames
tcga_ids <- pheno_df$bcr_patient_barcode;
tcga_ids <- substr(tcga_ids, 1, 15);
tcga_ids[duplicated(tcga_ids)];

#[1] "TCGA-06-0168-01" "TCGA-06-0145-01" "TCGA-06-0145-01" "TCGA-06-0145-01" "TCGA-06-0138-01"
#[6] "TCGA-06-0148-01" "TCGA-06-0148-01" "TCGA-06-0148-01" "TCGA-06-0154-01" "TCGA-06-0156-01"
#[11] "TCGA-06-0156-01" "TCGA-06-0137-01" "TCGA-06-0211-01" "TCGA-06-0176-01" "TCGA-06-0208-01"

# same samples but different portions, pick 01 because all the others are from 01 portion
# if it is missing the 01 portion, randomly pick a portion

dup_barcodes <- pheno_df$bcr_patient_barcode[tcga_ids %in% tcga_ids[duplicated(tcga_ids)]]
dup_barcodes_short <- substr(dup_barcodes, 1, 19);

# TCGA-06-0154, TCGA-06-0176 and missing 01 portion
# TODO?: look at the expression and pick the one with higher expression values

# both TCGA-06-0208 and TCGA-06-0211 samples have 01 portion, but from different vial
# pick vial A because all others are from vial A
keep_barcodes <- c(
  "TCGA-06-0168-01A-01", "TCGA-06-0145-01A-01", "TCGA-06-0138-01A-01", 
  "TCGA-06-0148-01A-01", "TCGA-06-0154-01A-02", "TCGA-06-0156-01A-01",
  "TCGA-06-0137-01A-01", "TCGA-06-0211-01A-01", "TCGA-06-0176-01A-02",
  "TCGA-06-0208-01A-01"
)

remove_barcodes <- dup_barcodes[!dup_barcodes_short %in% keep_barcodes]

remove_idx <- which(pheno_df$bcr_patient_barcode %in% remove_barcodes);
pheno_df <- pheno_df[-remove_idx, ]

rownames(pheno_df) <- substr(pheno_df$bcr_patient_barcode, 1, 15);

pheno_df[1:5, ]

# convert to all character for later merging
pheno_df[] <- lapply(pheno_df, as.character);

write.table(pheno_df, "gse68850_pheno.tsv", row.names = TRUE, 
  sep = "\t", quote = FALSE);