# run frma
# BiocManager::install("frma")
# BiocManager::install("affy")
# BiocManager::install("hthgu133afrmavecs")

library(io)
library(frma)
library(affy)
library(hthgu133afrmavecs)

setwd("./raw")
affy_batch_obj <- ReadAffy();
data_frma <- frma(affy_batch_obj);
data_frma_exprs <- exprs(data_frma);
dim(data_frma_exprs); # 22277   228

qdraw(GNUSE(data_frma, type = "plot", pch = ".", las = 2, cex.axis = 0.1),
  "../qc_boxplots.pdf")
exprs_qc <- GNUSE(data_frma, type = "stat");
qwrite(exprs_qc, "../qc.tsv")

qdraw(hist(data_frma_exprs, breaks = 100), "../exprs_hist.pdf")

# keeping only the GSM ids for the colnames
colname_lst <- colnames(data_frma_exprs);
colname_lst <- substr(colname_lst, 1, 10);

colnames(data_frma_exprs) <- colname_lst;

setwd("../");

# load pheno and keep exprs that can match to pheno
pheno_df <- qread("gse68850_pheno.tsv");
pheno_df[1:2,]

intersect_samples <- intersect(pheno_df$geo_accession, colnames(data_frma_exprs));

pheno_df <- pheno_df[match(intersect_samples, pheno_df$geo_accession), ]
data_frma_exprs <- data_frma_exprs[, match(intersect_samples, colnames(data_frma_exprs))];

stopifnot(colnames(data_frma_exprs) == pheno_df$geo_accession);

colnames(data_frma_exprs) <- rownames(pheno_df);
dim(data_frma_exprs) # 22277   213
data_frma_exprs[1:3, 1:6]

qwrite(data_frma_exprs, "gse68850_exprs.rds");