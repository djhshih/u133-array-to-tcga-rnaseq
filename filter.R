# Map microarray expression to RNA-seq expression
# microrray: Affymetrix U133 (fRMA normalized), OVCA
# RNA-seq: TCGA-OVCA from PanCanAtlas
# TODO Expand to other cancer types

library(io)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# for production use
#library(array2rnaseq)

# for development only
library(devtools)
load_all("~/projects/r/array2rnaseq")


db <- org.Hs.eg.db
marr <- qread("data/micro_exprs_matched.rds")
rseq <- qread("data/rna_seq_exprs_t_matched.rds")
dim(marr) # 12664   294
dim(rseq) # 12664   294

# sanity check for matched rows and columns
stopifnot(rownames(marr) == rownames(rseq))
stopifnot(colnames(marr) == colnames(rseq))

# mapping gene names
entrezs <- rownames(marr)
genes <- mapIds(db, keys = entrezs, keytype = "ENTREZID", column = "SYMBOL")
names(genes) <- genes
# check that there are no duplicated gene symbols
stopifnot(sum(duplicated(genes)) == 0)

# TODO move to array2rnaseq
rownames(marr) <- genes
rownames(rseq) <- genes
# Delete genes expressed 0 among all samples
uniq.num <- apply(rseq, 1, function(x) length(unique(x)))
uniq.lst <- which(uniq.num == 1)
marr <- marr[-uniq.lst, ]
rseq <- rseq[-uniq.lst, ]
genes <- genes[-uniq.lst]
entrezs <- entrezs[-uniq.lst]
dim(marr) # 12659  294
dim(rseq) # 12659   294


# Remove these genes which have high NA: pipeline deficiency?
rseq.p.na <- apply(rseq, 1, function(z) mean(is.na(z)))
hist(rseq.p.na, breaks = 100)
na.idx <- rseq.p.na <= 0.5
rseqt <- rseq[na.idx, ]
marrt <- marr[na.idx, ]
genes.f <- genes[na.idx]
entrezs.f <- entrezs[na.idx]
dim(marrt) # 12394   294
dim(rseqt) # 12394   294



marrt <- as.matrix(marrt)
rseqt <- as.matrix(rseqt)
smoothScatter(marrt, rseqt)
plot(marrt, rseqt, pch = ".")
# Mean expressed for each probes
rseqt.m <- rowMeans(rseqt)
marrt.m <- rowMeans(marrt)
plot(marrt.m, rseqt.m, pch = ".")
# Range expressed for each probes
marrt.lim <- c(min(marrt), max(marrt))
rseq.lim <- c(min(rseqt, na.rm = TRUE), max(rseqt, na.rm = TRUE))


# Calculate Mutual information for each probes
mutual_I <- mutual_info(marrt, rseqt)
hist(mutual_I, breaks = 100)


probes.info <- data.frame(
  genes = genes.f,
  entrezs = entrezs.f,
  marr_mean = marrt.m,
  rseq_mean = rseqt.m,
  mutual = mutual_I
)

# TODO output table of functional and non-functional probes

# Non-functional probes
# CONDITION:  mutual_I < 0.05
MI.cod <- quantile(mutual_I, 0.15)
nof.cond <- mutual_I < MI.cod
probes.nof <- subset(probes.info, nof.cond)
probes.nof <- subset(probes.info, nof.cond)
rseqt.nof <- rseqt[probes.nof$genes, ]
marrt.nof <- marrt[probes.nof$genes, ]
dim(probes.nof)[1] # 1859

# Functional probes
probes.f <- subset(probes.info, !nof.cond)
dim(probes.f)[1] # 10535

# write out probe annotation for all probes
qwrite(probes.f, "annot/probes.rds");

