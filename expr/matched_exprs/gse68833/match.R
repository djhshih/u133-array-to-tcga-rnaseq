# match the RNA-seq with the microarray

library(io)
library(hgu133plus2.db)

source("../../../R/clean.R")

rna_exprs <- qread("../../raw_exprs/rna-seq/rna_exprs.rds");
microarray_exprs <- qread("../../raw_exprs/u133/gse68833/gse68833_exprs.rds");

dim(rna_exprs) # 20531 11069
rna_exprs[1:3, 1:4]

dim(microarray_exprs) # 22277   134
microarray_exprs[1:3, 1:4]

# match samples
intersect_samples <- intersect(colnames(rna_exprs), colnames(microarray_exprs))

microarray_exprs <- microarray_exprs[, match(intersect_samples, colnames(microarray_exprs))]
rna_exprs <- rna_exprs[, match(intersect_samples, colnames(rna_exprs))]

dim(microarray_exprs) # 54675   163
dim(rna_exprs) # 20531   163


# match the gene names
# convert probe ids to entrez ids
microarray_exprs_m <- as.matrix(microarray_exprs)
microarray_genes <- rownames(microarray_exprs_m);
map <- AnnotationDbi::select(hgu133plus2.db, keys = microarray_genes, 
  keytype = "PROBEID", columns = "ENTREZID");
rownames(microarray_exprs_m) <- map$ENTREZID[match(microarray_genes, map$PROBEID)];
dim(microarray_exprs_m) # 54675   163

microarray_exprs_m <- collapse_duplicate_probe(microarray_exprs_m);
dim(microarray_exprs_m) # 21367   163

microarray_exprs <- data.frame(microarray_exprs_m)
dim(microarray_exprs)

# converted to matrix will change the colnames,
# so change back to the correct format
colnames(microarray_exprs) <- gsub("\\.", "-", colnames(microarray_exprs))

intersect_genes <- intersect(rownames(rna_exprs), rownames(microarray_exprs))

microarray_exprs <- microarray_exprs[match(intersect_genes, rownames(microarray_exprs)), ]
rna_exprs <- rna_exprs[match(intersect_genes, rownames(rna_exprs)), ]

dim(microarray_exprs) # 18404   163
dim(rna_exprs) # 18404   163

stopifnot(colnames(microarray_exprs) == colnames(rna_exprs))
stopifnot(rownames(microarray_exprs) == rownames(rna_exprs))

qwrite(microarray_exprs, "./gse68833_exprs-matched.rds");
qwrite(rna_exprs, "./rna_exprs-matched.rds");
