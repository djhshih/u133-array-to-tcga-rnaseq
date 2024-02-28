# match the RNA-seq with the microarray

library(io)
library(hthgu133a.db)

source("../../../R/clean.R")

rna_exprs <- qread("../../raw_exprs/rna-seq/rna_exprs.rds");
microarray_exprs <- qread("../../raw_exprs/u133/gse68793/gse68793_exprs.rds");

dim(rna_exprs) # 20531 11069
rna_exprs[1:3, 1:4]

dim(microarray_exprs) # 22277   134
microarray_exprs[1:3, 1:4]

# match samples
# first change the format of the sample names
colnames(microarray_exprs) <- substr(colnames(microarray_exprs), 1, 12);

intersect_samples <- intersect(colnames(rna_exprs), colnames(microarray_exprs))

microarray_exprs <- microarray_exprs[, match(intersect_samples, colnames(microarray_exprs))]
rna_exprs <- rna_exprs[, match(intersect_samples, colnames(rna_exprs))]

dim(microarray_exprs) # 22277   131
dim(rna_exprs) # 20531   131


# match the gene names
# convert probe ids to entrez ids
microarray_exprs_m <- as.matrix(microarray_exprs)
microarray_genes <- rownames(microarray_exprs_m);
map <- AnnotationDbi::select(hthgu133a.db, keys = microarray_genes, 
  keytype = "PROBEID", columns = "ENTREZID");
rownames(microarray_exprs_m) <- map$ENTREZID[match(microarray_genes, map$PROBEID)];
dim(microarray_exprs_m) # 22277   131

microarray_exprs_m <- collapse_duplicate_probe(microarray_exprs_m);
dim(microarray_exprs_m) # 13039   131

microarray_exprs <- data.frame(microarray_exprs_m)
dim(microarray_exprs)

# converted to matrix will change the colnames,
# so change back to the correct format
colnames(microarray_exprs) <- gsub("\\.", "-", colnames(microarray_exprs))

intersect_genes <- intersect(rownames(rna_exprs), rownames(microarray_exprs))

microarray_exprs <- microarray_exprs[match(intersect_genes, rownames(microarray_exprs)), ]
rna_exprs <- rna_exprs[match(intersect_genes, rownames(rna_exprs)), ]

dim(microarray_exprs) # 12728   131
dim(rna_exprs) # 12728   131

stopifnot(colnames(microarray_exprs) == colnames(rna_exprs))
stopifnot(rownames(microarray_exprs) == rownames(rna_exprs))

qwrite(microarray_exprs, "./gse68793_exprs-matched.rds");
qwrite(rna_exprs, "./rna_exprs-matched.rds");
