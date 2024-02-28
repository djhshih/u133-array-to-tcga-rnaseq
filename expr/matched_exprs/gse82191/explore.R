# quickly visualize how are the genes related in rna-seq vs u133

library(io)
library(ggplot2)

rseq <- qread("./rna_exprs-matched.rds");
marr <- qread("./gse82191_exprs-matched.rds")

# sanity check on the colnames and rownames
stopifnot(colnames(marr) == colnames(rseq))
stopifnot(rownames(marr) == rownames(rseq))

hist(log(as.matrix(rseq) + 1), breaks = 100);

# log transform rna-eq
rseq <- log(rseq + 1);

# explore
# per sample
.compare_sample <- function(target_sample, smooth = FALSE){
  # smooth = TRUE will plot the smoothScatter plot
  # otherwise will simply plot the scatter ggplot
  
  rseq_col <- rseq[, target_sample]
  marr_col <- marr[, target_sample]
  
  compare_df <- data.frame("rseq" = rseq_col, "marr" = marr_col);
  if (smooth){
    smoothScatter(
      compare_df$marr, compare_df$rseq, 
      main = target_sample, xlab = "microarray", ylab = "rna-seq"
    )
  }else{
    ggplot(compare_df, aes(x = marr, y = rseq)) + 
      geom_point(size = 0.1) + theme_classic() +
      labs(title = target_sample, x = "microarray", y = "rna-seq")
  }
}

.compare_sample("TCGA-09-1670")
.compare_sample("TCGA-09-1670", smooth = TRUE)

.compare_sample("TCGA-09-1662")
.compare_sample("TCGA-09-1662", smooth = TRUE)


# per gene
.compare_gene <- function(target_gene){
  rseq_row <- rseq[target_gene, ]
  marr_row <- marr[target_gene, ]
  
  compare_df <- data.frame(t(rseq_row), t(marr_row));
  colnames(compare_df) <- c("rseq", "marr")
  ggplot(compare_df, aes(x = marr, y = rseq)) +
    geom_point() + 
    theme_classic() +
    labs(title = target_gene, x = "microarray", y = "rna-seq")
}

# explore
.compare_gene("155060")
.compare_gene("54715")
