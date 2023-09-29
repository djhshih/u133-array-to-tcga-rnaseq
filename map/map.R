# Learn mapping from microarray to RNA-seq

library(io)
library(ggplot2)
library(dplyr)
library(scam)

# for production use
#library(array2rnaseq)

# for development only
library(devtools)
load_all("~/projects/r/array2rnaseq")

# read in data
probes <- qread("../probes/probes.rds");
rseq <- as.matrix(qread("../expr/rna_seq_exprs_t_matched.rds"));
marr <- as.matrix(qread("../expr/micro_exprs_matched.rds"));

# ensure that probes annotation are in correct order
stopifnot(probes$entrez == rownames(rseq))
stopifnot(probes$entrez == rownames(marr))

# annotate gene names
rownames(rseq) <- probes$gene;
rownames(marr) <- probes$gene;

probes.f <- probes[probes$keep, ];
marr.f <- marr[probes$keep, ];
rseq.f <- rseq[probes$keep, ];
J <- nrow(marr.f);

models <- select_models(marr.f, rseq.f);
models[1:5] <- "scam";
table(models)

maps <- array2rnaseq(marr.f, rseq.f, models=models);

preds <- predict.array2rnaseq(maps, marr.f);

# predictions
plot(marr.f, preds$mean, pch=".")

# calibration plot
plot(preds$mean, rseq.f, pch=".", col=models)
abline(a=0, b=1, col="grey")
cor(c(preds$mean), c(rseq.f), use="complete.obs")

# examine specific genes
gene <- 164;
plot(marr.f[gene, ], preds$mean[gene, ], pch=".")
points(marr.f[gene, ], preds$lower[gene, ], pch=".", col="blue")
points(marr.f[gene, ], preds$upper[gene, ], pch=".", col="blue")
points(marr.f[gene, ], rseq.f[gene, ], pch=".", col="red")

rs <- unlist(lapply(1:J,
  function(j) cor(rseq.f[j, ], preds$mean[j, ])
));
summary(rs)

fves <- unlist(lapply(1:J,
  function(j) fve(rseq.f[j, ], preds$mean[j, ])
));
summary(fves)


# TODO refactor and clean up garbage below


# Range expressed for each probes
marr.lim <- c(min(marr), max(marr))
rseq.lim <- c(min(rseq, na.rm = TRUE), max(rseq, na.rm = TRUE))

# Scatter plots of multiple genes fitted on the linear model
scatter.s(X = marrt.f.lin, Y = rseqt.f.lin, probes = probes.f.lin, pred = lin.pred,
               folder = "linear_plots", num.pic = 5, w = NULL,
               level = 0.95, label = FALSE, press = TRUE, view = TRUE
)
# num.pic: num.pic pictures drawn
# label: whether to named the picture(Y/N/S/HS/C)
# view: TRUE/FALSE (whether to view pic)
# press: TRUE/FALSE(whether see pic one by one )



# --- weighted least square ---

wlin.pred <- lapply(1:dim(marrt.f.lin)[1], function(i) {
  linear(i, X = marrt.f.lin, Y = rseqt.f.lin, level = 0.95, w = NULL)
}
)


# Calculate fev for all gene based on SCAM model
fev.lin <- unlist(lapply(1:length(wlin.pred), function(i) {fev.func(rseqt.f.lin[i, ], wlin.pred[[i]]$fit)}))


# Scatter plots of multiple genes fitted on the weighted linear model
scatter.s(X = marrt.f.lin, Y = rseqt.f.lin, probes = probes.f.lin, pred = wlin.pred,
               folder = "weighted_linear_plots", num.pic = 5, w = NULL,
               level = 0.95, label = FALSE, press = TRUE, view = TRUE
)



# --- SCAM model ---

# Nonlinear probes
probes.f.cur <- subset(probes.f, curve.cond)
rseqt.f.cur <- rseqt.f[curve_idx, ]
marrt.f.cur <- marrt.f[curve_idx, ]
dim(marrt.f.cur) # 10298   294


# Implement SCAM model for all genes
k.s <- c(seq(4, 20, 2))
scam.pred <- lapply(1:dim(marrt.f.cur)[1], function(i) {
  SCAM.model(i, k.s, X = marrt.f.cur, Y = rseqt.f.cur, level = 0.95)
  #print(i)
  }
)


# Calculate fev of scam model
fev.scam <- unlist(lapply(1:dim(marrt.f.cur)[1], function(i) 
  {fev.func(rseqt.f.cur[i, ], scam.pred[[i]]$fit)}
  )
)
hist(fev.scam, breaks = 100)



# Scatter plots of multiple genes fitted on SCAM model
scatter.s(X = marrt.f.cur, Y = rseqt.f.cur, probes = probes.f.cur, pred = scam.pred,
          folder = "scam_plots", num.pic = 5, w = NULL,
          level = 0.95, label = FALSE, press = TRUE, view = TRUE
)













# TODO
# 1. check the edf of non-conver probes, see does it equal to 1;


