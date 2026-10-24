# ==============================================================================
# 08_delong_comparison.R
# Paired DeLong test: radiogenomic vs genomic primary model (test set, n = 33)
# Reports AUC (95% CI) of each model, AUC difference (95% CI), Z and two-sided p.
# Manuscript values for comparison: Z = 1.11, p = 0.27, AUC difference CI -0.042 to 0.152.
# ==============================================================================
suppressPackageStartupMessages(library(pROC))
source("scripts/analysis_helpers.R")

# Paired DeLong test written out explicitly (DeLong et al., 1988; structural components),
# so that the AUC difference and its confidence interval are available.
# y: 1 = event, 0 = non-event; p1, p2: predicted probabilities of the event from two models
# evaluated on the SAME patients in the SAME order.
delong_paired <- function(y, p1, p2) {
  pos <- which(y == 1); neg <- which(y == 0)
  m <- length(pos); n <- length(neg)
  placement <- function(p) {
    v10 <- vapply(pos, function(i) mean((p[i] > p[neg]) + 0.5 * (p[i] == p[neg])), numeric(1))
    v01 <- vapply(neg, function(j) mean((p[pos] > p[j]) + 0.5 * (p[pos] == p[j])), numeric(1))
    list(auc = mean(v10), v10 = v10, v01 = v01)
  }
  a <- placement(p1); b <- placement(p2)
  s10 <- cov(cbind(a$v10, b$v10)); s01 <- cov(cbind(a$v01, b$v01))
  S   <- s10 / m + s01 / n                         # covariance matrix of (AUC1, AUC2)
  d   <- a$auc - b$auc
  se  <- sqrt(S[1, 1] + S[2, 2] - 2 * S[1, 2])
  z   <- d / se
  list(auc1 = a$auc, auc2 = b$auc, diff = d, se = se,
       ci_low = d - qnorm(0.975) * se, ci_high = d + qnorm(0.975) * se,
       z = z, p = 2 * pnorm(-abs(z)))
}

genomic <- read_event_preds(analysis_models[["Genomic (H2O, SC)"]])
rows <- list()
for (nm in c("Radiogenomic (TF, SB)", "Radiogenomic (TF, SC)")) {
  rg <- read_event_preds(analysis_models[[nm]])
  stopifnot(length(rg$y) == length(genomic$y), all(rg$y == genomic$y))   # same patients, same order
  res <- delong_paired(rg$y, rg$p, genomic$p)

  # cross-check with pROC (AUC direction fixed: event probability higher for events)
  r1 <- roc(rg$y,      rg$p,      levels = c(0, 1), direction = "<", quiet = TRUE)
  r2 <- roc(genomic$y, genomic$p, levels = c(0, 1), direction = "<", quiet = TRUE)
  tst <- roc.test(r1, r2, method = "delong", paired = TRUE)
  ci1 <- as.numeric(ci.auc(r1, method = "delong")); ci2 <- as.numeric(ci.auc(r2, method = "delong"))

  rows[[nm]] <- data.frame(
    Comparison  = paste(nm, "vs Genomic (H2O, SC)"),
    AUC_1       = round(res$auc1, 3), AUC_1_CI = sprintf("%.3f-%.3f", ci1[1], min(ci1[3], 1)),
    AUC_2       = round(res$auc2, 3), AUC_2_CI = sprintf("%.3f-%.3f", ci2[1], min(ci2[3], 1)),
    AUC_diff    = round(res$diff, 3),
    Diff_CI     = sprintf("%.3f to %.3f", res$ci_low, res$ci_high),
    Z           = round(res$z, 2), p_value = round(res$p, 3),
    pROC_Z      = round(as.numeric(tst$statistic), 2), pROC_p = round(tst$p.value, 3),
    n_events    = sum(rg$y == 1), n_nonevents = sum(rg$y == 0)
  )
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, out_path("delong_radiogenomic_vs_genomic.csv"), row.names = FALSE)
cat("Saved: ", out_path("delong_radiogenomic_vs_genomic.csv"), "\n")
cat("Reference values in the manuscript: Z = 1.11, p = 0.27, AUC difference CI -0.042 to 0.152.\n")
