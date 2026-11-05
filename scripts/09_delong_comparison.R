# ==============================================================================
# 09_delong_comparison.R
# Paired DeLong test between the primary models (output/primary_models.csv, script 08),
# on the test patients common to both models, matched by patient ID.
# Primary comparison: radiogenomic vs genomic. The comparisons with the radiomic model
# are reported as secondary.
# Reports AUC (95% CI) of each model, AUC difference (95% CI), Z and two-sided p.
# ==============================================================================
source("scripts/analysis_helpers.R")

# Paired DeLong test written out explicitly (DeLong et al., 1988; structural components),
# so that the AUC difference and its confidence interval are available.
# y: 1 = event, 0 = non-event; p1, p2: predicted probabilities of the event from two models
# for the SAME patients in the SAME order.
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

pm <- load_primary_models()
comparisons <- list(c("radiogenomic", "genomic",  "primary"),
                    c("radiogenomic", "radiomic", "secondary"),
                    c("genomic",      "radiomic", "secondary"))
rows <- list()
for (cmp in comparisons) {
  if (!all(cmp[1:2] %in% pm$Pipeline)) next
  m1 <- pm[pm$Pipeline == cmp[1], ]; m2 <- pm[pm$Pipeline == cmp[2], ]
  a  <- read_event_preds(m1$Preds_file); b <- read_event_preds(m2$Preds_file)
  d  <- pair_by_id(a, b)
  if (nrow(d) < max(nrow(a), nrow(b)))
    warning(m1$Label, " vs ", m2$Label, ": only ", nrow(d), " patients are common to both test sets.")
  res <- delong_paired(d$y, d$p_1, d$p_2)

  # cross-check with pROC (direction fixed: event probability higher for events)
  r1  <- pROC::roc(d$y, d$p_1, levels = c(0, 1), direction = "<", quiet = TRUE)
  r2  <- pROC::roc(d$y, d$p_2, levels = c(0, 1), direction = "<", quiet = TRUE)
  tst <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
  ci1 <- as.numeric(pROC::ci.auc(r1, method = "delong")); ci2 <- as.numeric(pROC::ci.auc(r2, method = "delong"))

  rows[[length(rows) + 1]] <- data.frame(
    Comparison = paste(m1$Label, "vs", m2$Label), Type = cmp[3],
    AUC_1 = round(res$auc1, 3), AUC_1_CI = sprintf("%.3f-%.3f", max(ci1[1], 0), min(ci1[3], 1)),
    AUC_2 = round(res$auc2, 3), AUC_2_CI = sprintf("%.3f-%.3f", max(ci2[1], 0), min(ci2[3], 1)),
    AUC_diff = round(res$diff, 3), Diff_CI = sprintf("%.3f to %.3f", res$ci_low, res$ci_high),
    Z = round(res$z, 2), p_value = signif(res$p, 3),
    pROC_Z = round(as.numeric(tst$statistic), 2), pROC_p = signif(tst$p.value, 3),
    n = nrow(d), n_events = sum(d$y == 1), n_nonevents = sum(d$y == 0), Event = event_class)
}
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, out_path("delong_comparisons.csv"), row.names = FALSE)
cat("Saved:", out_path("delong_comparisons.csv"), "\n")
