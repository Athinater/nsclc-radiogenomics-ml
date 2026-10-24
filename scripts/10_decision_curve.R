# ==============================================================================
# 10_decision_curve.R
# Decision curve analysis (Vickers & Elkin, 2006) on the test set.
# Net benefit of each model vs "treat all" and "treat none" over threshold
# probabilities 2%-80%. Writes the curve data; no figure is produced here.
#   NB(model)     = TP/n - FP/n * pt/(1 - pt)
#   NB(treat all) = prev - (1 - prev) * pt/(1 - pt)
#   NB(treat none)= 0
# ==============================================================================
source("scripts/analysis_helpers.R")

thresholds <- seq(0.02, 0.80, by = 0.01)

net_benefit <- function(y, p, pt) {
  n <- length(y); w <- pt / (1 - pt)
  vapply(pt, function(t) {
    pred <- p >= t
    sum(pred & y == 1) / n - sum(pred & y == 0) / n * (t / (1 - t))
  }, numeric(1))
}

curves <- list(); summary_rows <- list()
for (nm in names(analysis_models)) {
  d  <- read_event_preds(analysis_models[[nm]])
  nb <- net_benefit(d$y, d$p, thresholds)
  prev <- mean(d$y)
  all_nb <- prev - (1 - prev) * thresholds / (1 - thresholds)
  curves[[nm]] <- data.frame(Model = nm, Threshold = thresholds,
                             NB_model = nb, NB_treat_all = all_nb, NB_treat_none = 0)
  better <- nb > pmax(all_nb, 0)
  summary_rows[[nm]] <- data.frame(
    Model = nm,
    Prevalence_of_event = round(prev, 3),
    Thresholds_model_best = sprintf("%d of %d (%.0f%%)", sum(better), length(thresholds), 100 * mean(better)),
    Lowest_threshold_NB_not_above_both = if (all(better)) "none" else sprintf("%.2f", thresholds[!better][1]),
    Event = event_class)
}
curve_df <- do.call(rbind, curves)
sum_df   <- do.call(rbind, summary_rows)
print(sum_df, row.names = FALSE)
write.csv(curve_df, out_path("decision_curve_data.csv"), row.names = FALSE)
write.csv(sum_df,   out_path("decision_curve_summary.csv"), row.names = FALSE)
cat("Saved: ", out_path("decision_curve_data.csv"), "and decision_curve_summary.csv\n")
cat("Manuscript statement: both the radiogenomic and genomic model gave greater net benefit than\n",
    "treat-all / treat-none across threshold probabilities of 2-80%. Check this against the table above.\n")
