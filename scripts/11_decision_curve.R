# ==============================================================================
# 11_decision_curve.R
# Decision curve analysis (Vickers & Elkin, 2006) of the primary models on the test set.
# Event = ADC. Net benefit vs "treat all" and "treat none" over threshold probabilities
# 2%-80%.
#   NB(model)      = TP/n - FP/n * pt/(1 - pt)
#   NB(treat all)  = prev - (1 - prev) * pt/(1 - pt)
#   NB(treat none) = 0
# Outputs: output/decision_curve_data.csv, decision_curve_summary.csv,
#          output/figure4_decision_curve.png (600 dpi)
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))
source("scripts/analysis_helpers.R")

dca_pipelines <- c("radiogenomic", "genomic")     # models shown in Figure 4
thresholds    <- seq(0.02, 0.80, by = 0.01)

net_benefit <- function(y, p, pt) {
  n <- length(y)
  vapply(pt, function(t) {
    pred <- p >= t
    sum(pred & y == 1) / n - sum(pred & y == 0) / n * (t / (1 - t))
  }, numeric(1))
}

pm <- load_primary_models()
pm <- pm[pm$Pipeline %in% dca_pipelines, ]
curves <- list(); summary_rows <- list()
for (i in seq_len(nrow(pm))) {
  d    <- read_event_preds(pm$Preds_file[i])
  nb   <- net_benefit(d$y, d$p, thresholds)
  prev <- mean(d$y)
  all_nb <- prev - (1 - prev) * thresholds / (1 - thresholds)
  curves[[i]] <- data.frame(Model = pm$Label[i], Pipeline = pm$Pipeline[i], Threshold = thresholds,
                            NB_model = nb, NB_treat_all = all_nb, NB_treat_none = 0)
  better <- nb > pmax(all_nb, 0)
  summary_rows[[i]] <- data.frame(
    Model = pm$Label[i], Prevalence_of_event = round(prev, 3),
    Thresholds_model_best = sprintf("%d of %d (%.0f%%)", sum(better), length(thresholds), 100 * mean(better)),
    Range_model_best = if (any(better)) sprintf("%.2f-%.2f", min(thresholds[better]), max(thresholds[better])) else "none",
    Thresholds_not_best = if (all(better)) "none" else paste(sprintf("%.2f", thresholds[!better]), collapse = " "),
    Event = event_class)
}
curve_df <- do.call(rbind, curves)
sum_df   <- do.call(rbind, summary_rows)
print(sum_df, row.names = FALSE)
write.csv(curve_df, out_path("decision_curve_data.csv"), row.names = FALSE)
write.csv(sum_df,   out_path("decision_curve_summary.csv"), row.names = FALSE)

# ---- Figure ------------------------------------------------------------------------------
ref <- unique(curve_df[, c("Threshold", "NB_treat_all")])
plot_df <- rbind(
  data.frame(Strategy = curve_df$Model, Threshold = curve_df$Threshold, NB = curve_df$NB_model),
  data.frame(Strategy = "Treat all (ADC)", Threshold = ref$Threshold, NB = ref$NB_treat_all),
  data.frame(Strategy = "Treat none", Threshold = ref$Threshold, NB = 0))
lev <- c(pm$Label, "Treat all (ADC)", "Treat none")
plot_df$Strategy <- factor(plot_df$Strategy, levels = lev)
cols <- c(setNames(pipeline_colours[pm$Pipeline], pm$Label), "Treat all (ADC)" = "grey40", "Treat none" = "black")
ltys <- c(setNames(rep("solid", nrow(pm)), pm$Label), "Treat all (ADC)" = "dashed", "Treat none" = "dotted")
ymax <- max(curve_df$NB_model, ref$NB_treat_all, 0) + 0.05

p <- ggplot(plot_df, aes(Threshold, NB, colour = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols, breaks = lev) +
  scale_linetype_manual(values = ltys, breaks = lev) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.1), limits = range(thresholds)) +
  coord_cartesian(ylim = c(-0.1, ymax)) +
  labs(x = "Threshold probability (ADC)", y = "Net benefit", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1.2, "cm")) +
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))
fig_file <- out_path("figure4_decision_curve.png")
ggsave(fig_file, p, width = 16, height = 12, units = "cm", dpi = 600, bg = "white")
cat("Saved:", out_path("decision_curve_data.csv"), ",", out_path("decision_curve_summary.csv"), "and", fig_file, "\n")
