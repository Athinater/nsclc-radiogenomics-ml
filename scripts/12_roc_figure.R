# ==============================================================================
# 12_roc_figure.R
# ROC curves of the primary model of each pipeline (output/primary_models.csv, script 08)
# on the held-out test set, ADC as the positive class. AUC and 95% CI (DeLong, pROC) are
# recomputed from the same prediction files as Table 2.
# Outputs: output/figure3_roc_primary_models.png (600 dpi), output/roc_primary_models_auc.csv
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))
source("scripts/analysis_helpers.R")

pm <- load_primary_models()
curves <- list(); aucs <- list()
for (i in seq_len(nrow(pm))) {
  d  <- read_preds(pm$Preds_file[i])
  r  <- roc_adc(d$true_label, d$prob_SCC)
  a  <- auc_ci_adc(d$true_label, d$prob_SCC)
  co <- pROC::coords(r, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  co <- co[order(1 - co$specificity, co$sensitivity), ]
  curves[[i]] <- data.frame(Model = pm$Label[i], fpr = 1 - co$specificity, tpr = co$sensitivity)
  aucs[[i]] <- data.frame(Model = pm$Label[i], Pipeline = pm$Pipeline[i], AUC = round(a$AUC, 3),
                          CI_low = round(a$CI_Low, 3), CI_high = round(a$CI_High, 3),
                          Table_AUC = round(pm$Test_AUC[i], 3),
                          n_ADC = sum(d$true_label == 0), n_SCC = sum(d$true_label == 1))
}
curve_df <- do.call(rbind, curves)
auc_df   <- do.call(rbind, aucs)
if (any(abs(auc_df$AUC - auc_df$Table_AUC) > 0.0005)) warning("Recomputed AUC differs from primary_models.csv.")
print(auc_df, row.names = FALSE)
write.csv(auc_df, out_path("roc_primary_models_auc.csv"), row.names = FALSE)

models <- pm$Label
curve_df$Model <- factor(curve_df$Model, levels = models)
legend_labels <- setNames(ifelse(auc_df$AUC >= 1,
                                 sprintf("%s: AUC %.3f (95%% CI not estimable)", auc_df$Model, auc_df$AUC),
                                 sprintf("%s: AUC %.3f (95%% CI %.3f to %.3f)",
                                         auc_df$Model, auc_df$AUC, auc_df$CI_low, auc_df$CI_high)), auc_df$Model)
p <- ggplot(curve_df, aes(x = fpr, y = tpr, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_path(linewidth = 0.9) +     # empirical ROC (ties give diagonal segments)
  scale_colour_manual(values = setNames(pipeline_colours[pm$Pipeline], models), breaks = models,
                      labels = legend_labels[models]) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = 9), legend.key.width = unit(1.2, "cm"))
fig_file <- out_path("figure3_roc_primary_models.png")
ggsave(fig_file, p, width = 13, height = 15, units = "cm", dpi = 600, bg = "white")
cat("Saved:", fig_file, "\n")
cat("Sensitivity = proportion of ADC correctly identified; specificity = proportion of SCC correctly identified.\n")
