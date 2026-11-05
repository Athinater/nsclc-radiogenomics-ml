# ==============================================================================
# 16_figure3_roc_heatmap.R
# Figure 3 (combined):
#   (a)-(c) test-set ROC curves of the Table 2 models (one curve per scenario; for each
#           scenario the classifier with the highest CV-AUC), one panel per pipeline;
#           the primary model of the pipeline (highest CV-AUC) is drawn thicker and its
#           AUC (95% DeLong CI) is given under the panel title.
#   (d)     heatmap of sensitivity, specificity, MCC and F1 for all 36 pipeline x
#           classifier x scenario combinations (Table S1); the primary models are outlined.
# ADC is the positive class; fixed threshold 0.5. Everything is read from output/
# (scripts 01-08), so the figure, Table 2 and Table S1 use the same numbers.
# Output: output/figure3_combined.png (600 dpi)
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("scripts/analysis_helpers.R")

FS        <- 7          # text size (pt) of axes, legends and titles
FS_TILE   <- 5.5        # text size (pt) of the numbers in the heatmap cells
FIG_WIDTH <- 180        # mm
FIG_HEIGHT <- 180       # mm
pt <- ggplot2::.pt
scenario_colors <- c(S0 = "#2c3e50", SA = "#8e44ad", SB = "#16a085", SC = "#d35400")
heat_colours    <- c("#d73027", "#f7f7f7", "#1a9850")   # colour-blind-safe alternative: c("#b2182b", "#f7f7f7", "#2166ac")

t2 <- read.csv(out_path("table2_best_per_scenario.csv"), stringsAsFactors = FALSE)
s1 <- read.csv(out_path("table_S1_all_models.csv"), stringsAsFactors = FALSE)
scen_used <- intersect(scenarios_all, unique(s1$Scenario))
fmt_ci <- function(auc, lo, hi) if (auc >= 1 - 1e-12) sprintf("AUC %.3f (CI not estimable)", auc) else
  sprintf("AUC %.3f (95%% CI %.3f to %.3f)", auc, lo, hi)

# ---- (a)-(c) ROC panels ---------------------------------------------------------------------
make_roc_panel <- function(p) {
  d <- t2[t2$Pipeline == p, ]
  d <- d[order(match(d$Scenario, scenarios_all)), ]
  curves <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    pr <- read_preds(d$Preds_file[i])
    co <- pROC::coords(roc_adc(pr$true_label, pr$prob_SCC), x = "all",
                       ret = c("specificity", "sensitivity"), transpose = FALSE)
    co <- co[order(1 - co$specificity, co$sensitivity), ]
    data.frame(Scenario = d$Scenario[i], FPR = 1 - co$specificity, TPR = co$sensitivity,
               Primary = isTRUE(as.logical(d$Primary[i])))
  }))
  labels <- setNames(sprintf("%s (%s): AUC %.3f", d$Scenario, d$Classifier, d$Test_AUC), d$Scenario)
  pr_row <- d[as.logical(d$Primary), ]
  box <- sprintf("Primary: %s (%s)\n%s", pr_row$Scenario, pr_row$Classifier,
                 fmt_ci(pr_row$Test_AUC, pr_row$CI_Low, pr_row$CI_High))
  curves$Scenario <- factor(curves$Scenario, levels = d$Scenario)
  ggplot(curves, aes(FPR, TPR, colour = Scenario, group = Scenario)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.25) +
    geom_path(aes(linewidth = Primary)) +
    scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.1), guide = "none") +
    scale_colour_manual(values = scenario_colors[d$Scenario], labels = labels[d$Scenario], name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    coord_equal(clip = "off") +
    labs(title = pipeline_names[[p]], subtitle = box, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal(base_size = FS) +
    theme(plot.title = element_text(face = "bold", size = FS + 1, hjust = 0.5),
          plot.subtitle = element_text(face = "bold", size = FS - 0.5, hjust = 0.5,
                                       colour = scenario_colors[[pr_row$Scenario]]),
          legend.position = "bottom", legend.text = element_text(size = FS - 0.5),
          legend.key.size = unit(2.5, "mm"), legend.spacing.y = unit(0.3, "mm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
          plot.margin = unit(c(2, 3, 2, 3), "pt")) +
    guides(colour = guide_legend(ncol = 1))
}
roc_panels <- lapply(pipelines, make_roc_panel)

# ---- (d) heatmap of the threshold metrics ------------------------------------------------------
short <- c(radiomic = "Rad", radiogenomic = "Rg", genomic = "Gen")
s1$Column <- paste0(short[s1$Pipeline], "\n", s1$Classifier)
col_levels <- as.vector(t(outer(short[pipelines], names(classifiers), paste, sep = "\n")))
metrics <- c("Sensitivity", "Specificity", "MCC", "F1")
hm <- do.call(rbind, lapply(metrics, function(m) {
  v <- s1[[m]]
  data.frame(Metric = m, Column = s1$Column, Scenario = s1$Scenario,
             Value = if (m == "MCC") (v + 1) / 2 else v,          # colour scale 0-1 (MCC rescaled)
             Label = ifelse(is.na(v), "NA", sprintf("%.2f", v)),
             Primary = as.logical(s1$Primary))
}))
hm$Metric   <- factor(hm$Metric, levels = metrics)
hm$Column   <- factor(hm$Column, levels = col_levels)
hm$Scenario <- factor(hm$Scenario, levels = rev(scen_used))
hm$TextColour <- ifelse(!is.na(hm$Value) & (hm$Value > 0.78 | hm$Value < 0.18), "white", "grey20")
n_s <- length(scen_used)

p_heat <- ggplot(hm, aes(Column, Scenario, fill = Value)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = hm[hm$Primary, ], colour = "black", linewidth = 0.7, fill = NA, width = 0.96, height = 0.96) +
  geom_text(aes(label = Label, colour = TextColour), size = FS_TILE / pt) +
  geom_vline(xintercept = c(3.5, 6.5), colour = "grey25", linewidth = 0.5) +
  scale_fill_gradientn(colours = heat_colours, values = c(0, 0.5, 1), limits = c(0, 1),
                       na.value = "grey85", name = "Score") +
  scale_colour_identity() +
  scale_x_discrete(position = "top") +
  facet_wrap(~ Metric, ncol = 2) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = FS) +
  theme(strip.text = element_text(face = "bold", size = FS + 0.5, colour = "grey20"),
        strip.background = element_rect(fill = "grey95", colour = "grey80", linewidth = 0.2),
        axis.text.x = element_text(size = FS - 0.5, face = "bold", colour = "grey20"),
        axis.text.y = element_text(size = FS, face = "bold", colour = "grey30"),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(3, "mm"),
        legend.position = "right", legend.title = element_text(size = FS, face = "bold"),
        legend.text = element_text(size = FS - 0.5),
        legend.key.width = unit(3, "mm"), legend.key.height = unit(12, "mm"),
        plot.margin = unit(c(4, 3, 3, 3), "pt"))

# ---- assemble -------------------------------------------------------------------------------------
combined <- wrap_plots(A = roc_panels[[1]], B = roc_panels[[2]], C = roc_panels[[3]], D = p_heat,
                       design = "ABC\nDDD") +
  plot_layout(heights = c(0.72, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = FS + 2, face = "bold"))

fig_file <- out_path("figure3_combined.png")
ggsave(fig_file, combined, width = FIG_WIDTH, height = FIG_HEIGHT, units = "mm", dpi = 600, bg = "white")
cat("Saved:", fig_file, "\n")
