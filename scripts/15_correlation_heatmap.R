# ==============================================================================
# 15_correlation_heatmap.R
# Figure 5: heatmap of the Spearman correlations between the 10 most frequently selected
# radiomic features and the selected genes (output of script 14), with a colour bar for
# the subtype association of each gene. Asterisks mark pairs with BH q < 0.05.
# Output: output/figure5_correlation_heatmap.png (600 dpi)
# ==============================================================================
suppressPackageStartupMessages(library(ggplot2))
source("scripts/common.R")

res  <- read.csv(out_path("radiomic_gene_correlations_all_pairs.csv"), stringsAsFactors = FALSE)
summ <- read.csv(out_path("radiomic_gene_correlations_summary.csv"), stringsAsFactors = FALSE)

# Optional shorter display names for long feature names, e.g.
# label_map <- c("MORPHOLOGICAL_RadiusRoiNorm.MaxIntensityCoor.RoiCentroidCoor.Dist" = "RadiusRoiNorm (MaxIntensity-Centroid)")
label_map <- c()

group_levels  <- c("SCC-associated", "ADC-associated", "Not subtype-associated")
group_colours <- c("SCC-associated" = "#E69F00", "ADC-associated" = "#7B3294", "Not subtype-associated" = "grey65")

# ---- Order: features as in Table 4 (selection frequency), genes by subtype group -----------
feats <- summ$Feature
sel_file <- out_path("selected_genes.csv")
gene_order <- if (file.exists(sel_file)) read.csv(sel_file, stringsAsFactors = FALSE)$Gene else unique(res$Gene)
genes <- unique(res[, c("Gene", "Gene_group")])
genes <- genes[order(match(genes$Gene_group, group_levels), match(genes$Gene, gene_order)), ]

res$x <- match(res$Gene, genes$Gene)
res$y <- length(feats) + 1 - match(res$Feature, feats)          # first feature at the top
res$star <- ifelse(res$Significant, "*", "")
bar <- data.frame(x = seq_len(nrow(genes)), y = length(feats) + 1.1, Gene_group = genes$Gene_group)
bar$Gene_group <- factor(bar$Gene_group, levels = intersect(group_levels, bar$Gene_group))

feat_labels <- ifelse(feats %in% names(label_map), label_map[feats], feats)
# long names are split over two lines at the dot nearest to the middle
wrap_label <- function(x, max_chars = 40) {
  if (nchar(x) <= max_chars) return(x)
  dots <- gregexpr(".", x, fixed = TRUE)[[1]]
  if (dots[1] < 0) return(x)
  cut <- dots[which.min(abs(dots - nchar(x) / 2))]
  paste0(substr(x, 1, cut), "\n", substr(x, cut + 1, nchar(x)))
}
feat_labels <- vapply(feat_labels, wrap_label, character(1))
lim <- ceiling(max(abs(res$rho)) * 10) / 10

p <- ggplot(res, aes(x = x, y = y)) +
  geom_tile(aes(fill = rho), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 4.5, vjust = 0.75) +
  geom_tile(data = bar, aes(x = x, y = y), fill = group_colours[as.character(bar$Gene_group)],
            height = 0.5, inherit.aes = FALSE) +
  geom_point(data = bar, aes(x = x, y = y, colour = Gene_group), shape = 15, size = 0, alpha = 0,
             inherit.aes = FALSE) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "Spearman rho") +
  scale_colour_manual(values = group_colours, name = "Gene association", drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), order = 2),
         fill = guide_colourbar(order = 1)) +
  scale_x_continuous(breaks = seq_len(nrow(genes)), labels = genes$Gene, expand = c(0, 0)) +
  scale_y_continuous(breaks = rev(seq_along(feats)), labels = feat_labels, expand = c(0, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL,
       caption = sprintf("* BH-adjusted q < 0.05 (%d of %d pairs); n = %d patients", sum(res$Significant), nrow(res), 113)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
        axis.text.y = element_text(size = 8.5),
        legend.position = "right",
        plot.caption = element_text(size = 8, hjust = 0))

fig_file <- out_path("figure5_correlation_heatmap.png")
ggsave(fig_file, p, width = 22, height = 12, units = "cm", dpi = 600, bg = "white")
cat("Saved:", fig_file, "\n")
