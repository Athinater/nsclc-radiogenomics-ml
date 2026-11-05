# ==============================================================================
# 14_radiomic_gene_correlations.R
# Spearman correlations between the 10 most frequently selected radiomic features and
# the PCA-selected genes across the full cohort (10 features x number of genes pairs),
# Benjamini-Hochberg FDR across all pairs (significant at q < 0.05).
# Exploratory: this analysis does not feed back into model training.
# Gene expression: DESeq2 variance-stabilised counts of the full cohort.
# ==============================================================================
source("scripts/common.R")

# ---- Radiomic features ----------------------------------------------------------------
# Default: the 10 radiomic features selected most often in the rerun (final models of the
# feature-selection scenarios, ties broken by CV-fold frequency), from
# output/feature_selection_frequency.csv (script 08). Set use_rerun_features <- FALSE to use
# the list of the submitted manuscript instead.
use_rerun_features <- TRUE
manuscript_radiomic_features <- c(
  "GLSZM_SmallZoneEmphasis", "MedianIntensity", "50thIntensityPercentile", "GLCM_ClusterShade",
  "IBCoefficientOfVariation", "GLSZM_GreyLevelNonUniformity", "GLSZM_ZoneSizeVariance",
  "GLSZM_LargeZoneEmphasis", "NGTDM_Coarseness", "IHMinimumGreyLevel")

# Genes: those selected on the training set in the rerun (output/selected_genes.csv), SCC-associated
# first; without that file, the 18 genes of the submitted manuscript.
sel_file <- out_path("selected_genes.csv")
if (file.exists(sel_file) && "Association" %in% names(read.csv(sel_file, nrows = 1))) {
  sel   <- read.csv(sel_file, stringsAsFactors = FALSE)
  sel   <- sel[order(match(sel$Association, c("SCC-associated", "ADC-associated", "Not subtype-associated"))), ]
  genes <- sel$Gene
  gene_group <- setNames(sel$Association, sel$Gene)
} else {
  warning(sel_file, " not found; using the 18 genes of the manuscript.")
  genes <- manuscript_genes
  gene_group <- setNames(rep(c("SCC-associated", "ADC-associated"), c(10, 8)), manuscript_genes)
}
cat("Genes (", length(genes), "):", paste(genes, collapse = ", "), "\n")

# ---- Data ---------------------------------------------------------------------------
tbl <- load_counts(counts_file)
clinical_data <- read.csv(clinical_rg_file, sep = ";", header = TRUE, row.names = 1)
rownames(clinical_data) <- norm_id(rownames(clinical_data))
coldata <- read.csv(clinical_meta_file, header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1)
rownames(coldata) <- norm_id(rownames(coldata))
coldata <- coldata[coldata$Histology %in% c("Adenocarcinoma", "Squamous cell carcinoma"), , drop = FALSE]
coldata$Histology <- factor(coldata$Histology)

ids <- intersect(intersect(rownames(clinical_data), rownames(coldata)), colnames(tbl))
cat("Patients with radiomic features, clinical metadata and RNA-Seq:", length(ids), "(manuscript: 113)\n")
if (length(ids) != 113) warning("Cohort size differs from the 113 patients reported in the manuscript.")

strip_prefix <- function(x) sub("\\.IBSI\\..*$", "", sub("^(INTENSITY\\.BASED_|INTENSITY\\.HISTOGRAM_|X)", "", x))
freq_file <- out_path("feature_selection_frequency.csv")
if (use_rerun_features && file.exists(freq_file)) {
  fq <- read.csv(freq_file, stringsAsFactors = FALSE)
  fq <- fq[fq$Type == "radiomic", ]
  cols <- head(fq$Feature, 10)              # already ordered by frequency (script 08)
  cat("Radiomic features from the rerun (", freq_file, "):\n  ", paste(cols, collapse = ", "), "\n")
  same <- setequal(strip_prefix(cols), manuscript_radiomic_features)
  cat(if (same) "Same 10 features as in the manuscript.\n" else
      paste0("DIFFERENT from the manuscript list. Not in the manuscript: ",
             paste(setdiff(strip_prefix(cols), manuscript_radiomic_features), collapse = ", "), "\n"))
  key <- intersect(c("Final_models", "CV_folds"), names(fq))
  if (nrow(fq) > 10 && all(unlist(fq[10, key]) == unlist(fq[11, key])))
    warning("Tie at rank 10 of the selection frequency; the alphabetical order decides.")
  missing <- setdiff(cols, colnames(clinical_data))
  if (length(missing) > 0) stop("Features not found in ", clinical_rg_file, ": ", paste(missing, collapse = ", "))
} else {
  if (use_rerun_features) warning(freq_file, " not found; using the manuscript feature list.")
  cols <- vapply(manuscript_radiomic_features, function(f) {
    hit <- colnames(clinical_data)[colnames(clinical_data) == f | strip_prefix(colnames(clinical_data)) == f]
    if (length(hit) != 1) stop("Radiomic feature not found (or ambiguous) in the data: ", f)
    hit
  }, character(1))
}
feature_names <- strip_prefix(cols)
rad <- as.matrix(clinical_data[ids, cols, drop = FALSE]); colnames(rad) <- feature_names

# VST expression of the selected genes (full cohort)
dds <- DESeq2::DESeqDataSetFromMatrix(tbl[, ids], coldata[ids, , drop = FALSE], ~ Histology)
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)   # required for vst(blind = FALSE)
vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
missing_genes <- setdiff(genes, rownames(vsd))
if (length(missing_genes) > 0) stop("Genes not found in the expression data: ", paste(missing_genes, collapse = ", "))
expr <- t(vsd[genes, , drop = FALSE])

# ---- Spearman correlations ------------------------------------------------------------
grid <- expand.grid(Feature = feature_names, Gene = genes, stringsAsFactors = FALSE)
tests <- lapply(seq_len(nrow(grid)), function(i) {
  ct <- suppressWarnings(cor.test(rad[, grid$Feature[i]], expr[, grid$Gene[i]],
                                  method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
})
res <- cbind(grid, do.call(rbind, tests))
res$q_BH <- p.adjust(res$p, method = "BH")
res$Significant <- res$q_BH < 0.05
res$Gene_group <- unname(gene_group[res$Gene])

cat("\nPairs tested:", nrow(res), "| significant at q < 0.05:", sum(res$Significant), "\n")

summary_tbl <- do.call(rbind, lapply(feature_names, function(f) {
  r <- res[res$Feature == f, ]; s <- r[r$Significant, ]
  data.frame(Feature = f,
             Significant_pairs = sprintf("%d/%d", nrow(s), nrow(r)),
             Rho_range = sprintf("%.2f to %.2f", min(r$rho), max(r$rho)),
             Positive_with = paste(s$Gene[s$rho > 0], collapse = ", "),
             Negative_with = paste(s$Gene[s$rho < 0], collapse = ", "))
}))
print(summary_tbl, row.names = FALSE)

write.csv(res[order(res$q_BH), ], out_path("radiomic_gene_correlations_all_pairs.csv"), row.names = FALSE)
write.csv(summary_tbl, out_path("radiomic_gene_correlations_summary.csv"), row.names = FALSE)
cat("Saved correlation tables to", output_dir, "\n")
