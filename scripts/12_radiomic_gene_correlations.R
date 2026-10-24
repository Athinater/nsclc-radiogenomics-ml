# ==============================================================================
# 12_radiomic_gene_correlations.R
# Spearman correlations between the 10 most frequently selected radiomic features and
# the 18 PCA-selected genes across the full cohort (n = 113; 10 x 18 = 180 pairs),
# Benjamini-Hochberg FDR across all 180 tests (significant at q < 0.05).
# Exploratory: this analysis does not feed back into model training.
# Gene expression: DESeq2 variance-stabilized counts of the full cohort.
# ==============================================================================
suppressPackageStartupMessages({ library(DESeq2) })
source("scripts/common.R")

# The 10 radiomic features selected most often across all pipeline x scenario x classifier
# combinations (feature-selection results of scripts 01-04 and 07).
top_radiomic_features <- c(
  "NGTDM_Coarseness", "GLCM_ClusterShade", "IHMinimumGreyLevel", "IBCoefficientOfVariation",
  "GLSZM_SmallZoneEmphasis", "IH50thPercentile", "MedianIntensity", "IHMedian",
  "50thIntensityPercentile", "IBQuartileCoefficientOfDispersion")

genes <- manuscript_genes     # the 18 genes (first 10 SCC-associated, last 8 ADC-associated)

# ---- Data ---------------------------------------------------------------------------
tbl <- load_counts(counts_file)
clinical_data <- read.csv(clinical_rg_file, sep = ";", header = TRUE, row.names = 1)
coldata <- read.csv(clinical_meta_file, header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1)
rownames(coldata) <- toupper(trimws(rownames(coldata)))
coldata <- coldata[coldata$Histology %in% c("Adenocarcinoma", "Squamous cell carcinoma"), , drop = FALSE]

ids <- intersect(intersect(toupper(trimws(rownames(clinical_data))), rownames(coldata)), colnames(tbl))
cat("Patients with radiomic features, clinical metadata and RNA-Seq:", length(ids), "(manuscript: 113)\n")
if (length(ids) != 113) warning("Cohort size differs from the 113 patients reported in the manuscript.")
rownames(clinical_data) <- toupper(trimws(rownames(clinical_data)))

# Radiomic columns: match names with or without the INTENSITY.BASED_/INTENSITY.HISTOGRAM_/X prefix
strip_prefix <- function(x) sub("^(INTENSITY\\.BASED_|INTENSITY\\.HISTOGRAM_|X)", "", x)
cols <- vapply(top_radiomic_features, function(f) {
  hit <- colnames(clinical_data)[colnames(clinical_data) == f | strip_prefix(colnames(clinical_data)) == f]
  if (length(hit) != 1) stop("Radiomic feature not found (or ambiguous) in the data: ", f)
  hit
}, character(1))
rad <- as.matrix(clinical_data[ids, cols, drop = FALSE]); colnames(rad) <- top_radiomic_features

# VST expression of the 18 genes (full cohort)
dds <- DESeqDataSetFromMatrix(tbl[, ids], coldata[ids, , drop = FALSE], ~ Histology)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, quiet = TRUE)   # required for vst(blind = FALSE)
vsd <- vst(dds, blind = FALSE)
missing_genes <- setdiff(genes, rownames(vsd))
if (length(missing_genes) > 0) stop("Genes not found in the expression data: ", paste(missing_genes, collapse = ", "))
expr <- t(assay(vsd)[genes, , drop = FALSE])

# ---- Spearman correlations ------------------------------------------------------------
grid <- expand.grid(Feature = top_radiomic_features, Gene = genes, stringsAsFactors = FALSE)
tests <- lapply(seq_len(nrow(grid)), function(i) {
  ct <- suppressWarnings(cor.test(rad[, grid$Feature[i]], expr[, grid$Gene[i]],
                                  method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value)
})
res <- cbind(grid, do.call(rbind, tests))
res$q_BH <- p.adjust(res$p, method = "BH")
res$Significant <- res$q_BH < 0.05
res$Gene_group <- ifelse(res$Gene %in% genes[1:10], "SCC-associated", "ADC-associated")

cat("\nPairs tested:", nrow(res), "| significant at q < 0.05:", sum(res$Significant), "(manuscript: 180 / 34)\n")

summary_tbl <- do.call(rbind, lapply(top_radiomic_features, function(f) {
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
