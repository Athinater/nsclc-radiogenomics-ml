# ==============================================================================
# common.R -- shared configuration and helper functions
# Run every script from the repository root, e.g.  Rscript scripts/01_radiomic_h2o.R
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---- Paths (relative to the repository root) ---------------------------------
data_dir           <- "data"
output_dir         <- "output"
counts_file        <- file.path(data_dir, "GSE103584_raw_counts_GRCh38.p13_NCBI.tsv")  # GEO GSE103584 raw counts
clinical_rad_file  <- file.path(data_dir, "clinicaldatawithout.csv")  # radiomic features + histology (semicolon-separated)
clinical_rg_file   <- file.path(data_dir, "clinicaldata.csv")         # radiomic features + histology, patients with matched RNA-Seq
clinical_meta_file <- file.path(data_dir, "onlyclinical.csv")         # clinical metadata used as DESeq2 colData
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(f) file.path(output_dir, f)

# ---- Class coding --------------------------------------------------------------
# Histology is coded 0 = adenocarcinoma (ADC), 1 = squamous cell carcinoma (SCC).
# Sensitivity, specificity, MCC and F1 are computed with ADC as the POSITIVE class,
# as reported in the manuscript (Tables 2 and 3). All probabilities are written to
# the prediction files as both prob_SCC (= P(class 1)) and prob_ADC (= 1 - prob_SCC).
adc_positive_metrics <- function(pred01, truth01) {
  tp <- sum(pred01 == 0 & truth01 == 0)   # ADC predicted as ADC
  fn <- sum(pred01 == 1 & truth01 == 0)   # ADC predicted as SCC
  tn <- sum(pred01 == 1 & truth01 == 1)   # SCC predicted as SCC
  fp <- sum(pred01 == 0 & truth01 == 1)   # SCC predicted as ADC
  denom <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  list(Sensitivity = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
       Specificity = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp),
       MCC         = if (denom == 0) NA_real_ else (tp * tn - fp * fn) / denom,
       F1          = if ((2 * tp + fp + fn) == 0) NA_real_ else 2 * tp / (2 * tp + fp + fn),
       TP = tp, FN = fn, TN = tn, FP = fp)
}

write_preds <- function(file, truth01, prob_scc, scenario) {
  write.csv(data.frame(true_label = truth01,            # 0 = ADC, 1 = SCC
                       prob_SCC   = as.numeric(prob_scc),
                       prob_ADC   = 1 - as.numeric(prob_scc),
                       scenario   = scenario),
            file, row.names = FALSE)
}

# ---- RNA-Seq counts ---------------------------------------------------------------
# The GEO file contains Entrez GeneIDs only, so symbols are mapped with org.Hs.eg.db.
# Genes without a symbol are dropped; duplicated symbols are made unique.
# Pre-filter: >= 10 counts in >= 2 samples.
load_counts <- function(file) {
  tbl      <- fread(file, header = TRUE, data.table = FALSE)
  gene_ids <- tbl[, 1]
  symbols  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(gene_ids),
                                    column = "SYMBOL", keytype = "ENTREZID",
                                    multiVals = "first")
  tbl     <- as.matrix(tbl[, -1])
  keep    <- rowSums(tbl >= 10) >= 2
  tbl     <- tbl[keep, , drop = FALSE]
  symbols <- symbols[keep]
  valid   <- !is.na(symbols)
  cat("Genes dropped (no symbol):", sum(!valid), "of", length(valid), "after count filter\n")
  tbl     <- tbl[valid, , drop = FALSE]
  rownames(tbl) <- make.unique(unname(symbols[valid]))
  colnames(tbl) <- toupper(trimws(colnames(tbl)))
  cat("Count matrix:", nrow(tbl), "genes x", ncol(tbl), "samples\n")
  tbl
}

# ---- Consistency check of the selected genes against the manuscript -------------
manuscript_genes <- c("KRT5", "KRT6A", "KRT6C", "KRT16", "KRT13", "DSG3", "DSC3",
                      "SERPINB5", "CLCA2", "AKR1B10",                     # SCC-associated
                      "NAPSA", "CEACAM5", "CEACAM6", "MUC21", "PIGR", "SFTA2",
                      "SLC34A2", "C16orf89")                              # ADC-associated
check_selected_genes <- function(selected_genes) {
  if (setequal(selected_genes, manuscript_genes)) {
    cat("Selected genes match the 18 genes reported in the manuscript.\n")
  } else {
    warning("Selected genes differ from the manuscript list.\n  Missing: ",
            paste(setdiff(manuscript_genes, selected_genes), collapse = ", "),
            "\n  Extra:   ", paste(setdiff(selected_genes, manuscript_genes), collapse = ", "))
  }
}

# ---- Gene selection INSIDE a cross-validation fold ---------------------------------
# For the radiogenomic and genomic pipelines the 18 genes (top-10 |loadings| of PC1 and
# PC2 of the variance-stabilised training matrix) are re-selected within every fold using
# the fold-training samples only; the fold-validation samples are transformed with the
# fold-training dispersions (test-specific size factors), exactly as for the final test set.
fold_gene_matrices <- function(tbl_tr, coldata_tr, tbl_va, coldata_va, n_top = 10) {
  dds <- DESeqDataSetFromMatrix(tbl_tr, coldata_tr, ~ Histology)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, quiet = TRUE)
  vsd <- vst(dds, blind = FALSE)
  rot <- prcomp(t(assay(vsd)))$rotation
  genes <- unique(c(rownames(rot)[head(order(abs(rot[, "PC1"]), decreasing = TRUE), n_top)],
                    rownames(rot)[head(order(abs(rot[, "PC2"]), decreasing = TRUE), n_top)]))
  dds_va <- DESeqDataSetFromMatrix(tbl_va, coldata_va, ~ Histology)
  dds_va <- estimateSizeFactors(dds_va)
  dispersions(dds_va) <- dispersions(dds)[match(rownames(dds_va), rownames(dds))]
  vsd_va <- varianceStabilizingTransformation(dds_va, blind = FALSE)
  list(genes = genes,
       X_tr  = t(assay(vsd)[genes, , drop = FALSE]),
       X_va  = t(assay(vsd_va)[genes, , drop = FALSE]))
}

# ctx = list(ids        = sample IDs of the training rows of X_train_norm (same order),
#            tbl_train  = raw count matrix of the training samples,
#            coldata_train = colData of the training samples,
#            radio_cols = names of the radiomic columns of X_train_norm (character(0) for genomics))
build_fold_matrices <- function(X_train_norm, tr_idx, va_idx, ctx) {
  ids_tr <- ctx$ids[tr_idx]; ids_va <- ctx$ids[va_idx]
  g <- fold_gene_matrices(ctx$tbl_train[, ids_tr, drop = FALSE], ctx$coldata_train[ids_tr, , drop = FALSE],
                          ctx$tbl_train[, ids_va, drop = FALSE], ctx$coldata_train[ids_va, , drop = FALSE])
  genes  <- intersect(colnames(g$X_tr), colnames(g$X_va))
  G_tr   <- g$X_tr[, genes, drop = FALSE]
  G_va   <- g$X_va[, genes, drop = FALSE]
  mu     <- colMeans(G_tr)
  sdv    <- apply(G_tr, 2, sd); sdv[sdv == 0 | is.na(sdv)] <- 1
  G_tr   <- scale(G_tr, center = mu, scale = sdv)
  G_va   <- scale(G_va, center = mu, scale = sdv)
  colnames(G_tr) <- colnames(G_va) <- make.names(genes)
  if (length(ctx$radio_cols) > 0) {
    X_f_tr <- cbind(X_train_norm[tr_idx, ctx$radio_cols, drop = FALSE], G_tr)
    X_f_va <- cbind(X_train_norm[va_idx, ctx$radio_cols, drop = FALSE], G_va)
  } else {
    X_f_tr <- G_tr
    X_f_va <- G_va
  }
  list(X_f_tr = X_f_tr, X_f_va = X_f_va, genes = genes)
}
