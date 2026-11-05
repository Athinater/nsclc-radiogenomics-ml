# ==============================================================================
# common.R -- configuration, data preparation, feature selection and the nested
# cross-validation runner shared by the model scripts (01-07).
# Run every script from the repository root, e.g.  Rscript scripts/01_radiomic_h2o.R
#
# Packages are called with pkg:: so that this file can also be sourced by the
# analysis scripts. Required: caret, glmnet, mRMRe, pROC, data.table; for the
# radiogenomic/genomic pipelines also DESeq2, SummarizedExperiment, S4Vectors,
# AnnotationDbi and org.Hs.eg.db.
# ==============================================================================

# ---- Paths (relative to the repository root) ---------------------------------
data_dir           <- "data"
output_dir         <- "output"
counts_file        <- file.path(data_dir, "GSE103584_raw_counts_GRCh38.p13_NCBI.tsv")  # GEO GSE103584 raw counts
clinical_rad_file  <- file.path(data_dir, "clinicaldatawithout.csv")  # radiomic features + histology (semicolon-separated)
clinical_rg_file   <- file.path(data_dir, "clinicaldata.csv")         # radiomic features + histology, patients with matched RNA-Seq
clinical_meta_file <- file.path(data_dir, "onlyclinical.csv")         # clinical metadata used as DESeq2 colData
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(f) file.path(output_dir, f)

# ---- Analysis settings (fixed before any model is evaluated) -----------------
seed_value         <- 123
train_fraction     <- 0.70    # stratified train/test split
n_folds            <- 5       # outer cross-validation folds
inner_holdout      <- 0.20    # stratified inner hold-out for hyperparameter selection
decision_threshold <- 0.5     # same fixed threshold on P(SCC) for every classifier
n_top_genes        <- 10      # genes per principal component (PC1, PC2)
pc_variance        <- 0.90    # SC: components explaining >= 90% of the variance ...
pc_cap             <- 10      # ... capped at 10 components
n_fallback         <- 5       # LASSO fallback: top-5 features by |correlation| with the class
n_mrmr_min         <- 5       # mRMR selects max(#LASSO features, 5) features
n_force_radiomic   <- 5       # radiogenomic SA/SB/SC always contain the top-5 radiomic features

scenarios_all     <- c("S0", "SA", "SB", "SC")
scenarios_to_run  <- scenarios_all          # e.g. c("S0", "SA", "SB") to drop scenario SC everywhere
scenario_labels   <- c(S0 = "All features", SA = "LASSO", SB = "LASSO+mRMR",
                       SC = "LASSO+mRMR+PCA-loading features")
radiomic_pattern  <- "MORPHOLOGICAL_|INTENSITY.BASED_|INTENSITY.HISTOGRAM_|GLCM_|GLRLM_|GLSZM_|NGTDM_"

# ---- Class coding ----------------------------------------------------------------
# Histology is coded 0 = adenocarcinoma (ADC), 1 = squamous cell carcinoma (SCC).
# Models output P(SCC). ADC is the POSITIVE class for every reported metric (AUC,
# sensitivity, specificity, MCC, F1) and for the post-hoc analyses (DeLong, calibration,
# decision curves). Prediction files contain prob_SCC and prob_ADC = 1 - prob_SCC.
norm_id <- function(x) toupper(trimws(x))

# Predicted class (0 = ADC, 1 = SCC) from P(SCC) with the fixed threshold
classify <- function(p_scc) as.integer(p_scc >= decision_threshold)

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

# ROC with ADC as the positive class and a FIXED direction (P(ADC) higher for ADC).
# The direction is never chosen from the data, so a model worse than chance gives AUC < 0.5.
roc_adc <- function(truth01, p_scc) {
  pROC::roc(response = as.integer(truth01 == 0), predictor = 1 - as.numeric(p_scc),
            levels = c(0, 1), direction = "<", quiet = TRUE)
}
auc_adc <- function(truth01, p_scc) as.numeric(pROC::auc(roc_adc(truth01, p_scc)))
auc_ci_adc <- function(truth01, p_scc) {          # AUC with DeLong 95% CI, truncated to [0, 1]
  r  <- roc_adc(truth01, p_scc)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  list(AUC = as.numeric(pROC::auc(r)), CI_Low = max(ci[1], 0), CI_High = min(ci[3], 1))
}

# Test-set evaluation: AUC (DeLong CI) and threshold metrics, all with ADC positive
evaluate_test <- function(truth01, p_scc) {
  c(auc_ci_adc(truth01, p_scc), adc_positive_metrics(classify(p_scc), truth01))
}

write_preds <- function(file, ids, truth01, prob_scc, scenario) {
  stopifnot(length(ids) == length(truth01), length(truth01) == length(prob_scc))
  write.csv(data.frame(patient_id = norm_id(ids),
                       true_label = truth01,            # 0 = ADC, 1 = SCC
                       prob_SCC   = as.numeric(prob_scc),
                       prob_ADC   = 1 - as.numeric(prob_scc),
                       scenario   = scenario),
            file, row.names = FALSE)
}

# ---- Input data --------------------------------------------------------------------
read_clinical <- function(file) {
  d <- read.csv(file, sep = ";", header = TRUE, row.names = 1)
  rownames(d) <- norm_id(rownames(d))
  h   <- factor(d$Histology)
  lev <- levels(h)
  if (length(lev) != 2 || !grepl("^(0|adc|aden)", tolower(lev[1])))
    stop("Unexpected Histology coding in ", file, ": ", paste(lev, collapse = " / "),
         ". Expected 0 = ADC and 1 = SCC (or Adenocarcinoma / Squamous cell carcinoma).")
  list(data = d, y = as.integer(h) - 1L)   # 0 = ADC, 1 = SCC
}

# RNA-Seq counts. The GEO file contains Entrez GeneIDs only, so symbols are mapped with
# org.Hs.eg.db. Genes without a symbol are dropped; duplicated symbols are made unique.
# Pre-filter: >= 10 counts in >= 2 samples. This filter does not use the class labels and is
# applied once to all samples before the train/test split.
load_counts <- function(file) {
  tbl      <- data.table::fread(file, header = TRUE, data.table = FALSE)
  gene_ids <- tbl[, 1]
  symbols  <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(gene_ids),
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
  colnames(tbl) <- norm_id(colnames(tbl))
  cat("Count matrix:", nrow(tbl), "genes x", ncol(tbl), "samples\n")
  tbl
}

split_train_test <- function(y) {
  set.seed(seed_value)
  caret::createDataPartition(factor(y), p = train_fraction, list = FALSE)[, 1]
}

# Near-zero-variance filter and z-scaling with training statistics only
nzv_scale <- function(train_raw, test_raw) {
  train_raw <- as.data.frame(train_raw, check.names = FALSE)
  test_raw  <- as.data.frame(test_raw,  check.names = FALSE)
  nzv <- caret::nearZeroVar(train_raw)
  if (length(nzv) > 0) {
    cat("Near-zero-variance features removed:", paste(colnames(train_raw)[nzv], collapse = ", "), "\n")
    train_raw <- train_raw[, -nzv, drop = FALSE]
    test_raw  <- test_raw[,  -nzv, drop = FALSE]
  }
  mu  <- vapply(train_raw, mean, numeric(1))
  sdv <- vapply(train_raw, sd,   numeric(1))
  sdv[sdv == 0 | is.na(sdv)] <- 1
  list(train = scale(as.matrix(train_raw), center = mu, scale = sdv),
       test  = scale(as.matrix(test_raw),  center = mu, scale = sdv))
}

# ---- Gene selection and variance-stabilising transformation (DESeq2) -------------------
# Training samples: size factors and dispersions are estimated and the VST is computed exactly
# as DESeq2::vst(blind = FALSE) does it (dispersion trend fitted on 1000 genes with mean
# normalised count > 5); the 18 genes are the top-10 |loadings| of PC1 and PC2.
# New samples (validation fold / test set) get their own size factors and are transformed with
# the SAME training dispersion trend ("frozen" VST, see ?DESeq2::varianceStabilizingTransformation):
# their labels and their own dispersions are never used.

# The dispersion trend that DESeq2::vst() fits internally (same code as in DESeq2::vst)
vst_dispersion_trend <- function(dds, nsub = 1000) {
  base_mean <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
  sub <- dds[base_mean > 5, ]
  bm  <- base_mean[base_mean > 5]
  o   <- order(bm)
  sub <- sub[o[round(seq(from = 1, to = length(o), length.out = nsub))], ]
  sub <- DESeq2::estimateDispersionsGeneEst(sub, quiet = TRUE)
  sub <- DESeq2::estimateDispersionsFit(sub, fitType = "parametric", quiet = TRUE)
  DESeq2::dispersionFunction(sub)
}

vst_select_genes <- function(tbl_tr, coldata_tr, n_top = n_top_genes) {
  coldata_tr$Histology <- factor(make.names(coldata_tr$Histology))
  dds <- DESeq2::DESeqDataSetFromMatrix(tbl_tr, coldata_tr, ~ Histology)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
  ref <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))   # transformation of the original analysis
  suppressMessages(DESeq2::dispersionFunction(dds) <- vst_dispersion_trend(dds))
  vst_tr <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
  if (max(abs(vst_tr - ref)) > 1e-6)
    stop("Internal check failed: the frozen VST does not reproduce DESeq2::vst() on the training samples.")
  pca <- prcomp(t(vst_tr))
  rot <- pca$rotation
  pc1 <- rownames(rot)[head(order(abs(rot[, 1]), decreasing = TRUE), n_top)]
  pc2 <- rownames(rot)[head(order(abs(rot[, 2]), decreasing = TRUE), n_top)]
  list(dds = dds, vst = vst_tr, genes = unique(c(pc1, pc2)), top_pc1 = pc1, top_pc2 = pc2,
       var_explained = summary(pca)$importance[2, 1:2])
}

apply_vst <- function(tbl_new, dds_train) {
  dds_new <- DESeq2::DESeqDataSetFromMatrix(tbl_new, S4Vectors::DataFrame(row.names = colnames(tbl_new)), ~ 1)
  dds_new <- DESeq2::estimateSizeFactors(dds_new)                       # own size factors, no labels
  suppressMessages(DESeq2::dispersionFunction(dds_new) <- DESeq2::dispersionFunction(dds_train))
  SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds_new, blind = FALSE))
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

# ---- Pipeline preparation --------------------------------------------------------------
# Each prepare_* function returns a list with the scaled training/test matrices, the 0/1
# labels, the patient IDs, and the cross-validation folds with their fold-level matrices
# (computed once and reused by every scenario and classifier).
prepare_radiomic <- function() {
  cl    <- read_clinical(clinical_rad_file)
  feats <- grep(radiomic_pattern, colnames(cl$data), value = TRUE)
  idx   <- split_train_test(cl$y)
  pp    <- nzv_scale(cl$data[idx, feats], cl$data[-idx, feats])
  finish_prep(list(pipeline = "radiomic",
                   X_train = pp$train, X_test = pp$test,
                   y_train = cl$y[idx], y_test = cl$y[-idx],
                   ids_train = rownames(cl$data)[idx], ids_test = rownames(cl$data)[-idx],
                   force_cols = character(0)))
}

# Radiogenomic and genomic pipelines share the split, the genes and the CV folds.
# Returns a named list with one prepared pipeline per requested name.
prepare_omics <- function(pipelines = c("radiogenomic", "genomic")) {
  pipelines <- match.arg(pipelines, several.ok = TRUE)
  tbl  <- load_counts(counts_file)
  meta <- read.csv(clinical_meta_file, header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1)
  rownames(meta) <- norm_id(rownames(meta))
  meta <- meta[meta$Histology %in% c("Adenocarcinoma", "Squamous cell carcinoma"), , drop = FALSE]
  with_rna <- intersect(rownames(meta), colnames(tbl))

  cl  <- read_clinical(clinical_rg_file)
  idx <- split_train_test(cl$y)
  ids_train <- intersect(rownames(cl$data)[idx],  with_rna)
  ids_test  <- intersect(rownames(cl$data)[-idx], with_rna)
  lost <- setdiff(rownames(cl$data), c(ids_train, ids_test))
  if (length(lost) > 0) warning("Patients without RNA-Seq/metadata, excluded: ", paste(lost, collapse = ", "))
  y_train <- cl$y[match(ids_train, rownames(cl$data))]
  y_test  <- cl$y[match(ids_test,  rownames(cl$data))]
  meta_y  <- as.integer(meta[c(ids_train, ids_test), "Histology"] == "Squamous cell carcinoma")
  if (any(meta_y != c(y_train, y_test)))
    warning("Histology differs between ", clinical_rg_file, " and ", clinical_meta_file,
            " for: ", paste(c(ids_train, ids_test)[meta_y != c(y_train, y_test)], collapse = ", "))

  cat("\n=== Gene selection: DESeq2 VST + PCA on the training set only ===\n")
  gs <- vst_select_genes(tbl[, ids_train, drop = FALSE], meta[ids_train, , drop = FALSE])
  check_selected_genes(gs$genes)
  cat("PC1 top", n_top_genes, ":", paste(gs$top_pc1, collapse = ", "), "\n")
  cat("PC2 top", n_top_genes, ":", paste(gs$top_pc2, collapse = ", "), "\n")
  cat(sprintf("Variance explained: PC1 %.1f%%, PC2 %.1f%%\n", 100 * gs$var_explained[1], 100 * gs$var_explained[2]))
  # subtype association of each gene (training samples): mean VST expression at least 1 unit
  # (about 2-fold) higher in SCC or in ADC; otherwise "Not subtype-associated"
  y_tr_rna  <- cl$y[match(colnames(gs$vst), rownames(cl$data))]
  mean_adc  <- rowMeans(gs$vst[gs$genes, y_tr_rna == 0, drop = FALSE])
  mean_scc  <- rowMeans(gs$vst[gs$genes, y_tr_rna == 1, drop = FALSE])
  assoc     <- ifelse(mean_scc - mean_adc >= 1, "SCC-associated",
               ifelse(mean_adc - mean_scc >= 1, "ADC-associated", "Not subtype-associated"))
  write.csv(data.frame(Gene = gs$genes,
                       Top_PC1 = gs$genes %in% gs$top_pc1, Top_PC2 = gs$genes %in% gs$top_pc2,
                       Mean_VST_ADC = round(mean_adc, 3), Mean_VST_SCC = round(mean_scc, 3),
                       Association = assoc,
                       In_manuscript_list = gs$genes %in% manuscript_genes),
            out_path("selected_genes.csv"), row.names = FALSE)

  G_tr <- t(gs$vst[gs$genes, , drop = FALSE])
  G_te <- t(apply_vst(tbl[, ids_test, drop = FALSE], gs$dds)[gs$genes, , drop = FALSE])
  colnames(G_tr) <- colnames(G_te) <- make.names(gs$genes)
  radio_feats <- grep(radiomic_pattern, colnames(cl$data), value = TRUE)

  # fold-level gene matrices: genes re-selected inside every fold (shared by both pipelines)
  folds <- make_folds(y_train)
  fold_genes <- lapply(seq_along(folds), function(k) {
    va <- folds[[k]]; tr <- setdiff(seq_along(y_train), va)
    tryCatch(fold_gene_matrices(tbl[, ids_train[tr], drop = FALSE], meta[ids_train[tr], , drop = FALSE],
                                tbl[, ids_train[va], drop = FALSE]),
             error = function(e) { warning("DESeq2 failed in fold ", k, ": ", conditionMessage(e)); NULL })
  })

  out <- list()
  for (p in pipelines) {
    if (p == "radiogenomic") {
      raw_tr <- cbind(as.matrix(cl$data[ids_train, radio_feats, drop = FALSE]), G_tr)
      raw_te <- cbind(as.matrix(cl$data[ids_test,  radio_feats, drop = FALSE]), G_te)
    } else {
      raw_tr <- G_tr; raw_te <- G_te
    }
    pp <- nzv_scale(raw_tr, raw_te)
    radio_cols <- intersect(radio_feats, colnames(pp$train))
    prep <- list(pipeline = p, X_train = pp$train, X_test = pp$test,
                 y_train = y_train, y_test = y_test, ids_train = ids_train, ids_test = ids_test,
                 force_cols = radio_cols, genes = gs$genes)
    prep$cv <- lapply(seq_along(folds), function(k) {
      va <- folds[[k]]; tr <- setdiff(seq_along(y_train), va)
      g  <- fold_genes[[k]]
      if (is.null(g)) return(list(tr = tr, va = va, X_tr = NULL, X_va = NULL))
      X_tr <- g$X_tr; X_va <- g$X_va
      if (length(radio_cols) > 0) {      # radiomic columns re-scaled with fold-training statistics
        r <- fold_scale(pp$train[tr, radio_cols, drop = FALSE], pp$train[va, radio_cols, drop = FALSE])
        X_tr <- cbind(r$X_tr, X_tr)
        X_va <- cbind(r$X_va, X_va)
      }
      list(tr = tr, va = va, X_tr = X_tr, X_va = X_va, genes = g$genes)
    })
    out[[p]] <- finish_prep(prep)
  }
  out
}

make_folds <- function(y) {
  set.seed(seed_value)
  caret::createFolds(factor(y), k = n_folds, list = TRUE)   # stratified by class
}

# Genes re-selected from the fold-training samples; fold-validation samples transformed with
# the fold-training dispersion trend; z-scaled with fold-training statistics.
fold_gene_matrices <- function(tbl_tr, coldata_tr, tbl_va) {
  gs   <- vst_select_genes(tbl_tr, coldata_tr)
  G_tr <- t(gs$vst[gs$genes, , drop = FALSE])
  G_va <- t(apply_vst(tbl_va, gs$dds)[gs$genes, , drop = FALSE])
  mu   <- colMeans(G_tr)
  sdv  <- apply(G_tr, 2, sd); sdv[sdv == 0 | is.na(sdv)] <- 1
  G_tr <- scale(G_tr, center = mu, scale = sdv)
  G_va <- scale(G_va, center = mu, scale = sdv)
  colnames(G_tr) <- colnames(G_va) <- make.names(gs$genes)
  list(X_tr = G_tr, X_va = G_va, genes = gs$genes)
}

finish_prep <- function(prep) {
  stopifnot(nrow(prep$X_train) == length(prep$y_train), nrow(prep$X_test) == length(prep$y_test),
            length(prep$ids_test) == length(prep$y_test))
  cat(sprintf("\n[%s] training n = %d (ADC %d, SCC %d) | test n = %d (ADC %d, SCC %d) | features = %d\n",
              prep$pipeline, length(prep$y_train), sum(prep$y_train == 0), sum(prep$y_train == 1),
              length(prep$y_test), sum(prep$y_test == 0), sum(prep$y_test == 1), ncol(prep$X_train)))
  writeLines(sort(prep$ids_test), out_path(paste0(prep$pipeline, "_test_ids.txt")))
  if (is.null(prep$cv)) {                   # radiomic: rows of the fold, re-scaled with fold-training statistics
    folds <- make_folds(prep$y_train)
    prep$cv <- lapply(folds, function(va) {
      tr <- setdiff(seq_along(prep$y_train), va)
      r  <- fold_scale(prep$X_train[tr, , drop = FALSE], prep$X_train[va, , drop = FALSE])
      list(tr = tr, va = va, X_tr = r$X_tr, X_va = r$X_va)
    })
  }
  prep
}

# z-scaling with the statistics of the fold-training patients only
fold_scale <- function(X_tr, X_va) {
  mu  <- colMeans(X_tr)
  sdv <- apply(X_tr, 2, sd); sdv[sdv == 0 | is.na(sdv)] <- 1
  list(X_tr = scale(X_tr, center = mu, scale = sdv), X_va = scale(X_va, center = mu, scale = sdv))
}

# ---- Feature selection (identical for every classifier) ------------------------------------
# S0 = all features; SA = LASSO (cv.glmnet, alpha = 1, 5-fold, lambda.min; fallback: top-5 by
# |correlation|); SB = SA + mRMR (max(#SA, 5) features); SC = SB + the features with the largest
# absolute loading on the first n principal components explaining >= 90% of the variance (n <= 10).
# Radiogenomic pipeline: the 5 radiomic features most correlated with the class are added to SA
# (and hence to SB and SC), so that every radiogenomic model contains imaging features.
top_cor <- function(X, y, cols = colnames(X), n) {
  cols <- intersect(cols, colnames(X))
  if (length(cols) == 0) return(character(0))
  r <- abs(cor(X[, cols, drop = FALSE], y)[, 1])
  names(sort(r, decreasing = TRUE))[seq_len(min(n, sum(!is.na(r))))]
}

select_features <- function(X, y, scenario, force_cols = character(0), seed = seed_value) {
  if (scenario == "S0") return(colnames(X))
  set.seed(seed)
  lc <- glmnet::cv.glmnet(as.matrix(X), y, family = "binomial", alpha = 1, nfolds = 5)
  co <- as.matrix(coef(lc, s = lc$lambda.min))
  f  <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
  if (length(f) == 0) f <- top_cor(X, y, n = n_fallback)
  if (length(force_cols) > 0) f <- unique(c(f, top_cor(X, y, force_cols, n_force_radiomic)))
  if (scenario == "SA") return(f)
  n_mrmr <- min(max(length(f), n_mrmr_min), ncol(X))
  dm <- mRMRe::mRMR.data(data = data.frame(X, response = as.ordered(y)))
  mr <- mRMRe::mRMR.classic(data = dm, target_indices = ncol(X) + 1, feature_count = n_mrmr)
  f  <- unique(c(f, colnames(X)[unlist(mr@filters)]))
  if (scenario == "SB") return(f)
  pc   <- prcomp(X, scale. = FALSE)
  pvar <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  nc   <- min(which(pvar >= pc_variance)[1], pc_cap)
  pcon <- apply(abs(pc$rotation[, seq_len(nc), drop = FALSE]), 1, max)
  unique(c(f, names(sort(pcon, decreasing = TRUE))[seq_len(nc)]))
}

# ---- Hyperparameter selection and refit -------------------------------------------------------
# configs: list of configurations (each with a $label); fit_predict(X_a, y_a, X_b, cfg) fits a
# model on (X_a, y_a) and returns P(SCC) for X_b.
# The configurations are compared on a stratified inner 20% hold-out of the data passed in
# (inner AUC); the best one (first in grid order if tied) is refitted on ALL the data passed in
# and used to predict X_new. Used identically inside every CV fold and for the final model.
tune_refit_predict <- function(X, y, X_new, configs, fit_predict) {
  set.seed(seed_value)
  va <- caret::createDataPartition(factor(y), p = inner_holdout, list = FALSE)[, 1]
  inner <- vapply(configs, function(cfg) {
    p <- tryCatch(fit_predict(X[-va, , drop = FALSE], y[-va], X[va, , drop = FALSE], cfg),
                  error = function(e) { cat("    ", cfg$label, "failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(p)) NA_real_ else auc_adc(y[va], p)
  }, numeric(1))
  if (all(is.na(inner))) stop("no configuration could be fitted")
  best <- which.max(inner)
  list(prob = fit_predict(X, y, X_new, configs[[best]]), label = configs[[best]]$label,
       inner_auc = inner[best])
}

# ---- Nested cross-validation and final evaluation for one scenario ---------------------------
run_nested <- function(prep, classifier, scenario, make_configs, fit_predict, run_tag) {
  cat(sprintf("\n=== %s | %s | %s (%s) ===\n", prep$pipeline, classifier, scenario, scenario_labels[[scenario]]))
  cv_aucs   <- rep(NA_real_, length(prep$cv))
  fold_rows <- list()
  for (k in seq_along(prep$cv)) {
    fd <- prep$cv[[k]]
    if (is.null(fd$X_tr)) next
    y_tr  <- prep$y_train[fd$tr]; y_va <- prep$y_train[fd$va]
    feats <- select_features(fd$X_tr, y_tr, scenario, prep$force_cols, seed = seed_value + k)
    res   <- tryCatch(tune_refit_predict(fd$X_tr[, feats, drop = FALSE], y_tr, fd$X_va[, feats, drop = FALSE],
                                         make_configs(length(feats)), fit_predict),
                      error = function(e) { cat("  fold", k, "failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(res)) cv_aucs[k] <- auc_adc(y_va, res$prob)
    cat(sprintf("  fold %d: %3d features | %-34s | fold AUC %.3f\n", k, length(feats),
                if (is.null(res)) "-" else res$label, cv_aucs[k]))
    fold_rows[[k]] <- data.frame(Pipeline = prep$pipeline, Classifier = classifier, Scenario = scenario,
                                 Fold = k, Features = length(feats),
                                 Selected_config = if (is.null(res)) NA_character_ else res$label,
                                 Inner_AUC = if (is.null(res)) NA_real_ else res$inner_auc,
                                 Fold_AUC = cv_aucs[k], Feature_list = paste(feats, collapse = ";"))
  }

  # Final model: features selected and configuration tuned on the full training set, refitted on
  # all training patients, evaluated once on the test set.
  feats <- select_features(prep$X_train, prep$y_train, scenario, prep$force_cols, seed = seed_value)
  res   <- tune_refit_predict(prep$X_train[, feats, drop = FALSE], prep$y_train,
                              prep$X_test[, feats, drop = FALSE], make_configs(length(feats)), fit_predict)
  ev <- evaluate_test(prep$y_test, res$prob)
  preds_file <- sprintf("%s_preds_%s.csv", run_tag, scenario)
  write_preds(out_path(preds_file), prep$ids_test, prep$y_test, res$prob, scenario)
  valid <- cv_aucs[!is.na(cv_aucs)]
  cv    <- if (length(valid) > 0) mean(valid) else NA_real_
  row <- data.frame(Pipeline = prep$pipeline, Classifier = classifier, Scenario = scenario,
                    Features = length(feats), Valid_Folds = length(valid),
                    CV_AUC = cv, CV_AUC_SD = if (length(valid) > 1) sd(valid) else NA_real_,
                    Test_AUC = ev$AUC, CI_Low = ev$CI_Low, CI_High = ev$CI_High,
                    Overfitting_Gap = cv - ev$AUC,
                    Sensitivity = ev$Sensitivity, Specificity = ev$Specificity, MCC = ev$MCC, F1 = ev$F1,
                    TP = ev$TP, FN = ev$FN, TN = ev$TN, FP = ev$FP,
                    Threshold = decision_threshold, Selected_config = res$label,
                    Run_tag = run_tag, Preds_file = preds_file,
                    Feature_list = paste(feats, collapse = ";"), stringsAsFactors = FALSE)
  cat(sprintf("  FINAL: %d features | %s | CV-AUC %.3f | test AUC %.3f [%.3f-%.3f] | Sens %.3f Spec %.3f MCC %.3f F1 %.3f\n",
              row$Features, row$Selected_config, row$CV_AUC, row$Test_AUC, row$CI_Low, row$CI_High,
              row$Sensitivity, row$Specificity, row$MCC, row$F1))
  list(result = row, folds = do.call(rbind, fold_rows))
}

# grids: named list (one entry per scenario) of functions n_features -> list of configurations
run_all_scenarios <- function(prep, classifier, grids, fit_predict, run_tag,
                              scenarios = scenarios_to_run, cleanup = function() invisible(NULL)) {
  out <- lapply(scenarios, function(s) {
    r <- run_nested(prep, classifier, s, grids[[s]], fit_predict, run_tag)
    cleanup()
    r
  })
  results <- do.call(rbind, lapply(out, `[[`, "result"))
  folds   <- do.call(rbind, lapply(out, `[[`, "folds"))
  write.csv(results, out_path(paste0(run_tag, "_results.csv")), row.names = FALSE)
  write.csv(folds,   out_path(paste0(run_tag, "_cv_folds.csv")), row.names = FALSE)
  extra <- if (exists("software_note", mode = "function")) software_note() else character(0)
  writeLines(c(extra, capture.output(sessionInfo())), out_path(paste0(run_tag, "_sessionInfo.txt")))
  cat("\n========== SUMMARY:", run_tag, "==========\n")
  print(results[, c("Scenario", "Features", "CV_AUC", "Test_AUC", "CI_Low", "CI_High",
                    "Sensitivity", "Specificity", "MCC", "F1", "Selected_config")],
        row.names = FALSE, digits = 3)
  cat("Saved:", out_path(paste0(run_tag, "_results.csv")), "and", out_path(paste0(run_tag, "_cv_folds.csv")), "\n")
  invisible(results)
}
