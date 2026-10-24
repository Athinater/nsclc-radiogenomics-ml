# ==============================================================================
# GENOMIC PIPELINE: GENES ONLY | H2O | NESTED CV
# Gene selection: PCA loadings (PC1 top10 + PC2 top10) on the training set only;
# re-selected inside every cross-validation fold. Test set never used for selection.
# Scenarios: S0, SA, SB, SC | Metrics: AUC, CI, Gap, Sens, Spec, MCC, F1 (ADC = positive class)
# ==============================================================================

# 1. LOAD LIBRARIES
library(caret)
library(glmnet)
library(pROC)
library(mRMRe)
library(h2o)
library(DESeq2)
library(data.table)
source("scripts/common.R")
run_tag <- "genomic_h2o"


# ==============================================================================
# 2. LOAD RAW COUNT DATA (Entrez GeneID -> Symbol mapping via org.Hs.eg.db)
# ==============================================================================
tbl <- load_counts(counts_file)
# ==============================================================================
# 3. LOAD CLINICAL METADATA (DESeq2 colData)
# ==============================================================================
coldata_full <- read.csv(
  clinical_meta_file,
  header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1
)
rownames(coldata_full) <- toupper(trimws(rownames(coldata_full)))
coldata_full <- coldata_full[coldata_full$Histology %in% c("Adenocarcinoma", "Squamous cell carcinoma"), ]

common_samples <- intersect(rownames(coldata_full), colnames(tbl))
coldata_full   <- coldata_full[common_samples, , drop = FALSE]
tbl_filtered   <- tbl[, common_samples, drop = FALSE]

# ==============================================================================
# 4. LOAD CLINICAL DATA (for the split)
# ==============================================================================
clinical_data <- read.csv(
  clinical_rg_file,
  sep = ";", header = TRUE, row.names = 1
)
response_variable <- "Histology"
clinical_data[[response_variable]] <- as.factor(clinical_data[[response_variable]])
levels(clinical_data[[response_variable]]) <- make.names(levels(clinical_data[[response_variable]]))

# ==============================================================================
# 5. STRATIFIED SPLIT (70/30)
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(clinical_data[[response_variable]], p = 0.7, list = FALSE)

train_ids <- toupper(trimws(rownames(clinical_data)[trainIndex]))
test_ids  <- toupper(trimws(rownames(clinical_data)[-trainIndex]))

train_rna_ids <- intersect(train_ids, colnames(tbl_filtered))
test_rna_ids  <- intersect(test_ids,  colnames(tbl_filtered))

tbl_train     <- tbl_filtered[, train_rna_ids, drop = FALSE]
tbl_test      <- tbl_filtered[, test_rna_ids,  drop = FALSE]
coldata_train <- coldata_full[train_rna_ids, , drop = FALSE]
coldata_test  <- coldata_full[test_rna_ids,  , drop = FALSE]

trainResp <- clinical_data[[response_variable]][trainIndex][
  match(train_rna_ids, toupper(trimws(rownames(clinical_data)[trainIndex])))]
trainResp <- droplevels(trainResp)

testResp <- clinical_data[[response_variable]][-trainIndex][
  match(test_rna_ids, toupper(trimws(rownames(clinical_data)[-trainIndex])))]
testResp <- droplevels(testResp)

cat("Train samples:", length(train_rna_ids), "\n")
cat("Test  samples:", length(test_rna_ids),  "\n")
cat("Train class dist:\n"); print(table(trainResp))
cat("Test  class dist:\n"); print(table(testResp))

# ==============================================================================
# 6. GENE SELECTION: PCA LOADINGS ON TRAINING SET ONLY
# ==============================================================================
cat("\n=== Gene Selection: PCA on Training Set Only ===\n")

dds_train_gs <- DESeqDataSetFromMatrix(tbl_train, coldata_train, ~ Histology)
dds_train_gs <- estimateSizeFactors(dds_train_gs)
dds_train_gs <- estimateDispersions(dds_train_gs, quiet = TRUE)
vsd_train_gs <- vst(dds_train_gs, blind = FALSE)

pca_gs   <- prcomp(t(assay(vsd_train_gs)))
load_gs  <- pca_gs$rotation

top_pc1 <- rownames(load_gs)[head(order(abs(load_gs[, "PC1"]), decreasing = TRUE), 10)]
top_pc2 <- rownames(load_gs)[head(order(abs(load_gs[, "PC2"]), decreasing = TRUE), 10)]
selected_genes <- unique(c(top_pc1, top_pc2))
check_selected_genes(selected_genes)

cat("PC1 top 10:", paste(top_pc1, collapse = ", "), "\n")
cat("PC2 top 10:", paste(top_pc2, collapse = ", "), "\n")
cat("Total selected genes:", length(selected_genes), "\n")

var_exp <- summary(pca_gs)$importance[2, 1:2]
cat(sprintf("Variance explained - PC1: %.1f%%  PC2: %.1f%%\n",
            var_exp[1]*100, var_exp[2]*100))

# ==============================================================================
# 7. VST NORMALIZATION - TRAIN & TEST
#    VST fit on training set, applied to test set via size factors only
# ==============================================================================

# Training VST (already done above - reuse)
vst_train_mat <- t(assay(vsd_train_gs)[selected_genes, , drop = FALSE])  # samples x genes

# Test VST: test-specific size factors + train dispersions (Manual VST)
dds_test_vst <- DESeqDataSetFromMatrix(tbl_test, coldata_test, ~ Histology)
dds_test_vst <- estimateSizeFactors(dds_test_vst)
dispersions(dds_test_vst) <- dispersions(dds_train_gs)[
  match(rownames(dds_test_vst), rownames(dds_train_gs))
]
vsd_test_vst <- varianceStabilizingTransformation(dds_test_vst, blind = FALSE)
vst_test_mat <- t(assay(vsd_test_vst)[selected_genes, , drop = FALSE])   # samples x genes

cat("\nVST train matrix:", nrow(vst_train_mat), "samples x", ncol(vst_train_mat), "genes\n")
cat("VST test  matrix:", nrow(vst_test_mat),  "samples x", ncol(vst_test_mat),  "genes\n")

# ==============================================================================
# 8. ALIGN & NORMALIZE
# ==============================================================================
tr_common <- intersect(rownames(vst_train_mat), train_rna_ids)
te_common <- intersect(rownames(vst_test_mat),  test_rna_ids)

X_train_raw <- vst_train_mat[tr_common, , drop = FALSE]
X_test_raw  <- vst_test_mat[te_common,  , drop = FALSE]

trainResp <- trainResp[match(tr_common, train_rna_ids)]
testResp  <- testResp[match(te_common,  test_rna_ids)]
trainResp <- droplevels(trainResp)
testResp  <- droplevels(testResp)

y_train <- as.numeric(trainResp) - 1

# Remove near-zero variance (train statistics only)
nzv <- nearZeroVar(as.data.frame(X_train_raw))
if (length(nzv) > 0) {
  cat("Removing", length(nzv), "near-zero variance genes\n")
  X_train_raw <- X_train_raw[, -nzv, drop = FALSE]
  X_test_raw  <- X_test_raw[,  -nzv, drop = FALSE]
}

# Scale (train mean/sd only)
train_means <- colMeans(X_train_raw)
train_sds   <- apply(X_train_raw, 2, sd)
train_sds[train_sds == 0 | is.na(train_sds)] <- 1

X_train_norm <- scale(X_train_raw, center = train_means, scale = train_sds)
X_test_norm  <- scale(X_test_raw,  center = train_means, scale = train_sds)

cat("Final train:", nrow(X_train_norm), "samples x", ncol(X_train_norm), "genes\n")
cat("Final test: ", nrow(X_test_norm),  "samples x", ncol(X_test_norm),  "genes\n")

# ==============================================================================
# 9. SCENARIO FEATURE SETS (computed on full training set)
# ==============================================================================

# S0: all selected genes (from PCA)
feats_S0 <- colnames(X_train_norm)
cat("\nS0 features (all PCA genes):", length(feats_S0), "\n")

# SA: LASSO on S0 genes
cv_lasso <- cv.glmnet(as.matrix(X_train_norm), y_train, family = "binomial", alpha = 1, nfolds = 5)
co_SA    <- coef(cv_lasso, s = cv_lasso$lambda.min)
feats_SA <- rownames(co_SA)[co_SA[, 1] != 0]
feats_SA <- feats_SA[feats_SA != "(Intercept)"]
if (length(feats_SA) == 0)
  feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[, 1], decreasing = TRUE)[1:10])
cat("SA features (LASSO, lambda.min =", round(cv_lasso$lambda.min, 5), "):", length(feats_SA), "\n")

# SB: LASSO + mRMR
n_mrmr <- max(length(feats_SA), 5)
d_mrmr <- mRMR.data(data = data.frame(X_train_norm, response = as.ordered(y_train)))
m_mrmr <- mRMR.classic(data = d_mrmr,
                       target_indices = ncol(X_train_norm) + 1,
                       feature_count  = n_mrmr)
feats_SB <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_mrmr@filters)]))
cat("SB features (LASSO + mRMR):", length(feats_SB), "\n")

# SC: LASSO + mRMR + PCA components
pca_sc   <- prcomp(X_train_norm, scale. = FALSE)
pca_var  <- cumsum(pca_sc$sdev^2) / sum(pca_sc$sdev^2)
n_comp   <- min(which(pca_var >= 0.90)[1], ncol(X_train_norm) - 1)
cat("SC PCA components (>=90% variance):", n_comp,
    "| Variance explained:", round(pca_var[n_comp] * 100, 1), "%\n")
pca_contrib <- apply(abs(pca_sc$rotation[, 1:n_comp, drop = FALSE]), 1, max)
pca_top     <- names(sort(pca_contrib, decreasing = TRUE)[1:n_comp])
feats_SC    <- unique(c(feats_SB, pca_top))
cat("SC features (LASSO + mRMR + PCA):", length(feats_SC), "\n")

# ==============================================================================
# 10. NESTED CV WITH H2O
# ==============================================================================
h2o.init(strict_version_check = FALSE)

run_h2o_genomics <- function(scenario_type, scenario_name, hidden_list) {
  cat(sprintf("\n=== Nested CV: %s ===\n", scenario_name))
  
  set.seed(123)
  folds_idx <- createFolds(trainResp, k = 5, list = TRUE)
  cv_aucs   <- numeric(5)
  
  for (k in 1:5) {
    val_idx <- folds_idx[[k]]
    tr_idx  <- setdiff(seq_along(trainResp), val_idx)
    
    # Fold train/val IDs for PCA-based gene selection
    fold_tr_ids  <- tr_common[tr_idx]
    fold_val_ids <- tr_common[val_idx]
    
    fold_tr_resp  <- trainResp[tr_idx]
    fold_val_resp <- trainResp[val_idx]
    
    # ---- PCA gene selection INSIDE fold (on 4/5 training data) ----
    fold_tbl_tr     <- tbl_train[, fold_tr_ids, drop = FALSE]
    fold_coldata_tr <- coldata_train[fold_tr_ids, , drop = FALSE]
    
    fold_dds <- tryCatch({
      d <- DESeqDataSetFromMatrix(fold_tbl_tr, fold_coldata_tr, ~ Histology)
      d <- estimateSizeFactors(d)
      d <- estimateDispersions(d, quiet = TRUE)
      d
    }, error = function(e) {
      cat(sprintf("  Fold %d DESeq2 error: %s\n", k, e$message)); NULL
    })
    
    if (is.null(fold_dds)) { cv_aucs[k] <- NA; next }
    
    fold_vsd      <- vst(fold_dds, blind = FALSE)
    fold_pca      <- prcomp(t(assay(fold_vsd)))
    fold_load     <- fold_pca$rotation
    fold_top_pc1  <- rownames(fold_load)[head(order(abs(fold_load[,"PC1"]), decreasing=TRUE), 10)]
    fold_top_pc2  <- rownames(fold_load)[head(order(abs(fold_load[,"PC2"]), decreasing=TRUE), 10)]
    fold_genes    <- unique(c(fold_top_pc1, fold_top_pc2))
    
    # VST matrices for this fold
    fold_vst_tr <- t(assay(fold_vsd)[fold_genes, , drop = FALSE])
    
    fold_tbl_val     <- tbl_train[, fold_val_ids, drop = FALSE]
    fold_coldata_val <- coldata_train[fold_val_ids, , drop = FALSE]
    fold_dds_val     <- DESeqDataSetFromMatrix(fold_tbl_val, fold_coldata_val, ~ Histology)
    fold_dds_val     <- estimateSizeFactors(fold_dds_val)
    dispersions(fold_dds_val) <- dispersions(fold_dds)[
      match(rownames(fold_dds_val), rownames(fold_dds))
    ]
    fold_vsd_val     <- varianceStabilizingTransformation(fold_dds_val, blind = FALSE)
    
    avail_genes  <- intersect(fold_genes, rownames(assay(fold_vsd_val)))
    fold_vst_val <- t(assay(fold_vsd_val)[avail_genes, , drop = FALSE])
    fold_vst_tr  <- fold_vst_tr[, avail_genes, drop = FALSE]
    
    # Normalize (fold train stats)
    f_means <- colMeans(fold_vst_tr)
    f_sds   <- apply(fold_vst_tr, 2, sd)
    f_sds[f_sds == 0 | is.na(f_sds)] <- 1
    X_f_tr <- scale(fold_vst_tr, center = f_means, scale = f_sds)
    X_f_va <- scale(fold_vst_val, center = f_means, scale = f_sds)
    
    y_f_tr <- as.numeric(fold_tr_resp) - 1
    
    # Feature selection per scenario
    f_fold <- switch(scenario_type,
                     "S0" = colnames(X_f_tr),
                     "SA" = {
                       lc <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
                       co <- coef(lc, s=lc$lambda.min)
                       fs <- rownames(co)[co[,1] != 0]; fs <- fs[fs != "(Intercept)"]
                       if (length(fs) == 0) colnames(X_f_tr) else fs
                     },
                     {
                       lc <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
                       co <- coef(lc, s=lc$lambda.min)
                       fl <- rownames(co)[co[,1] != 0]; fl <- fl[fl != "(Intercept)"]
                       if (length(fl) == 0) fl <- colnames(X_f_tr)
                       nm <- min(max(length(fl), 5), ncol(X_f_tr))
                       dm <- mRMR.data(data=data.frame(X_f_tr, response=as.ordered(y_f_tr)))
                       mr <- mRMR.classic(data=dm, target_indices=ncol(X_f_tr)+1, feature_count=nm)
                       fs <- unique(c(fl, colnames(X_f_tr)[unlist(mr@filters)]))
                       if (scenario_type == "SC" && nrow(X_f_tr) > ncol(X_f_tr)) {
                         pc   <- prcomp(X_f_tr, scale.=FALSE)
                         pvar <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
                         nc   <- min(which(pvar >= 0.90)[1], ncol(X_f_tr)-1)
                         pcon <- apply(abs(pc$rotation[,1:nc,drop=FALSE]),1,max)
                         fs   <- unique(c(fs, names(sort(pcon,decreasing=TRUE)[1:nc])))
                       }
                       fs
                     }
    )
    
    f_fold <- intersect(f_fold, intersect(colnames(X_f_tr), colnames(X_f_va)))
    if (length(f_fold) == 0) { cv_aucs[k] <- NA; next }
    
    cat(sprintf("  Fold %d: %d genes -> %d features\n", k, length(fold_genes), length(f_fold)))
    
    tr_h2o_f <- as.h2o(cbind(as.data.frame(X_f_tr[, f_fold, drop=FALSE]),
                             Histology = fold_tr_resp))
    va_h2o_f <- as.h2o(cbind(as.data.frame(X_f_va[, f_fold, drop=FALSE]),
                             Histology = fold_val_resp))
    
    m_f <- tryCatch(
      h2o.deeplearning(
        x = f_fold, y = "Histology",
        training_frame   = tr_h2o_f,
        validation_frame = va_h2o_f,
        epochs           = 40,
        hidden           = hidden_list[[1]],
        balance_classes  = FALSE,
        seed             = 123
      ), error = function(e) NULL
    )
    
    cv_aucs[k] <- if (!is.null(m_f)) { a <- h2o.auc(m_f, valid=TRUE); h2o.rm(m_f); a } else NA
    cat(sprintf("  Fold %d AUC: %s\n", k,
                if(is.na(cv_aucs[k])) "SKIPPED" else round(cv_aucs[k], 4)))
  }
  
  valid_aucs  <- cv_aucs[!is.na(cv_aucs)]
  mean_cv_auc <- if (length(valid_aucs) > 0) mean(valid_aucs) else NA
  cat(sprintf("  Mean CV AUC (%d/5 valid folds): %s\n",
              length(valid_aucs),
              if(is.na(mean_cv_auc)) "NA" else round(mean_cv_auc, 4)))
  
  # ===== FINAL MODEL ON FULL TRAINING SET =====
  final_features <- switch(scenario_type,
                           "S0" = feats_S0,
                           "SA" = feats_SA,
                           "SB" = feats_SB,
                           "SC" = feats_SC)
  final_features <- intersect(make.names(final_features),
                              intersect(make.names(colnames(X_train_norm)),
                                        make.names(colnames(X_test_norm))))
  
  tr_f_df <- as.data.frame(X_train_norm[, final_features, drop = FALSE])
  ts_f_df <- as.data.frame(X_test_norm[,  final_features, drop = FALSE])
  
  set.seed(123)
  val_grid_idx <- createDataPartition(trainResp, p = 0.20, list = FALSE)
  tr_h2o <- as.h2o(cbind(tr_f_df[-val_grid_idx, ], Histology = trainResp[-val_grid_idx]))
  va_h2o <- as.h2o(cbind(tr_f_df[ val_grid_idx, ], Histology = trainResp[ val_grid_idx]))
  ts_h2o <- as.h2o(cbind(ts_f_df, Histology = testResp))
  
  best_m <- NULL; best_val_auc <- -Inf
  
  for (hid in hidden_list) {
    for (idr in c(0.0, 0.1)) {
      tryCatch({
        cand <- h2o.deeplearning(
          x = final_features, y = "Histology",
          training_frame      = tr_h2o,
          validation_frame    = va_h2o,
          epochs              = 40,
          activation          = "RectifierWithDropout",
          hidden              = hid,
          input_dropout_ratio = idr,
          balance_classes     = FALSE,
          reproducible        = TRUE,
          seed                = 123
        )
        v <- h2o.auc(cand, valid = TRUE)
        cat(sprintf("  hidden=%s dropout=%.1f -> val AUC=%.4f\n",
                    paste(hid, collapse="x"), idr, v))
        if (v > best_val_auc) {
          if (!is.null(best_m)) h2o.rm(best_m)
          best_val_auc <- v; best_m <- cand
        } else h2o.rm(cand)
      }, error = function(e) cat(sprintf("  SKIP %s/%.1f: %s\n",
                                         paste(hid,collapse="x"), idr, e$message)))
    }
  }
  
  if (is.null(best_m)) { h2o.removeAll(); return(NULL) }
  
  pred   <- as.data.frame(h2o.predict(best_m, newdata = ts_h2o))
  resp_l <- levels(testResp)
  colnames(pred) <- c("predict", resp_l[1], resp_l[2])
  
  roc_obj <- roc(testResp, pred[[resp_l[2]]], quiet = TRUE)
  ci_auc  <- ci.auc(roc_obj)

  # Metrics with ADC as the positive class (as reported in the manuscript).
  # Class labels come from H2O's default decision threshold.
  truth01 <- as.numeric(testResp) - 1          # 0 = ADC, 1 = SCC
  m <- adc_positive_metrics(as.integer(pred$predict == resp_l[2]), truth01)
  write_preds(out_path(paste0(run_tag, "_preds_", scenario_type, ".csv")),
              truth01, pred[[resp_l[2]]], scenario_name)
  
  h2o.removeAll()
  
  return(list(
    Name            = scenario_name,
    Features        = length(final_features),
    Valid_Folds     = length(valid_aucs),
    CV_AUC          = mean_cv_auc,
    Test_AUC        = as.numeric(roc_obj$auc),
    Overfitting_Gap = mean_cv_auc - as.numeric(roc_obj$auc),
    CI_Low          = ci_auc[1],
    CI_High         = ci_auc[3],
    Sensitivity     = m$Sensitivity,
    Specificity     = m$Specificity,
    MCC             = m$MCC,
    F1              = m$F1,
    TP = m$TP, FN = m$FN, TN = m$TN, FP = m$FP
  ))
}

# ==============================================================================
# 11. RUN ALL SCENARIOS
# ==============================================================================
results_list <- list()
results_list[[1]] <- run_h2o_genomics("S0", "GENES_S0_ALL",            list(c(25,25), c(20,20)))
results_list[[2]] <- run_h2o_genomics("SA", "GENES_SA_LASSO",          list(c(28,28),   c(25,25)))
results_list[[3]] <- run_h2o_genomics("SB", "GENES_SB_LASSO_mRMR",     list(c(20,20), c(18,18)))
results_list[[4]] <- run_h2o_genomics("SC", "GENES_SC_LASSO_mRMR_PCA", list(c(26,26), c(18,18)))

# ==============================================================================
# 12. FINAL SUMMARY
# ==============================================================================
valid_results <- Filter(Negate(is.null), results_list)

if (length(valid_results) > 0) {
  final_results <- do.call(rbind, lapply(valid_results, as.data.frame))
  write.csv(final_results, out_path(paste0(run_tag, "_results.csv")), row.names = FALSE)
  
  cat("\n\n========== FINAL SUMMARY ==========\n")
  cat(sprintf("Gene selection: PCA loadings (PC1 top10 + PC2 top10) on training set only\n"))
  cat(sprintf("Selected genes (%d): %s\n\n", length(selected_genes),
              paste(selected_genes, collapse=", ")))
  
  cat(sprintf("%-35s %5s %6s %7s %7s %8s %7s %7s %6s %6s\n",
              "Scenario", "Feats", "Folds", "CV-AUC", "T-AUC", "Gap",
              "Sens", "Spec", "MCC", "F1"))
  cat(strrep("-", 105), "\n")
  
  for (i in 1:nrow(final_results)) {
    row <- final_results[i, ]
    cat(sprintf("%-35s %5d %6s %7.4f %7.4f %8.4f %7.3f %7.3f %6.3f %6.3f\n",
                row$Name, row$Features,
                paste0(row$Valid_Folds, "/5"),
                row$CV_AUC, row$Test_AUC, row$Overfitting_Gap,
                row$Sensitivity, row$Specificity, row$MCC, row$F1))
  }
} else {
  cat("\nWARNING: All scenarios returned NULL.\n")
}

# ==============================================================================
# 13. GENE SELECTION REPORT
# ==============================================================================
cat("\n=== GENE SELECTION REPORT ===\n")
cat("Method: estimateSizeFactors + estimateDispersions (train) -> vst(blind=FALSE) -> prcomp -> PC1/PC2 top10\n")
cat("VST test: estimateSizeFactors (test) + dispersions(train) -> varianceStabilizingTransformation\n")
cat("PC1 top 10:\n"); print(top_pc1)
cat("PC2 top 10:\n"); print(top_pc2)
cat(sprintf("Total unique genes: %d\n", length(selected_genes)))
cat(sprintf("Variance explained: PC1=%.1f%%, PC2=%.1f%%\n",
            var_exp[1]*100, var_exp[2]*100))

cat(sprintf("\nAll outputs written to: %s\n", output_dir))
