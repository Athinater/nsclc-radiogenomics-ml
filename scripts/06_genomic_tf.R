# ==============================================================================
# GENOMICS PIPELINE: GENES ONLY | TENSORFLOW | NESTED CV
# Gene selection: PCA loadings (PC1 top10 + PC2 top10) on TRAINING SET ONLY
# Zero leakage: test set never seen during gene selection
# Scenarios: S0, SA, SB, SC | Metrics: AUC, CI, Gap, MCC, F1, Sens, Spec
# ==============================================================================

# 1. LOAD LIBRARIES
library(keras3)
library(tensorflow)
library(caret)
library(pROC)
library(glmnet)
library(mRMRe)
library(DESeq2)
library(data.table)
source("scripts/common.R")
run_tag <- "genomic_tf"

# ---- Reproducibility ----
set.seed(123)
tf$random$set_seed(123)

# ==============================================================================
# 2. LOAD RAW COUNT DATA
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
y_train   <- as.numeric(trainResp) - 1

testResp <- clinical_data[[response_variable]][-trainIndex][
  match(test_rna_ids, toupper(trimws(rownames(clinical_data)[-trainIndex])))]
testResp <- droplevels(testResp)
y_test   <- as.numeric(testResp) - 1

cat("Train samples:", length(train_rna_ids), "\n")
cat("Test  samples:", length(test_rna_ids),  "\n")
cat("Train class dist:\n"); print(table(trainResp))

# ==============================================================================
# 6. GENE SELECTION: PCA LOADINGS ON TRAINING SET ONLY
# ==============================================================================
cat("\n=== Gene Selection: PCA on Training Set Only ===\n")

dds_train_gs <- DESeqDataSetFromMatrix(tbl_train, coldata_train, ~ Histology)
dds_train_gs <- estimateSizeFactors(dds_train_gs)
dds_train_gs <- estimateDispersions(dds_train_gs, quiet = TRUE)
vsd_train_gs <- vst(dds_train_gs, blind = FALSE)

pca_gs  <- prcomp(t(assay(vsd_train_gs)))
load_gs <- pca_gs$rotation

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
# ==============================================================================
vst_train_mat <- t(assay(vsd_train_gs)[selected_genes, , drop = FALSE])

dds_test_vst <- DESeqDataSetFromMatrix(tbl_test, coldata_test, ~ Histology)
dds_test_vst <- estimateSizeFactors(dds_test_vst)
dispersions(dds_test_vst) <- dispersions(dds_train_gs)[
  match(rownames(dds_test_vst), rownames(dds_train_gs))
]
vsd_test_vst <- varianceStabilizingTransformation(dds_test_vst, blind = FALSE)
vst_test_mat <- t(assay(vsd_test_vst)[selected_genes, , drop = FALSE])

# ==============================================================================
# 8. ALIGN & NORMALIZE
# ==============================================================================
tr_common <- intersect(rownames(vst_train_mat), train_rna_ids)
te_common <- intersect(rownames(vst_test_mat),  test_rna_ids)

X_train_raw <- vst_train_mat[tr_common, , drop = FALSE]
X_test_raw  <- vst_test_mat[te_common,  , drop = FALSE]

trainResp <- trainResp[match(tr_common, train_rna_ids)]; trainResp <- droplevels(trainResp)
testResp  <- testResp[match(te_common,  test_rna_ids)];  testResp  <- droplevels(testResp)
y_train   <- as.numeric(trainResp) - 1
y_test    <- as.numeric(testResp)  - 1

nzv <- nearZeroVar(as.data.frame(X_train_raw))
if (length(nzv) > 0) {
  X_train_raw <- X_train_raw[, -nzv, drop = FALSE]
  X_test_raw  <- X_test_raw[,  -nzv, drop = FALSE]
}

train_means <- colMeans(X_train_raw)
train_sds   <- apply(X_train_raw, 2, sd)
train_sds[train_sds == 0 | is.na(train_sds)] <- 1
X_train_norm <- scale(X_train_raw, center = train_means, scale = train_sds)
X_test_norm  <- scale(X_test_raw,  center = train_means, scale = train_sds)

cat("Final train:", nrow(X_train_norm), "x", ncol(X_train_norm), "\n")
cat("Final test: ", nrow(X_test_norm),  "x", ncol(X_test_norm),  "\n")

# ==============================================================================
# 9. SCENARIO FEATURE SETS
# ==============================================================================
feats_S0 <- colnames(X_train_norm)
cat("\nS0 features:", length(feats_S0), "\n")

cv_lasso <- cv.glmnet(as.matrix(X_train_norm), y_train, family = "binomial", alpha = 1, nfolds = 5)
co_SA    <- coef(cv_lasso, s = cv_lasso$lambda.min)
feats_SA <- rownames(co_SA)[co_SA[, 1] != 0]
feats_SA <- feats_SA[feats_SA != "(Intercept)"]
if (length(feats_SA) == 0)
  feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[, 1], decreasing = TRUE)[1:5])
cat("SA features (LASSO):", length(feats_SA), "\n")

n_mrmr <- max(length(feats_SA), 5)
d_mrmr <- mRMR.data(data = data.frame(X_train_norm, response = as.ordered(y_train)))
m_mrmr <- mRMR.classic(data = d_mrmr, target_indices = ncol(X_train_norm) + 1,
                        feature_count = n_mrmr)
feats_SB <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_mrmr@filters)]))
cat("SB features (LASSO+mRMR):", length(feats_SB), "\n")

pca_sc      <- prcomp(X_train_norm, scale. = FALSE)
pca_var     <- cumsum(pca_sc$sdev^2) / sum(pca_sc$sdev^2)
n_comp      <- min(which(pca_var >= 0.90)[1], ncol(X_train_norm) - 1)
pca_contrib <- apply(abs(pca_sc$rotation[, 1:n_comp, drop = FALSE]), 1, max)
pca_top     <- names(sort(pca_contrib, decreasing = TRUE)[1:n_comp])
feats_SC    <- unique(c(feats_SB, pca_top))
cat("SC features (LASSO+mRMR+PCA):", length(feats_SC), "\n")

# ==============================================================================
# 10. HELPER FUNCTIONS
# ==============================================================================
compute_class_weights <- function(y) {
  n     <- length(y)
  n_neg <- sum(y == 0); n_pos <- sum(y == 1)
  list("0" = n / (2 * n_neg), "1" = n / (2 * n_pos))
}

build_model <- function(input_dim, hidden_units, dropout_rate, l1_reg, l2_reg) {
  model <- keras_model_sequential() %>%
    layer_dense(units = hidden_units[1], activation = "relu",
                input_shape = input_dim,
                kernel_regularizer = regularizer_l1_l2(l1 = l1_reg, l2 = l2_reg)) %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(units = hidden_units[2], activation = "relu",
                kernel_regularizer = regularizer_l1_l2(l1 = l1_reg, l2 = l2_reg)) %>%
    layer_dropout(rate = dropout_rate) %>%
    layer_dense(units = 1, activation = "sigmoid")

  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 1e-3),
    loss      = "binary_crossentropy",
    metrics   = list(metric_auc(name = "auc"))
  )
  return(model)
}

# ==============================================================================
# 11. NESTED CV FUNCTION
# ==============================================================================
run_tf_genomics <- function(scenario_type, scenario_name, hidden_list) {
  cat(sprintf("\n=== Nested CV: %s ===\n", scenario_name))

  set.seed(123)
  tf$random$set_seed(123)
  folds   <- createFolds(y_train, k = 5, list = TRUE)
  cv_aucs <- numeric(5)
  oof_probs  <- numeric(0)   # out-of-fold predicted probabilities
  oof_labels <- numeric(0)   # corresponding true labels
  
  for (k in 1:5) {
    val_idx <- folds[[k]]
    tr_idx  <- setdiff(seq_along(y_train), val_idx)

    fold_tr_ids  <- tr_common[tr_idx]
    fold_val_ids <- tr_common[val_idx]
    fold_tr_resp <- trainResp[tr_idx]
    fold_val_resp <- trainResp[val_idx]

    # PCA gene selection inside fold
    fold_tbl_tr     <- tbl_train[, fold_tr_ids, drop = FALSE]
    fold_coldata_tr <- coldata_train[fold_tr_ids, , drop = FALSE]

    fold_dds <- tryCatch({
      d <- DESeqDataSetFromMatrix(fold_tbl_tr, fold_coldata_tr, ~ Histology)
      d <- estimateSizeFactors(d)
      d <- estimateDispersions(d, quiet = TRUE)
      d
    }, error = function(e) { cat(sprintf("  Fold %d DESeq2 error\n", k)); NULL })

    if (is.null(fold_dds)) { cv_aucs[k] <- NA; next }

    fold_vsd   <- vst(fold_dds, blind = FALSE)
    fold_pca   <- prcomp(t(assay(fold_vsd)))
    fold_load  <- fold_pca$rotation
    fold_genes <- unique(c(
      rownames(fold_load)[head(order(abs(fold_load[,"PC1"]), decreasing=TRUE), 10)],
      rownames(fold_load)[head(order(abs(fold_load[,"PC2"]), decreasing=TRUE), 10)]
    ))

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

    f_means <- colMeans(fold_vst_tr)
    f_sds   <- apply(fold_vst_tr, 2, sd); f_sds[f_sds == 0 | is.na(f_sds)] <- 1
    X_f_tr  <- scale(fold_vst_tr, center = f_means, scale = f_sds)
    X_f_va  <- scale(fold_vst_val, center = f_means, scale = f_sds)
    y_f_tr  <- as.numeric(fold_tr_resp) - 1
    y_f_va  <- as.numeric(fold_val_resp) - 1

    f_fold <- switch(scenario_type,
      "S0" = colnames(X_f_tr),
      "SA" = {
        lc <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
        co <- coef(lc, s=lc$lambda.min)
        fs <- rownames(co)[co[,1]!=0]; fs <- fs[fs!="(Intercept)"]
        if (length(fs)==0) colnames(X_f_tr) else fs
      },
      {
        lc <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
        co <- coef(lc, s=lc$lambda.min)
        fl <- rownames(co)[co[,1]!=0]; fl <- fl[fl!="(Intercept)"]
        if (length(fl)==0) fl <- colnames(X_f_tr)
        nm <- min(max(length(fl),5), ncol(X_f_tr))
        dm <- mRMR.data(data=data.frame(X_f_tr, response=as.ordered(y_f_tr)))
        mr <- mRMR.classic(data=dm, target_indices=ncol(X_f_tr)+1, feature_count=nm)
        fs <- unique(c(fl, colnames(X_f_tr)[unlist(mr@filters)]))
        if (scenario_type=="SC" && nrow(X_f_tr) > ncol(X_f_tr)) {
          pc   <- prcomp(X_f_tr, scale.=FALSE)
          pvar <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
          nc   <- min(which(pvar>=0.90)[1], ncol(X_f_tr)-1)
          pcon <- apply(abs(pc$rotation[,1:nc,drop=FALSE]),1,max)
          fs   <- unique(c(fs, names(sort(pcon,decreasing=TRUE)[1:nc])))
        }
        fs
      }
    )

    f_fold <- intersect(f_fold, intersect(colnames(X_f_tr), colnames(X_f_va)))
    if (length(f_fold) == 0) { cv_aucs[k] <- NA; next }

    X_cv_tr <- as.matrix(X_f_tr[, f_fold, drop=FALSE])
    X_cv_va <- as.matrix(X_f_va[, f_fold, drop=FALSE])

    m_cv <- build_model(ncol(X_cv_tr), hidden_list[[1]], 0.1, 0, 0)
    m_cv %>% fit(X_cv_tr, y_f_tr, epochs=40, batch_size=8,
                 class_weight=compute_class_weights(y_f_tr), verbose=0)
    fold_probs  <- as.vector(predict(m_cv, X_cv_va, verbose=0))
    cv_aucs[k]  <- as.numeric(roc(y_f_va, fold_probs, quiet=TRUE)$auc)
    oof_probs   <- c(oof_probs,  fold_probs)
    oof_labels  <- c(oof_labels, y_f_va)
    cat(sprintf("  Fold %d: %d genes -> %d feats | AUC=%.4f\n",
                k, length(fold_genes), length(f_fold), cv_aucs[k]))
  }

  valid_aucs  <- cv_aucs[!is.na(cv_aucs)]
  mean_cv_auc <- mean(valid_aucs)
  cat(sprintf("  Mean CV AUC (%d/5 folds): %.4f\n", length(valid_aucs), mean_cv_auc))
  oof_df <- data.frame(label = oof_labels, prob = oof_probs)
  write.csv(oof_df, out_path(paste0(run_tag, "_oof_preds_", scenario_type, ".csv")), row.names=FALSE)
  
  # Final model
  final_f <- intersect(
    switch(scenario_type, "S0"=feats_S0, "SA"=feats_SA, "SB"=feats_SB, "SC"=feats_SC),
    intersect(colnames(X_train_norm), colnames(X_test_norm))
  )

  X_tr_f <- as.matrix(X_train_norm[, final_f, drop=FALSE])
  X_ts_f <- as.matrix(X_test_norm[,  final_f, drop=FALSE])

  set.seed(123); tf$random$set_seed(123)
  val_grid_idx <- createDataPartition(y_train, p=0.20, list=FALSE)
  X_tr_grid <- X_tr_f[-val_grid_idx, ]; y_tr_grid <- y_train[-val_grid_idx]
  X_va_grid <- X_tr_f[ val_grid_idx, ]; y_va_grid <- y_train[ val_grid_idx]

  best_val_auc <- -Inf; best_model <- NULL

  for (cfg in hidden_list) {
    for (dr in c(0.1, 0.2)) {
      m_cand <- build_model(ncol(X_tr_grid), cfg, dr, 1e-4, 1e-4)
      m_cand %>% fit(X_tr_grid, y_tr_grid, epochs=40, batch_size=8,
                     class_weight=compute_class_weights(y_tr_grid), verbose=0)
      val_auc <- as.numeric(roc(y_va_grid,
                                 as.vector(predict(m_cand, X_va_grid, verbose=0)),
                                 quiet=TRUE)$auc)
      cat(sprintf("  hidden=%s dropout=%.1f -> val AUC=%.4f\n",
                  paste(cfg,collapse="x"), dr, val_auc))
      if (val_auc > best_val_auc) { best_val_auc <- val_auc; best_model <- m_cand }
    }
  }

  probs   <- as.vector(predict(best_model, X_ts_f, verbose=0))
  roc_obj <- roc(y_test, probs, quiet=TRUE)
  ci_auc  <- ci.auc(roc_obj)

  # Fixed decision threshold of 0.5 on P(SCC); Sens/Spec/MCC/F1 use ADC as the positive class.
  m <- adc_positive_metrics(as.integer(probs >= 0.5), y_test)   # 0 = ADC (positive), 1 = SCC
  write_preds(out_path(paste0(run_tag, "_preds_", scenario_type, ".csv")),
              y_test, probs, scenario_name)
  
  
  return(list(
    Name            = scenario_name,
    Features        = length(final_f),
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
# 12. RUN ALL SCENARIOS
# ==============================================================================
results_list <- list()
results_list[[1]] <- run_tf_genomics("S0", "GENES_S0_ALL",            list(c(16,16), c(10,10)))
results_list[[2]] <- run_tf_genomics("SA", "GENES_SA_LASSO",          list(c(8,8),   c(5,5)))
results_list[[3]] <- run_tf_genomics("SB", "GENES_SB_LASSO_mRMR",     list(c(16,16), c(10,10)))
results_list[[4]] <- run_tf_genomics("SC", "GENES_SC_LASSO_mRMR_PCA", list(c(16,16), c(10,10)))

# ==============================================================================
# 13. FINAL SUMMARY
# ==============================================================================
valid_results <- Filter(Negate(is.null), results_list)
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
  row <- final_results[i,]
  cat(sprintf("%-35s %5d %6s %7.4f %7.4f %8.4f %7.3f %7.3f %6.3f %6.3f\n",
              row$Name, row$Features, paste0(row$Valid_Folds,"/5"),
              row$CV_AUC, row$Test_AUC, row$Overfitting_Gap,
              row$Sensitivity, row$Specificity, row$MCC, row$F1))
}