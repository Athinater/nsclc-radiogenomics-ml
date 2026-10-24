# ==============================================================================
# RADIOGENOMIC PIPELINE: NESTED CV | TENSORFLOW
# Gene selection: PCA loadings (PC1 top10 + PC2 top10) on the training set only;
# re-selected inside every cross-validation fold. Test set never used for selection.
# Scenarios: S0, SA, SB, SC | Metrics: AUC, CI, Gap, Sens, Spec, MCC, F1 (ADC = positive class)
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
run_tag <- "radiogenomic_tf"

# ---- Reproducibility ----
set.seed(123)
tf$random$set_seed(123)


# ==============================================================================
# 2. LOAD RAW COUNT DATA (for DESeq2)
# ==============================================================================
tbl <- load_counts(counts_file)

# ==============================================================================
# 3. LOAD CLINICAL + RADIOMIC DATA
# ==============================================================================
clinical_data <- read.csv(
  clinical_rg_file,
  sep = ";", header = TRUE, row.names = 1
)

response_variable <- "Histology"
clinical_data[[response_variable]] <- as.factor(clinical_data[[response_variable]])
levels(clinical_data[[response_variable]]) <- make.names(levels(clinical_data[[response_variable]]))

imaging_features <- colnames(clinical_data)[grepl(
  "MORPHOLOGICAL_|INTENSITY.BASED_|INTENSITY.HISTOGRAM_|GLCM_|GLRLM_|GLSZM_|NGTDM_",
  colnames(clinical_data)
)]

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
# 4. STRATIFIED SPLIT (70/30)
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(clinical_data[[response_variable]], p = 0.7, list = FALSE)

train_ids <- toupper(trimws(rownames(clinical_data)[trainIndex]))
test_ids  <- toupper(trimws(rownames(clinical_data)[-trainIndex]))

# ==============================================================================
# 5. GENE SELECTION: PCA LOADINGS ON TRAINING SET ONLY
# ==============================================================================
cat("\n=== Gene Selection: PCA on Training Set Only ===\n")

train_rna_ids <- intersect(train_ids, colnames(tbl_filtered))
test_rna_ids  <- intersect(test_ids,  colnames(tbl_filtered))

coldata_train <- coldata_full[train_rna_ids, , drop = FALSE]
tbl_train     <- tbl_filtered[, train_rna_ids, drop = FALSE]
coldata_test  <- coldata_full[test_rna_ids,  , drop = FALSE]
tbl_test      <- tbl_filtered[, test_rna_ids,  drop = FALSE]

# DESeq2 + VST on training set only
dds_train <- DESeqDataSetFromMatrix(tbl_train, coldata_train, ~ Histology)
dds_train <- estimateSizeFactors(dds_train)
dds_train <- estimateDispersions(dds_train, quiet = TRUE)
vsd_train <- vst(dds_train, blind = FALSE)

# PCA on VST training matrix -> top loadings
pca_gs  <- prcomp(t(assay(vsd_train)))
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

# VST matrices (train: from dds_train; test: Manual VST with train dispersions)
vst_matrix_tr <- t(assay(vsd_train)[selected_genes, , drop = FALSE])

dds_test_vst <- DESeqDataSetFromMatrix(tbl_test, coldata_test, ~ Histology)
dds_test_vst <- estimateSizeFactors(dds_test_vst)
dispersions(dds_test_vst) <- dispersions(dds_train)[
  match(rownames(dds_test_vst), rownames(dds_train))
]
vsd_test      <- varianceStabilizingTransformation(dds_test_vst, blind = FALSE)
vst_matrix_te <- t(assay(vsd_test)[selected_genes, , drop = FALSE])

cat("Train RNA samples:", nrow(vst_matrix_tr), "\n")
cat("Test  RNA samples:", nrow(vst_matrix_te), "\n")

# ==============================================================================
# 6. COMBINE RADIOMIC + GENOMIC FEATURES
# ==============================================================================
trainData_radio <- clinical_data[trainIndex, imaging_features]
testData_radio  <- clinical_data[-trainIndex, imaging_features]

# Train alignment
tr_radio_ids            <- toupper(trimws(rownames(trainData_radio)))
tr_common               <- intersect(tr_radio_ids, rownames(vst_matrix_tr))
trainData_radio_aligned <- trainData_radio[match(tr_common, tr_radio_ids), ]
vst_tr_aligned          <- vst_matrix_tr[tr_common, ]
trainData_raw           <- cbind(trainData_radio_aligned, vst_tr_aligned)
trainResp               <- clinical_data[[response_variable]][trainIndex][match(tr_common, tr_radio_ids)]
trainResp               <- droplevels(trainResp)
y_train                 <- as.numeric(trainResp) - 1

# Test alignment
te_radio_ids           <- toupper(trimws(rownames(testData_radio)))
te_common              <- intersect(te_radio_ids, rownames(vst_matrix_te))
testData_radio_aligned <- testData_radio[match(te_common, te_radio_ids), ]
vst_te_aligned         <- vst_matrix_te[te_common, ]
testData_raw           <- cbind(testData_radio_aligned, vst_te_aligned)
testResp               <- clinical_data[[response_variable]][-trainIndex][match(te_common, te_radio_ids)]
testResp               <- droplevels(testResp)
y_test                 <- as.numeric(testResp) - 1

cat("Final train samples (radio+genomic):", nrow(trainData_raw), "\n")
cat("Final test  samples (radio+genomic):", nrow(testData_raw),  "\n")

# ==============================================================================
# 7. PREPROCESSING (NZV & NORMALIZATION)
# ==============================================================================
nzv <- nearZeroVar(trainData_raw)
if(length(nzv) > 0) {
  trainData_raw <- trainData_raw[, -nzv]
  testData_raw  <- testData_raw[,  -nzv]
}

train_means <- apply(trainData_raw, 2, mean)
train_sds   <- apply(trainData_raw, 2, sd)
train_sds[train_sds == 0 | is.na(train_sds)] <- 1
X_train_norm <- scale(trainData_raw, center = train_means, scale = train_sds)
X_test_norm  <- scale(testData_raw,  center = train_means, scale = train_sds)

# ==============================================================================
# 8. PRE-COMPUTE FEATURES FOR FINAL MODEL (FULL TRAINING SET)
# ==============================================================================

# Helper: top-N radiomic features by absolute correlation with the class
get_top_radio <- function(X, y, radio_feats, n = 5) {
  radio_avail <- intersect(radio_feats, colnames(X))
  if(length(radio_avail) == 0) return(character(0))
  cors <- abs(cor(X[, radio_avail, drop = FALSE], y)[, 1])
  names(sort(cors, decreasing = TRUE))[1:min(n, length(radio_avail))]
}

radio_clean    <- make.names(imaging_features)
radio_in_train <- intersect(radio_clean, colnames(X_train_norm))

# S0: Radiomic + PCA genomic genes
genes_clean         <- make.names(selected_genes)
feats_S0_candidates <- c(radio_clean, genes_clean)
feats_S0            <- intersect(feats_S0_candidates, colnames(X_train_norm))
cat("S0 features (radiomic + PCA genomic):", length(feats_S0), "\n")

# SA: LASSO + guaranteed top-5 radiomic
cv_l     <- cv.glmnet(as.matrix(X_train_norm), y_train, family = "binomial", alpha = 1, nfolds = 5)
co_SA    <- coef(cv_l, s = cv_l$lambda.min)
feats_SA <- rownames(co_SA)[co_SA[,1] != 0]
feats_SA <- feats_SA[feats_SA != "(Intercept)"]
if(length(feats_SA) == 0)
  feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[,1], decreasing=TRUE)[1:5])
top5_radio_ext <- get_top_radio(X_train_norm, y_train, radio_in_train, n = 5)
feats_SA <- unique(c(feats_SA, top5_radio_ext))
cat("SA features (LASSO + top5 radiomic guaranteed):", length(feats_SA), "\n")

# SB: LASSO + mRMR + guaranteed top-5 radiomic
n_mrmr_ext <- max(length(feats_SA), 5)
d_mrmr_ext <- mRMR.data(data = data.frame(X_train_norm, response = as.ordered(y_train)))
m_res_ext  <- mRMR.classic(data = d_mrmr_ext,
                           target_indices = ncol(X_train_norm) + 1,
                           feature_count  = n_mrmr_ext)
feats_SB <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_res_ext@filters)]))
cat("SB features (LASSO + mRMR + top5 radiomic guaranteed):", length(feats_SB), "\n")

# SC: LASSO + mRMR + PCA
pca_ext  <- prcomp(X_train_norm, scale. = FALSE)
pca_var  <- cumsum(pca_ext$sdev^2) / sum(pca_ext$sdev^2)
n_comp   <- min(which(pca_var >= 0.90)[1], 10)
cat("SC PCA components (>=90% variance):", n_comp,
    "| Variance explained:", round(pca_var[n_comp] * 100, 1), "%\n")

pca_contrib <- apply(abs(pca_ext$rotation[, 1:n_comp, drop = FALSE]), 1, max)
pca_top_ext <- names(sort(pca_contrib, decreasing = TRUE)[1:n_comp])
feats_SC    <- unique(c(feats_SB, pca_top_ext))
cat("SC features (LASSO + mRMR + PCA):", length(feats_SC), "\n")

# ==============================================================================
# 9. HELPER FUNCTIONS
# ==============================================================================

compute_class_weights <- function(y) {
  counts <- table(y)
  n      <- length(y)
  list(`0` = as.numeric((1 / counts["0"]) * (n / 2)),
       `1` = as.numeric((1 / counts["1"]) * (n / 2)))
}

build_model <- function(n_features, u, dr, l1, l2) {
  model <- keras_model_sequential() %>%
    layer_dense(units = u[1], activation = "relu",
                input_shape = c(n_features),
                kernel_regularizer = regularizer_l1_l2(l1 = l1, l2 = l2)) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dr) %>%
    layer_dense(units = u[2], activation = "relu",
                kernel_regularizer = regularizer_l1_l2(l1 = l1, l2 = l2)) %>%
    layer_dropout(rate = dr) %>%
    layer_dense(units = 1, activation = "sigmoid")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss      = "binary_crossentropy",
    metrics   = list(metric_auc(name = "auc"))
  )
  return(model)
}

# ==============================================================================
# 10. NESTED CV FUNCTION (NO LEAKAGE)
# ==============================================================================

# Context for gene re-selection inside every cross-validation fold (see common.R)
fold_ctx <- list(ids = tr_common, tbl_train = tbl_train, coldata_train = coldata_train,
                 radio_cols = radio_in_train)

run_tf_scenario_unbiased <- function(scenario_type, scenario_name, hidden_list) {
  cat(sprintf("\nRUNNING TF: %s\n", scenario_name))
  
  set.seed(123)
  folds   <- createFolds(y_train, k = 5, list = TRUE)
  cv_aucs <- numeric(5)
  
  for(k in 1:5) {
    val_idx <- folds[[k]]
    tr_idx  <- setdiff(seq_along(y_train), val_idx)
    fm <- build_fold_matrices(X_train_norm, tr_idx, val_idx, fold_ctx)   # genes re-selected within the fold
    X_f_tr <- fm$X_f_tr;  y_f_tr <- y_train[tr_idx]
    X_f_va <- fm$X_f_va;  y_f_va <- y_train[val_idx]
    
    f_fold <- NULL
    
    if(scenario_type == "S0") {
      f_fold <- colnames(X_f_tr)
    }
    
    if(scenario_type == "SA") {
      l_cv   <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family = "binomial", alpha = 1, nfolds = 5)
      co     <- coef(l_cv, s = l_cv$lambda.min)
      f_fold <- rownames(co)[co[,1] != 0]
      f_fold <- f_fold[f_fold != "(Intercept)"]
      if(length(f_fold) == 0)
        f_fold <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
      top5_radio_fold <- get_top_radio(X_f_tr, y_f_tr, radio_in_train, n = 5)
      f_fold <- unique(c(f_fold, top5_radio_fold))
    }
    
    if(scenario_type %in% c("SB", "SC")) {
      l_cv    <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family = "binomial", alpha = 1, nfolds = 5)
      co      <- coef(l_cv, s = l_cv$lambda.min)
      f_lasso <- rownames(co)[co[,1] != 0]
      f_lasso <- f_lasso[f_lasso != "(Intercept)"]
      if(length(f_lasso) == 0)
        f_lasso <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
      top5_radio_fold <- get_top_radio(X_f_tr, y_f_tr, radio_in_train, n = 5)
      f_lasso <- unique(c(f_lasso, top5_radio_fold))
      
      n_mrmr_fold <- max(length(f_lasso), 5)
      d_m    <- mRMR.data(data = data.frame(X_f_tr, response = as.ordered(y_f_tr)))
      m_r    <- mRMR.classic(data = d_m,
                             target_indices = ncol(X_f_tr) + 1,
                             feature_count  = n_mrmr_fold)
      f_fold <- unique(c(f_lasso, colnames(X_f_tr)[unlist(m_r@filters)]))
      
      if(scenario_type == "SC") {
        pc       <- prcomp(X_f_tr, scale. = FALSE)
        pc_var   <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
        n_comp_f <- min(which(pc_var >= 0.90)[1], 10)
        pc_contrib_f <- apply(abs(pc$rotation[, 1:n_comp_f, drop = FALSE]), 1, max)
        f_fold   <- unique(c(f_fold, names(sort(pc_contrib_f, decreasing=TRUE)[1:n_comp_f])))
      }
    }
    
    if(length(f_fold) == 0)
      f_fold <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
    
    X_cv_tr_mat <- as.matrix(X_f_tr[, f_fold, drop = FALSE])
    X_cv_va_mat <- as.matrix(X_f_va[, f_fold, drop = FALSE])
    
    m_cv <- build_model(ncol(X_cv_tr_mat), hidden_list[[1]], 0.1, 0, 0)
    m_cv %>% fit(
      X_cv_tr_mat, y_f_tr,
      epochs       = 40,
      batch_size   = 8,
      class_weight = compute_class_weights(y_f_tr),
      verbose      = 0
    )
    cv_aucs[k] <- as.numeric(roc(y_f_va,
                                 as.vector(predict(m_cv, X_cv_va_mat, verbose = 0)),
                                 quiet = TRUE)$auc)
  }
  mean_cv_auc <- mean(cv_aucs)
  cat(sprintf("  Mean CV AUC: %.4f\n", mean_cv_auc))
  
  final_f <- switch(scenario_type,
                    "S0" = feats_S0,
                    "SA" = feats_SA,
                    "SB" = feats_SB,
                    "SC" = feats_SC)
  
  X_tr_f <- as.matrix(X_train_norm[, final_f, drop = FALSE])
  X_ts_f <- as.matrix(X_test_norm[,  final_f, drop = FALSE])
  
  set.seed(123)
  val_grid_idx <- createDataPartition(y_train, p = 0.20, list = FALSE)
  X_tr_grid    <- X_tr_f[-val_grid_idx, ]
  X_va_grid    <- X_tr_f[ val_grid_idx, ]
  y_tr_grid    <- y_train[-val_grid_idx]
  y_va_grid    <- y_train[ val_grid_idx]
  
  best_val_auc <- -Inf
  best_model   <- NULL
  
  for(cfg in hidden_list) {
    for(dr in c(0.1, 0.2)) {
      m_candidate <- build_model(ncol(X_tr_grid), cfg, dr, 1e-4, 1e-4)
      m_candidate %>% fit(
        X_tr_grid, y_tr_grid,
        epochs       = 40,
        batch_size   = 8,
        class_weight = compute_class_weights(y_tr_grid),
        verbose      = 0
      )
      val_probs <- as.vector(predict(m_candidate, X_va_grid, verbose = 0))
      val_auc   <- as.numeric(roc(y_va_grid, val_probs, quiet = TRUE)$auc)
      cat(sprintf("  hidden=%s, dropout=%.1f -> val AUC=%.4f\n",
                  paste(cfg, collapse="x"), dr, val_auc))
      
      if(val_auc > best_val_auc) {
        best_val_auc <- val_auc
        best_model   <- m_candidate
      }
    }
  }
  
  probs   <- as.vector(predict(best_model, X_ts_f, verbose = 0))
  roc_obj <- roc(y_test, probs, quiet = TRUE)
  ci_auc  <- ci.auc(roc_obj)
  
  # Fixed decision threshold of 0.5 on P(SCC); Sens/Spec/MCC/F1 use ADC as the positive class.
  m <- adc_positive_metrics(as.integer(probs >= 0.5), y_test)   # 0 = ADC (positive), 1 = SCC
  write_preds(out_path(paste0(run_tag, "_preds_", scenario_type, ".csv")),
              y_test, probs, scenario_name)
  
  return(list(
    Name            = scenario_name,
    Features        = length(final_f),
    CV_AUC_NoLeak   = mean_cv_auc,
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
results_list[[1]] <- run_tf_scenario_unbiased("S0", "SCENARIO_0_ALL",            list(c(160,160), c(100,100), c(70,70)))
results_list[[2]] <- run_tf_scenario_unbiased("SA", "SCENARIO_A_LASSO",          list(c(10,10), c(8,8),   c(5,5)))
results_list[[3]] <- run_tf_scenario_unbiased("SB", "SCENARIO_B_LASSO_mRMR",     list(c(20,20),  c(17,17), c(12,12)))
results_list[[4]] <- run_tf_scenario_unbiased("SC", "SCENARIO_C_LASSO_mRMR_PCA", list(c(25,25), c(15,15), c(10,10)))

# ==============================================================================
# 12. FINAL RESULTS & SUMMARY
# ==============================================================================
final_results_tf <- do.call(rbind, lapply(results_list, as.data.frame))
write.csv(final_results_tf, out_path(paste0(run_tag, "_results.csv")), row.names = FALSE)

cat("\n\n========== FINAL SUMMARY ==========\n")
for(i in 1:nrow(final_results_tf)){
  row <- final_results_tf[i,]
  cat(sprintf("%30s: %2d feats | CV-AUC: %.4f | Test-AUC: %.4f [%.4f-%.4f] | Gap: %.4f | MCC: %.3f | F1: %.3f\n",
              row$Name, row$Features, row$CV_AUC_NoLeak, row$Test_AUC,
              row$CI_Low, row$CI_High, row$Overfitting_Gap, row$MCC, row$F1))
}

cat("\n=== GENE SELECTION REPORT ===\n")
cat("Method: estimateSizeFactors + estimateDispersions (train) -> vst(blind=FALSE) -> prcomp -> PC1/PC2 top10\n")
cat("PC1 top 10:", paste(top_pc1, collapse=", "), "\n")
cat("PC2 top 10:", paste(top_pc2, collapse=", "), "\n")
cat(sprintf("Total selected genes: %d\n", length(selected_genes)))
cat(sprintf("Variance explained: PC1=%.1f%%, PC2=%.1f%%\n",
            var_exp[1]*100, var_exp[2]*100))

cat(sprintf("\nAll outputs written to: %s\n", output_dir))
