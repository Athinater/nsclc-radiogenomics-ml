# ==============================================================================
# RADIOGENOMICS PIPELINE: NESTED CV | H2O
# Gene selection: PCA loadings (PC1 top10 + PC2 top10) on TRAINING SET ONLY
# Zero leakage: test set never seen during gene selection
# Scenarios: S0, SA, SB, SC | Metrics: AUC, CI, Gap, MCC, F1, Sens, Spec
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
run_tag <- "radiogenomic_h2o"

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

# Load clinical metadata (DESeq2 colData)
coldata_full <- read.csv(
  clinical_meta_file,
  header = TRUE, sep = ";", stringsAsFactors = FALSE, row.names = 1
)
rownames(coldata_full) <- toupper(trimws(rownames(coldata_full)))
coldata_full <- coldata_full[coldata_full$Histology %in% c("Adenocarcinoma", "Squamous cell carcinoma"), ]

# Samples common to the counts and the clinical metadata
common_samples <- intersect(rownames(coldata_full), colnames(tbl))
coldata_full   <- coldata_full[common_samples, , drop = FALSE]
tbl_filtered   <- tbl[, common_samples, drop = FALSE]

# ==============================================================================
# 4. STRATIFIED SPLIT (70/30) -- based on clinical_data (radiomic samples)
# ==============================================================================
set.seed(123)
trainIndex <- createDataPartition(clinical_data[[response_variable]], p = 0.7, list = FALSE)

trainResp <- clinical_data[trainIndex,  response_variable]
testResp  <- clinical_data[-trainIndex, response_variable]

# Sample IDs for train/test
train_ids <- rownames(clinical_data)[trainIndex]
test_ids  <- rownames(clinical_data)[-trainIndex]

# ==============================================================================
# 5. GENE SELECTION: PCA LOADINGS ON TRAINING SET ONLY
# ==============================================================================
cat("\n=== Gene Selection: PCA on Training Set Only ===\n")

train_ids_upper <- toupper(trimws(train_ids))
test_ids_upper  <- toupper(trimws(test_ids))
train_rna_ids   <- intersect(train_ids_upper, colnames(tbl_filtered))
test_rna_ids    <- intersect(test_ids_upper,  colnames(tbl_filtered))

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
vst_matrix_tr <- t(assay(vsd_train)[selected_genes, , drop = FALSE])  # samples x genes

dds_test_vst <- DESeqDataSetFromMatrix(tbl_test, coldata_test, ~ Histology)
dds_test_vst <- estimateSizeFactors(dds_test_vst)
dispersions(dds_test_vst) <- dispersions(dds_train)[
  match(rownames(dds_test_vst), rownames(dds_train))
]
vsd_test      <- varianceStabilizingTransformation(dds_test_vst, blind = FALSE)
vst_matrix_te <- t(assay(vsd_test)[selected_genes, , drop = FALSE])   # samples x genes

cat("Train RNA samples:", nrow(vst_matrix_tr), "\n")
cat("Test  RNA samples:", nrow(vst_matrix_te), "\n")

# ==============================================================================
# 6. COMBINE RADIOMIC + GENOMIC FEATURES
# ==============================================================================
# Radiomic features
trainData_radio <- clinical_data[trainIndex, imaging_features]
testData_radio  <- clinical_data[-trainIndex, imaging_features]

# Align genomic with radiomic features (by sample ID)
# Train
tr_radio_ids <- toupper(trimws(rownames(trainData_radio)))
tr_common    <- intersect(tr_radio_ids, rownames(vst_matrix_tr))
trainData_radio_aligned <- trainData_radio[match(tr_common, tr_radio_ids), ]
vst_tr_aligned          <- vst_matrix_tr[tr_common, ]
trainData_raw <- cbind(trainData_radio_aligned, vst_tr_aligned)
trainResp_al  <- clinical_data[[response_variable]][trainIndex][match(tr_common, tr_radio_ids)]
trainResp_al  <- droplevels(trainResp_al)

# Test
te_radio_ids <- toupper(trimws(rownames(testData_radio)))
te_common    <- intersect(te_radio_ids, rownames(vst_matrix_te))
testData_radio_aligned <- testData_radio[match(te_common, te_radio_ids), ]
vst_te_aligned         <- vst_matrix_te[te_common, ]
testData_raw  <- cbind(testData_radio_aligned, vst_te_aligned)
testResp_al   <- clinical_data[[response_variable]][-trainIndex][match(te_common, te_radio_ids)]
testResp_al   <- droplevels(testResp_al)

trainResp <- trainResp_al
testResp  <- testResp_al
y_train   <- as.numeric(trainResp) - 1

cat("Final train samples (radio+genomic):", nrow(trainData_raw), "\n")
cat("Final test  samples (radio+genomic):", nrow(testData_raw),  "\n")

# ==============================================================================
# 7. PREPROCESSING (NZV & NORMALIZATION)
# ==============================================================================
nzv <- nearZeroVar(trainData_raw)
if(length(nzv) > 0){
  trainData_raw <- trainData_raw[, -nzv]
  testData_raw  <- testData_raw[, -nzv]
}

train_means <- apply(trainData_raw, 2, mean)
train_sds   <- apply(trainData_raw, 2, sd)
train_sds[train_sds == 0 | is.na(train_sds)] <- 1

X_train_norm <- scale(trainData_raw, center = train_means, scale = train_sds)
X_test_norm  <- scale(testData_raw,  center = train_means, scale = train_sds)

# ==============================================================================
# 8. PRE-COMPUTATION OF FEATURES FOR FINAL MODEL (FULL TRAINING SET)
# ==============================================================================

# Helper: top-N radiomic features by absolute correlation with the class
# Guarantees that the radiogenomic scenarios always contain radiomic features
get_top_radio <- function(X, y, radio_feats, n = 5) {
  radio_avail <- intersect(radio_feats, colnames(X))
  if(length(radio_avail) == 0) return(character(0))
  cors <- abs(cor(X[, radio_avail, drop = FALSE], y)[, 1])
  names(sort(cors, decreasing = TRUE))[1:min(n, length(radio_avail))]
}

radio_clean    <- make.names(imaging_features)
radio_in_train <- intersect(radio_clean, colnames(X_train_norm))

# S0: Radiomic features + selected genomic genes
genes_clean         <- make.names(selected_genes)
feats_S0_candidates <- c(radio_clean, genes_clean)
feats_S0 <- intersect(feats_S0_candidates, colnames(X_train_norm))
cat("S0 features (radiomic + PCA genomic):", length(feats_S0), "\n")

# SA: LASSO with lambda.min + guaranteed top-5 radiomic
cv_l     <- cv.glmnet(as.matrix(X_train_norm), y_train, family = "binomial", alpha = 1, nfolds = 5)
co_SA    <- coef(cv_l, s = cv_l$lambda.min)
feats_SA <- rownames(co_SA)[co_SA[,1] != 0]
feats_SA <- feats_SA[feats_SA != "(Intercept)"]
if(length(feats_SA) == 0)
  feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[,1], decreasing=TRUE)[1:5])
# Guaranteed top-5 radiomic features
top5_radio_ext <- get_top_radio(X_train_norm, y_train, radio_in_train, n = 5)
feats_SA <- unique(c(feats_SA, top5_radio_ext))
cat("SA features (LASSO + top5 radiomic guaranteed =", round(cv_l$lambda.min, 5), "):", length(feats_SA), "\n")

# SB: LASSO + mRMR + guaranteed top-5 radiomic
n_mrmr_ext <- max(length(feats_SA), 5)
cat("SB mRMR feature count (= LASSO count):", n_mrmr_ext, "\n")

d_mrmr_ext <- mRMR.data(data = data.frame(X_train_norm, response = as.ordered(y_train)))
m_res_ext  <- mRMR.classic(data = d_mrmr_ext,
                            target_indices = ncol(X_train_norm) + 1,
                            feature_count  = n_mrmr_ext)
feats_SB <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_res_ext@filters)]))
# top-5 radiomic features are already included through feats_SA
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
# 9. NESTED CV FUNCTION (NO LEAKAGE)
# ==============================================================================
# Context for gene re-selection inside every cross-validation fold (see common.R)
fold_ctx <- list(ids = tr_common, tbl_train = tbl_train, coldata_train = coldata_train,
                 radio_cols = radio_in_train)

h2o.init(strict_version_check = FALSE)

run_h2o_scenario_no_leakage <- function(scenario_type, scenario_name, hidden_list) {
  cat(sprintf("\n--- Starting Nested CV for: %s ---\n", scenario_name))

  set.seed(123)
  folds_idx <- createFolds(trainResp, k = 5, list = TRUE)
  cv_aucs   <- numeric(5)

  for(k in 1:5) {
    val_idx <- folds_idx[[k]]; tr_idx <- setdiff(seq_along(trainResp), val_idx)
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
      # Guaranteed top-5 radiomic features (selected within the fold)
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
      # Guaranteed top-5 radiomic features (selected within the fold)
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

    tr_h2o_f <- as.h2o(cbind(as.data.frame(X_f_tr[, f_fold, drop=FALSE]),
                              Histology = trainResp[tr_idx]))
    va_h2o_f <- as.h2o(cbind(as.data.frame(X_f_va[, f_fold, drop=FALSE]),
                              Histology = trainResp[val_idx]))

    m_f <- h2o.deeplearning(
      x = setdiff(colnames(tr_h2o_f), "Histology"), y = "Histology",
      training_frame   = tr_h2o_f,
      validation_frame = va_h2o_f,
      epochs           = 40,
      hidden           = hidden_list[[1]],
      balance_classes  = FALSE,
      seed             = 123
    )
    cv_aucs[k] <- h2o.auc(m_f, valid = TRUE)
    h2o.rm(m_f)
  }
  mean_cv_auc <- mean(cv_aucs)
  cat(sprintf("  Mean CV AUC: %.4f\n", mean_cv_auc))

  # FINAL TRAINING & EXTERNAL TEST
  final_features <- switch(scenario_type,
                            "S0" = feats_S0,
                            "SA" = feats_SA,
                            "SB" = feats_SB,
                            "SC" = feats_SC)

  tr_f_df <- as.data.frame(X_train_norm[, final_features, drop=FALSE])
  ts_f_df <- as.data.frame(X_test_norm[,  final_features, drop=FALSE])
  colnames(tr_f_df) <- make.names(colnames(tr_f_df))
  colnames(ts_f_df) <- make.names(colnames(ts_f_df))

  set.seed(123)
  val_grid_idx <- createDataPartition(trainResp, p = 0.20, list = FALSE)
  tr_grid_df   <- tr_f_df[-val_grid_idx, ]
  va_grid_df   <- tr_f_df[ val_grid_idx, ]
  tr_grid_resp <- trainResp[-val_grid_idx]
  va_grid_resp <- trainResp[ val_grid_idx]

  tr_h2o <- as.h2o(cbind(tr_grid_df, Histology = tr_grid_resp))
  va_h2o <- as.h2o(cbind(va_grid_df, Histology = va_grid_resp))
  ts_h2o <- as.h2o(cbind(ts_f_df,    Histology = testResp))

  predictor_cols <- colnames(tr_grid_df)

  best_m       <- NULL
  best_val_auc <- -Inf

  for(hid in hidden_list) {
    for(idr in c(0.0, 0.1)) {
      tryCatch({
        candidate <- h2o.deeplearning(
          x                   = predictor_cols,
          y                   = "Histology",
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
        val_auc <- h2o.auc(candidate, valid = TRUE)
        cat(sprintf("  hidden=%s, dropout=%.1f -> val AUC=%.4f\n",
                    paste(hid, collapse="x"), idr, val_auc))
        if(val_auc > best_val_auc) {
          if(!is.null(best_m)) h2o.rm(best_m)
          best_val_auc <- val_auc
          best_m       <- candidate
        } else {
          h2o.rm(candidate)
        }
      }, error = function(e) {
        cat(sprintf("  SKIP hidden=%s, dropout=%.1f: %s\n",
                    paste(hid, collapse="x"), idr, e$message))
      })
    }
  }

  if(is.null(best_m)) stop(sprintf("No models built successfully for %s", scenario_name))

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
# 10. RUN ALL SCENARIOS
# ==============================================================================
results_list <- list()
results_list[[1]] <- run_h2o_scenario_no_leakage("S0", "SCENARIO_0_ALL",            list(c(100,100), c(70,70)))
results_list[[2]] <- run_h2o_scenario_no_leakage("SA", "SCENARIO_A_LASSO",          list(c(8,8),     c(5,5)))
results_list[[3]] <- run_h2o_scenario_no_leakage("SB", "SCENARIO_B_LASSO_mRMR",     list(c(20,20),   c(12,12)))
results_list[[4]] <- run_h2o_scenario_no_leakage("SC", "SCENARIO_C_LASSO_mRMR_PCA", list(c(25,25),   c(15,15)))

# ==============================================================================
# 11. FINAL RESULTS & SUMMARY
# ==============================================================================
final_results <- do.call(rbind, lapply(results_list, as.data.frame))
write.csv(final_results, out_path(paste0(run_tag, "_results.csv")), row.names = FALSE)

cat("\n\n========== FINAL SUMMARY (UNBIASED) ==========\n")
for(i in 1:nrow(final_results)){
  row <- final_results[i, ]
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
