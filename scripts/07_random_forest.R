# ==============================================================================
# RANDOM FOREST PIPELINE: NESTED CV - RADIOMICS | GENOMICS | RADIOGENOMICS
# Gene selection: PCA loadings (PC1 top10 + PC2 top10) on TRAINING SET ONLY
# Zero leakage: test set never seen during gene selection
# Scenarios: S0, SA, SB, SC | Metrics: AUC, CI, Gap, MCC, F1, Sens, Spec
# ==============================================================================

library(caret)
library(glmnet)
library(pROC)
library(mRMRe)
library(randomForest)
library(DESeq2)
library(data.table)
source("scripts/common.R")

# ==============================================================================
# NESTED CV FUNCTION (shared across all datasets)
# ==============================================================================
run_rf_scenario <- function(scenario_type, scenario_name,
                            X_train_norm, X_test_norm,
                            trainResp, testResp, y_train,
                            feats_S0, feats_SA, feats_SB, feats_SC,
                            radio_in_train = NULL, fold_ctx = NULL) {
  
  cat(sprintf("\n--- Starting Nested CV for: %s ---\n", scenario_name))
  
  get_top_radio <- function(X, y, radio_feats, n = 5) {
    if (is.null(radio_feats)) return(character(0))
    radio_avail <- intersect(radio_feats, colnames(X))
    if (length(radio_avail) == 0) return(character(0))
    cors <- abs(cor(X[, radio_avail, drop=FALSE], y)[, 1])
    names(sort(cors, decreasing=TRUE))[1:min(n, length(radio_avail))]
  }
  
  # ---- PART A: NESTED 5-FOLD CV ----
  set.seed(123)
  folds_idx <- createFolds(trainResp, k=5, list=TRUE)
  cv_aucs   <- numeric(5)
  
  for (k in 1:5) {
    val_idx <- folds_idx[[k]]
    tr_idx  <- setdiff(seq_along(trainResp), val_idx)
    y_f_tr  <- y_train[tr_idx]
    y_f_va  <- y_train[val_idx]
    if (is.null(fold_ctx)) {
      X_f_tr <- X_train_norm[tr_idx, ]
      X_f_va <- X_train_norm[val_idx, ]
    } else {
      # genes re-selected within the fold (see common.R)
      fm <- build_fold_matrices(X_train_norm, tr_idx, val_idx, fold_ctx)
      X_f_tr <- fm$X_f_tr
      X_f_va <- fm$X_f_va
    }
    
    f_fold <- NULL
    
    if (scenario_type == "S0") {
      f_fold <- colnames(X_f_tr)
    }
    
    if (scenario_type == "SA") {
      l_cv   <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
      co     <- coef(l_cv, s=l_cv$lambda.min)
      f_fold <- rownames(co)[co[,1] != 0]
      f_fold <- f_fold[f_fold != "(Intercept)"]
      if (length(f_fold)==0)
        f_fold <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
      if (!is.null(radio_in_train)) {
        top5_radio_fold <- get_top_radio(X_f_tr, y_f_tr, radio_in_train, n=5)
        f_fold <- unique(c(f_fold, top5_radio_fold))
      }
    }
    
    if (scenario_type %in% c("SB","SC")) {
      l_cv    <- cv.glmnet(as.matrix(X_f_tr), y_f_tr, family="binomial", alpha = 1, nfolds = 5)
      co      <- coef(l_cv, s=l_cv$lambda.min)
      f_lasso <- rownames(co)[co[,1] != 0]
      f_lasso <- f_lasso[f_lasso != "(Intercept)"]
      if (length(f_lasso)==0)
        f_lasso <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
      if (!is.null(radio_in_train)) {
        top5_radio_fold <- get_top_radio(X_f_tr, y_f_tr, radio_in_train, n=5)
        f_lasso <- unique(c(f_lasso, top5_radio_fold))
      }
      n_mrmr_fold <- max(length(f_lasso), 5)
      d_m    <- mRMR.data(data=data.frame(X_f_tr, response=as.ordered(y_f_tr)))
      m_r    <- mRMR.classic(data=d_m, target_indices=ncol(X_f_tr)+1,
                             feature_count=n_mrmr_fold)
      f_fold <- unique(c(f_lasso, colnames(X_f_tr)[unlist(m_r@filters)]))
      
      if (scenario_type == "SC") {
        pc       <- prcomp(X_f_tr, scale.=FALSE)
        pc_var   <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
        n_comp_f <- min(which(pc_var >= 0.90)[1], 10)
        pc_contrib_f <- apply(abs(pc$rotation[,1:n_comp_f,drop=FALSE]), 1, max)
        f_fold   <- unique(c(f_fold, names(sort(pc_contrib_f, decreasing=TRUE)[1:n_comp_f])))
      }
    }
    
    if (length(f_fold)==0)
      f_fold <- names(sort(abs(cor(X_f_tr, y_f_tr))[,1], decreasing=TRUE)[1:5])
    
    df_tr_cv <- data.frame(X_f_tr[, f_fold, drop=FALSE],
                           Histology=factor(y_f_tr, levels=c(0,1)))
    colnames(df_tr_cv) <- make.names(colnames(df_tr_cv))
    df_va_cv <- as.data.frame(X_f_va[, f_fold, drop=FALSE])
    colnames(df_va_cv) <- make.names(colnames(df_va_cv))
    
    set.seed(123)
    rf_cv <- randomForest(Histology ~ ., data=df_tr_cv,
                          ntree=500, mtry=max(1, floor(sqrt(length(f_fold)))),
                          maxnodes=20, classwt=c("0"=1, "1"=2), importance=FALSE)
    
    prob_va    <- predict(rf_cv, newdata=df_va_cv, type="prob")[,"1"]
    cv_aucs[k] <- as.numeric(roc(y_f_va, prob_va, quiet=TRUE)$auc)
  }
  
  mean_cv_auc <- mean(cv_aucs)
  cat(sprintf("  Mean CV AUC: %.4f\n", mean_cv_auc))
  
  # ---- PART B: GRID SEARCH ----
  final_features <- switch(scenario_type,
                           "S0"=feats_S0, "SA"=feats_SA,
                           "SB"=feats_SB, "SC"=feats_SC)
  
  tr_f_df <- as.data.frame(X_train_norm[, final_features, drop=FALSE])
  ts_f_df <- as.data.frame(X_test_norm[,  final_features, drop=FALSE])
  colnames(tr_f_df) <- make.names(colnames(tr_f_df))
  colnames(ts_f_df) <- make.names(colnames(ts_f_df))
  
  set.seed(123)
  val_grid_idx  <- createDataPartition(trainResp, p=0.20, list=FALSE)
  tr_grid_df    <- tr_f_df[-val_grid_idx, ]
  va_grid_df    <- tr_f_df[ val_grid_idx, ]
  y_tr_grid     <- y_train[-val_grid_idx]
  y_va_grid     <- y_train[ val_grid_idx]
  
  n_feats       <- ncol(tr_grid_df)
  mtry_grid     <- unique(c(max(1, floor(sqrt(n_feats))), max(1, floor(n_feats/3))))
  maxnodes_grid <- c(10, 20, 30)
  
  best_val_auc <- -Inf
  best_rf      <- NULL
  
  for (mt in mtry_grid) {
    for (mn in maxnodes_grid) {
      set.seed(123)
      cand <- randomForest(
        x=tr_grid_df, y=factor(y_tr_grid, levels=c(0,1)),
        ntree=500, mtry=mt, maxnodes=mn,
        classwt=c("0"=1, "1"=2), importance=FALSE
      )
      prob_val <- predict(cand, newdata=va_grid_df, type="prob")[,"1"]
      val_auc  <- as.numeric(roc(y_va_grid, prob_val, quiet=TRUE)$auc)
      cat(sprintf("  mtry=%d maxnodes=%d -> val AUC=%.4f\n", mt, mn, val_auc))
      if (val_auc > best_val_auc) { best_val_auc <- val_auc; best_rf <- cand }
    }
  }
  
  probs   <- predict(best_rf, newdata=ts_f_df, type="prob")[,"1"]
  roc_obj <- roc(testResp, probs, quiet=TRUE)
  ci_auc  <- ci.auc(roc_obj)
  
  # Metrics with ADC as the positive class (as reported in the manuscript); fixed threshold 0.5
  truth01 <- as.numeric(testResp) - 1          # 0 = ADC, 1 = SCC
  m <- adc_positive_metrics(as.integer(probs >= 0.5), truth01)
  write_preds(out_path(paste0("rf_preds_", scenario_name, ".csv")), truth01, probs, scenario_name)
  
  return(list(
    Name=scenario_name, Features=length(final_features),
    CV_AUC=mean_cv_auc, Test_AUC=as.numeric(roc_obj$auc),
    Overfitting_Gap=mean_cv_auc - as.numeric(roc_obj$auc),
    CI_Low=as.numeric(ci_auc[1]), CI_High=as.numeric(ci_auc[3]),
    Sensitivity=m$Sensitivity,
    Specificity=m$Specificity,
    MCC=m$MCC, F1=m$F1, TP=m$TP, FN=m$FN, TN=m$TN, FP=m$FP, rf_model=best_rf
  ))
}

# ==============================================================================
# HELPER: Radiomics feature precomputation
# ==============================================================================
precompute_features_radio <- function(X_train_norm, y_train, min_mrmr=5) {
  feats_S0 <- colnames(X_train_norm)
  cat("S0 features:", length(feats_S0), "\n")
  
  cv_l     <- cv.glmnet(as.matrix(X_train_norm), y_train, family="binomial", alpha = 1, nfolds = 5)
  co_SA    <- coef(cv_l, s=cv_l$lambda.min)
  feats_SA <- rownames(co_SA)[co_SA[,1] != 0]
  feats_SA <- feats_SA[feats_SA != "(Intercept)"]
  if (length(feats_SA)==0)
    feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[,1], decreasing=TRUE)[1:min_mrmr])
  cat("SA features:", length(feats_SA), "\n")
  
  n_mrmr_ext <- max(length(feats_SA), min_mrmr)
  d_mrmr_ext <- mRMR.data(data=data.frame(X_train_norm, response=as.ordered(y_train)))
  m_res_ext  <- mRMR.classic(data=d_mrmr_ext, target_indices=ncol(X_train_norm)+1,
                             feature_count=n_mrmr_ext)
  feats_SB   <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_res_ext@filters)]))
  cat("SB features:", length(feats_SB), "\n")
  
  pca_ext     <- prcomp(X_train_norm, scale.=FALSE)
  pca_var     <- cumsum(pca_ext$sdev^2)/sum(pca_ext$sdev^2)
  n_comp      <- min(which(pca_var >= 0.90)[1], 10)
  pca_contrib <- apply(abs(pca_ext$rotation[,1:n_comp,drop=FALSE]), 1, max)
  feats_SC    <- unique(c(feats_SB, names(sort(pca_contrib, decreasing=TRUE)[1:n_comp])))
  cat("SC features:", length(feats_SC), "\n")
  
  list(S0=feats_S0, SA=feats_SA, SB=feats_SB, SC=feats_SC)
}

# ==============================================================================
# HELPER: Genomic/Radiogenomic feature precomputation (PCA-based genes)
# ==============================================================================
precompute_features_pca <- function(X_train_norm, y_train,
                                    selected_genes, radio_in_train=NULL, min_mrmr=5) {
  
  get_top_radio <- function(X, y, radio_feats, n=5) {
    if (is.null(radio_feats)) return(character(0))
    radio_avail <- intersect(radio_feats, colnames(X))
    if (length(radio_avail)==0) return(character(0))
    cors <- abs(cor(X[, radio_avail, drop=FALSE], y)[, 1])
    names(sort(cors, decreasing=TRUE))[1:min(n, length(radio_avail))]
  }
  
  genes_clean <- make.names(selected_genes)
  
  # S0: PCA genes (+ all radiomic features for the radiogenomic pipeline)
  feats_S0_cand <- if (!is.null(radio_in_train)) c(radio_in_train, genes_clean) else genes_clean
  feats_S0 <- intersect(feats_S0_cand, colnames(X_train_norm))
  cat("S0 features:", length(feats_S0), "\n")
  
  # SA: LASSO + guaranteed top-5 radiomic
  cv_l     <- cv.glmnet(as.matrix(X_train_norm), y_train, family="binomial", alpha = 1, nfolds = 5)
  co_SA    <- coef(cv_l, s=cv_l$lambda.min)
  feats_SA <- rownames(co_SA)[co_SA[,1] != 0]
  feats_SA <- feats_SA[feats_SA != "(Intercept)"]
  if (length(feats_SA)==0)
    feats_SA <- names(sort(abs(cor(X_train_norm, y_train))[,1], decreasing=TRUE)[1:min_mrmr])
  if (!is.null(radio_in_train))
    feats_SA <- unique(c(feats_SA, get_top_radio(X_train_norm, y_train, radio_in_train, n=5)))
  cat("SA features:", length(feats_SA), "\n")
  
  # SB: LASSO + mRMR
  n_mrmr_ext <- max(length(feats_SA), min_mrmr)
  d_mrmr_ext <- mRMR.data(data=data.frame(X_train_norm, response=as.ordered(y_train)))
  m_res_ext  <- mRMR.classic(data=d_mrmr_ext, target_indices=ncol(X_train_norm)+1,
                             feature_count=n_mrmr_ext)
  feats_SB   <- unique(c(feats_SA, colnames(X_train_norm)[unlist(m_res_ext@filters)]))
  cat("SB features:", length(feats_SB), "\n")
  
  # SC: LASSO + mRMR + PCA
  pca_ext     <- prcomp(X_train_norm, scale.=FALSE)
  pca_var     <- cumsum(pca_ext$sdev^2)/sum(pca_ext$sdev^2)
  n_comp      <- min(which(pca_var >= 0.90)[1], ncol(X_train_norm)-1)
  pca_contrib <- apply(abs(pca_ext$rotation[,1:n_comp,drop=FALSE]), 1, max)
  feats_SC    <- unique(c(feats_SB, names(sort(pca_contrib, decreasing=TRUE)[1:n_comp])))
  cat("SC features:", length(feats_SC), "\n")
  
  list(S0=feats_S0, SA=feats_SA, SB=feats_SB, SC=feats_SC)
}

collect_df <- function(res_list) {
  do.call(rbind, lapply(res_list, function(r) {
    data.frame(Name=r$Name, Features=r$Features,
               CV_AUC=round(r$CV_AUC,4), Test_AUC=round(r$Test_AUC,4),
               Overfitting_Gap=round(r$Overfitting_Gap,4),
               CI_Low=round(r$CI_Low,4), CI_High=round(r$CI_High,4),
               Sensitivity=round(r$Sensitivity,4), Specificity=round(r$Specificity,4),
               MCC=round(r$MCC,4), F1=round(r$F1,4),
               TP=r$TP, FN=r$FN, TN=r$TN, FP=r$FP, stringsAsFactors=FALSE)
  }))
}

# ==============================================================================
# PATHS & CONFIG
# ==============================================================================
path_rad          <- clinical_rad_file
path_rg           <- clinical_rg_file
path_counts       <- counts_file
path_clinical_rna <- clinical_meta_file
response_variable <- "Histology"

# ==============================================================================
# LOAD RAW COUNTS (shared by the genomic and radiogenomic pipelines)
# ==============================================================================
tbl <- load_counts(path_counts)

coldata_full <- read.csv(path_clinical_rna, header=TRUE, sep=";",
                         stringsAsFactors=FALSE, row.names=1)
rownames(coldata_full) <- toupper(trimws(rownames(coldata_full)))
coldata_full <- coldata_full[coldata_full$Histology %in%
                               c("Adenocarcinoma","Squamous cell carcinoma"), ]
common_samples <- intersect(rownames(coldata_full), colnames(tbl))
coldata_full   <- coldata_full[common_samples, , drop=FALSE]
tbl_filtered   <- tbl[, common_samples, drop=FALSE]

# ==============================================================================
# DATASET 1: RADIOMICS (gene selection not applicable)
# ==============================================================================
cat("\n\n========== DATASET 1: RADIOMICS ==========\n")

clinical_data <- read.csv(path_rad, sep=";", header=TRUE, row.names=1)
clinical_data[[response_variable]] <- as.factor(clinical_data[[response_variable]])
levels(clinical_data[[response_variable]]) <- make.names(levels(clinical_data[[response_variable]]))

imaging_features <- colnames(clinical_data)[grepl(
  "MORPHOLOGICAL_|INTENSITY.BASED_|INTENSITY.HISTOGRAM_|GLCM_|GLRLM_|GLSZM_|NGTDM_",
  colnames(clinical_data)
)]

set.seed(123)
trainIndex    <- createDataPartition(clinical_data[[response_variable]], p=0.7, list=FALSE)
trainData_raw <- clinical_data[trainIndex,  imaging_features]
testData_raw  <- clinical_data[-trainIndex, imaging_features]
trainResp     <- clinical_data[trainIndex,  response_variable]
testResp      <- clinical_data[-trainIndex, response_variable]
y_train       <- as.numeric(trainResp) - 1

nzv <- nearZeroVar(trainData_raw)
if (length(nzv)>0) { trainData_raw <- trainData_raw[,-nzv]; testData_raw <- testData_raw[,-nzv] }
train_means <- apply(trainData_raw,2,mean); train_sds <- apply(trainData_raw,2,sd)
train_sds[train_sds==0|is.na(train_sds)] <- 1
X_train_norm <- scale(trainData_raw, center=train_means, scale=train_sds)
X_test_norm  <- scale(testData_raw,  center=train_means, scale=train_sds)

feats <- precompute_features_radio(X_train_norm, y_train)

results_list <- list()
results_list[[1]] <- run_rf_scenario("S0","RAD_S0_ALL",            X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC)
results_list[[2]] <- run_rf_scenario("SA","RAD_SA_LASSO",          X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC)
results_list[[3]] <- run_rf_scenario("SB","RAD_SB_LASSO_mRMR",     X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC)
results_list[[4]] <- run_rf_scenario("SC","RAD_SC_LASSO_mRMR_PCA", X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC)

final_rad <- collect_df(results_list)
write.csv(final_rad, out_path("rf_radiomic_results.csv"), row.names=FALSE)
cat("\n--- Radiomics RF Results ---\n"); print(final_rad)

# ==============================================================================
# SHARED GENE SELECTION: PCA ON TRAINING SET ONLY
# (used by both Dataset 2 Radiogenomics and Dataset 3 Genomics)
# ==============================================================================
cat("\n\n=== Gene Selection: PCA on Training Set Only ===\n")

# Use the same split as the radiogenomics dataset (path_rg)
clinical_data_rg <- read.csv(path_rg, sep=";", header=TRUE, row.names=1)
clinical_data_rg[[response_variable]] <- as.factor(clinical_data_rg[[response_variable]])
levels(clinical_data_rg[[response_variable]]) <- make.names(levels(clinical_data_rg[[response_variable]]))

set.seed(123)
trainIndex_rg <- createDataPartition(clinical_data_rg[[response_variable]], p=0.7, list=FALSE)
train_ids_rg  <- toupper(trimws(rownames(clinical_data_rg)[trainIndex_rg]))
test_ids_rg   <- toupper(trimws(rownames(clinical_data_rg)[-trainIndex_rg]))

train_rna_ids <- intersect(train_ids_rg, colnames(tbl_filtered))
test_rna_ids  <- intersect(test_ids_rg,  colnames(tbl_filtered))

tbl_train_gs     <- tbl_filtered[, train_rna_ids, drop=FALSE]
tbl_test_gs      <- tbl_filtered[, test_rna_ids,  drop=FALSE]
coldata_train_gs <- coldata_full[train_rna_ids, , drop=FALSE]
coldata_test_gs  <- coldata_full[test_rna_ids,  , drop=FALSE]

# DESeq2: fit size factors + dispersions on TRAINING SET ONLY
dds_gs <- DESeqDataSetFromMatrix(tbl_train_gs, coldata_train_gs, ~ Histology)
dds_gs <- estimateSizeFactors(dds_gs)
dds_gs <- estimateDispersions(dds_gs, quiet=TRUE)
vsd_gs <- vst(dds_gs, blind=FALSE)

# PCA loadings -> selected genes
pca_gs  <- prcomp(t(assay(vsd_gs)))
load_gs <- pca_gs$rotation
top_pc1 <- rownames(load_gs)[head(order(abs(load_gs[,"PC1"]), decreasing=TRUE), 10)]
top_pc2 <- rownames(load_gs)[head(order(abs(load_gs[,"PC2"]), decreasing=TRUE), 10)]
selected_genes <- unique(c(top_pc1, top_pc2))
check_selected_genes(selected_genes)

cat("PC1 top 10:", paste(top_pc1, collapse=", "), "\n")
cat("PC2 top 10:", paste(top_pc2, collapse=", "), "\n")
cat("Total selected genes:", length(selected_genes), "\n")

# VST matrices: train from dds_gs; test with test-specific size factors and training dispersions
vst_matrix_tr <- t(assay(vsd_gs)[selected_genes, , drop=FALSE])

dds_test_vst <- DESeqDataSetFromMatrix(tbl_test_gs, coldata_test_gs, ~ Histology)

# test-specific size factors
dds_test_vst <- estimateSizeFactors(dds_test_vst)   # test-specific size factors

# dispersions from the training set
dispersions(dds_test_vst) <- dispersions(dds_gs)[match(rownames(dds_test_vst), rownames(dds_gs))]   # training dispersions

# VST with the training dispersions
vsd_test_vst <- varianceStabilizingTransformation(dds_test_vst, blind=FALSE)
vst_matrix_te <- t(assay(vsd_test_vst)[selected_genes, , drop=FALSE])

cat("VST train:", nrow(vst_matrix_tr), "x", ncol(vst_matrix_tr), "\n")
cat("VST test: ", nrow(vst_matrix_te), "x", ncol(vst_matrix_te), "\n")

# Imaging features (same for both RG datasets)
imaging_features <- colnames(clinical_data_rg)[grepl(
  "MORPHOLOGICAL_|INTENSITY.BASED_|INTENSITY.HISTOGRAM_|GLCM_|GLRLM_|GLSZM_|NGTDM_",
  colnames(clinical_data_rg)
)]

# ==============================================================================
# DATASET 2: RADIOGENOMICS
# ==============================================================================
cat("\n\n========== DATASET 2: RADIOGENOMICS ==========\n")

trainData_radio <- clinical_data_rg[trainIndex_rg, imaging_features]
testData_radio  <- clinical_data_rg[-trainIndex_rg, imaging_features]

tr_radio_ids <- toupper(trimws(rownames(trainData_radio)))
tr_common    <- intersect(tr_radio_ids, rownames(vst_matrix_tr))
trainData_raw <- cbind(trainData_radio[match(tr_common, tr_radio_ids), ],
                       vst_matrix_tr[tr_common, ])
trainResp <- clinical_data_rg[[response_variable]][trainIndex_rg][match(tr_common, tr_radio_ids)]
trainResp <- droplevels(trainResp)

te_radio_ids <- toupper(trimws(rownames(testData_radio)))
te_common    <- intersect(te_radio_ids, rownames(vst_matrix_te))
testData_raw <- cbind(testData_radio[match(te_common, te_radio_ids), ],
                      vst_matrix_te[te_common, ])
testResp <- clinical_data_rg[[response_variable]][-trainIndex_rg][match(te_common, te_radio_ids)]
testResp <- droplevels(testResp)
y_train  <- as.numeric(trainResp) - 1

nzv <- nearZeroVar(trainData_raw)
if (length(nzv)>0) { trainData_raw <- trainData_raw[,-nzv]; testData_raw <- testData_raw[,-nzv] }
train_means <- apply(trainData_raw,2,mean); train_sds <- apply(trainData_raw,2,sd)
train_sds[train_sds==0|is.na(train_sds)] <- 1
X_train_norm <- scale(trainData_raw, center=train_means, scale=train_sds)
X_test_norm  <- scale(testData_raw,  center=train_means, scale=train_sds)

radio_clean    <- make.names(imaging_features)
radio_in_train <- intersect(radio_clean, colnames(X_train_norm))

# Context for gene re-selection inside every cross-validation fold (see common.R)
fold_ctx_rg <- list(ids = tr_common, tbl_train = tbl_train_gs, coldata_train = coldata_train_gs,
                    radio_cols = radio_in_train)

feats <- precompute_features_pca(X_train_norm, y_train, selected_genes,
                                 radio_in_train=radio_in_train)

results_list <- list()
results_list[[1]] <- run_rf_scenario("S0","RG_S0_ALL",            X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, radio_in_train, fold_ctx_rg)
results_list[[2]] <- run_rf_scenario("SA","RG_SA_LASSO",          X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, radio_in_train, fold_ctx_rg)
results_list[[3]] <- run_rf_scenario("SB","RG_SB_LASSO_mRMR",     X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, radio_in_train, fold_ctx_rg)
results_list[[4]] <- run_rf_scenario("SC","RG_SC_LASSO_mRMR_PCA", X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, radio_in_train, fold_ctx_rg)

final_rg <- collect_df(results_list)
write.csv(final_rg, out_path("rf_radiogenomic_results.csv"), row.names=FALSE)
cat("\n--- Radiogenomics RF Results ---\n"); print(final_rg)

# ==============================================================================
# DATASET 3: GENOMICS
# ==============================================================================
cat("\n\n========== DATASET 3: GENOMICS ==========\n")

# (fold_ctx_gen is defined below, once the genomic sample order is known)
tr_common_g <- intersect(train_ids_rg, rownames(vst_matrix_tr))
te_common_g <- intersect(test_ids_rg,  rownames(vst_matrix_te))
fold_ctx_gen <- list(ids = tr_common_g, tbl_train = tbl_train_gs, coldata_train = coldata_train_gs,
                     radio_cols = character(0))

trainData_raw <- as.data.frame(vst_matrix_tr[tr_common_g, , drop=FALSE])
trainResp     <- clinical_data_rg[[response_variable]][trainIndex_rg][
  match(tr_common_g, train_ids_rg)]
trainResp <- droplevels(trainResp)

testData_raw <- as.data.frame(vst_matrix_te[te_common_g, , drop=FALSE])
testResp     <- clinical_data_rg[[response_variable]][-trainIndex_rg][
  match(te_common_g, test_ids_rg)]
testResp <- droplevels(testResp)
y_train  <- as.numeric(trainResp) - 1

nzv <- nearZeroVar(trainData_raw)
if (length(nzv)>0) { trainData_raw <- trainData_raw[,-nzv]; testData_raw <- testData_raw[,-nzv] }
train_means <- apply(trainData_raw,2,mean); train_sds <- apply(trainData_raw,2,sd)
train_sds[train_sds==0|is.na(train_sds)] <- 1
X_train_norm <- scale(trainData_raw, center=train_means, scale=train_sds)
X_test_norm  <- scale(testData_raw,  center=train_means, scale=train_sds)

feats <- precompute_features_pca(X_train_norm, y_train, selected_genes,
                                 radio_in_train=NULL, min_mrmr=3)

results_list <- list()
results_list[[1]] <- run_rf_scenario("S0","GEN_S0_ALL",            X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, fold_ctx = fold_ctx_gen)
results_list[[2]] <- run_rf_scenario("SA","GEN_SA_LASSO",          X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, fold_ctx = fold_ctx_gen)
results_list[[3]] <- run_rf_scenario("SB","GEN_SB_LASSO_mRMR",     X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, fold_ctx = fold_ctx_gen)
results_list[[4]] <- run_rf_scenario("SC","GEN_SC_LASSO_mRMR_PCA", X_train_norm,X_test_norm,trainResp,testResp,y_train,feats$S0,feats$SA,feats$SB,feats$SC, fold_ctx = fold_ctx_gen)

final_gen <- collect_df(results_list)
write.csv(final_gen, out_path("rf_genomic_results.csv"), row.names=FALSE)
cat("\n--- Genomics RF Results ---\n"); print(final_gen)

# ==============================================================================
# COMBINED SUMMARY
# ==============================================================================
all_results <- rbind(final_rad, final_rg, final_gen)
write.csv(all_results, out_path("rf_all_results.csv"), row.names=FALSE)

cat("\n\n========== FINAL SUMMARY - RANDOM FOREST ALL DATASETS ==========\n")
cat(sprintf("Gene selection: PCA loadings (PC1 top10 + PC2 top10) on training set only\n"))
cat(sprintf("Selected genes (%d): %s\n\n", length(selected_genes),
            paste(selected_genes, collapse=", ")))

for (i in seq_len(nrow(all_results))) {
  r <- all_results[i,]
  cat(sprintf(
    "%25s: %3d feats | CV-AUC: %.4f | Test-AUC: %.4f [%.4f-%.4f] | Gap: %+.4f | Sens: %.3f | Spec: %.3f | MCC: %.3f | F1: %.3f\n",
    r$Name, r$Features, r$CV_AUC, r$Test_AUC,
    r$CI_Low, r$CI_High, r$Overfitting_Gap,
    r$Sensitivity, r$Specificity, r$MCC, r$F1
  ))
}

cat("\nMethodological notes:\n")
cat("  Gene selection: estimateSizeFactors + estimateDispersions (train) -> vst(blind=FALSE) -> prcomp -> PC1/PC2 top10\n")
cat("  VST test set:   estimateSizeFactors (test-specific) + dispersions from train -> varianceStabilizingTransformation\n")
cat("  PC1 top10 + PC2 top10 loadings -> selected genes\n")
cat("  Package:    randomForest (Breiman 2001)\n")
cat("  Threshold:  fixed 0.5\n")
cat("  ntree=500 | seed=123\n")