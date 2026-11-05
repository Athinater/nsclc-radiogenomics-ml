# ==============================================================================
# 04_radiogenomic_tf.R -- radiogenomic pipeline (radiomic features + 18 genes), TensorFlow
# Genes: DESeq2 VST + PCA loadings on the training set; re-selected inside every CV fold.
# Outputs: output/radiogenomic_tf_results.csv, _cv_folds.csv, _preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_tf.R")

prep <- prepare_omics("radiogenomic")[["radiogenomic"]]

hidden <- list(S0 = list(c(160, 160), c(100, 100), c(70, 70)),
               SA = list(c(10, 10),   c(8, 8),     c(5, 5)),
               SB = list(c(20, 20),   c(17, 17),   c(12, 12)),
               SC = list(c(25, 25),   c(15, 15),   c(10, 10)))
run_all_scenarios(prep, "TF", lapply(hidden, tf_grid), tf_fit_predict, run_tag = "radiogenomic_tf")
