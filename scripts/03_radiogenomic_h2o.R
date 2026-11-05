# ==============================================================================
# 03_radiogenomic_h2o.R -- radiogenomic pipeline (radiomic features + 18 genes), H2O
# Genes: DESeq2 VST + PCA loadings on the training set; re-selected inside every CV fold.
# Outputs: output/radiogenomic_h2o_results.csv, _cv_folds.csv, _preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_h2o.R")

prep <- prepare_omics("radiogenomic")[["radiogenomic"]]
h2o_start()

hidden <- list(S0 = list(c(100, 100), c(70, 70)),
               SA = list(c(8, 8),     c(5, 5)),
               SB = list(c(20, 20),   c(12, 12)),
               SC = list(c(25, 25),   c(15, 15)))
run_all_scenarios(prep, "H2O", lapply(hidden, h2o_grid), h2o_fit_predict,
                  run_tag = "radiogenomic_h2o", cleanup = h2o_cleanup)
