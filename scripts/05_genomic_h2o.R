# ==============================================================================
# 05_genomic_h2o.R -- genomic pipeline (18 genes), H2O deep learning
# Genes: DESeq2 VST + PCA loadings on the training set; re-selected inside every CV fold.
# Outputs: output/genomic_h2o_results.csv, _cv_folds.csv, _preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_h2o.R")

prep <- prepare_omics("genomic")[["genomic"]]
h2o_start()

hidden <- list(S0 = list(c(25, 25), c(20, 20)),
               SA = list(c(28, 28), c(25, 25)),
               SB = list(c(20, 20), c(18, 18)),
               SC = list(c(26, 26), c(18, 18)))
run_all_scenarios(prep, "H2O", lapply(hidden, h2o_grid), h2o_fit_predict,
                  run_tag = "genomic_h2o", cleanup = h2o_cleanup)
