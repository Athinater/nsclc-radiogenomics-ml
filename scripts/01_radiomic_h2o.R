# ==============================================================================
# 01_radiomic_h2o.R -- radiomic pipeline, H2O deep learning
# Nested 5-fold CV (inner 20% hold-out for tuning, refit) + one evaluation on the test set.
# Outputs: output/radiomic_h2o_results.csv, radiomic_h2o_cv_folds.csv, radiomic_h2o_preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_h2o.R")

prep <- prepare_radiomic()
h2o_start()

hidden <- list(S0 = list(c(100, 100), c(70, 70)),
               SA = list(c(8, 8),     c(5, 5)),
               SB = list(c(20, 20),   c(12, 12)),
               SC = list(c(25, 25),   c(15, 15)))
run_all_scenarios(prep, "H2O", lapply(hidden, h2o_grid), h2o_fit_predict,
                  run_tag = "radiomic_h2o", cleanup = h2o_cleanup)
