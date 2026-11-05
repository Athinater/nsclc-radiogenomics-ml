# ==============================================================================
# 02_radiomic_tf.R -- radiomic pipeline, TensorFlow/Keras
# Nested 5-fold CV (inner 20% hold-out for tuning, refit) + one evaluation on the test set.
# Outputs: output/radiomic_tf_results.csv, radiomic_tf_cv_folds.csv, radiomic_tf_preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_tf.R")

prep <- prepare_radiomic()

hidden <- list(S0 = list(c(100, 100), c(70, 70)),
               SA = list(c(8, 8),     c(5, 5)),
               SB = list(c(20, 20),   c(12, 12)),
               SC = list(c(25, 25),   c(15, 15)))
run_all_scenarios(prep, "TF", lapply(hidden, tf_grid), tf_fit_predict, run_tag = "radiomic_tf")
