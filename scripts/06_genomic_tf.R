# ==============================================================================
# 06_genomic_tf.R -- genomic pipeline (18 genes), TensorFlow/Keras
# Same architecture as scripts 02 and 04 (including batch normalisation).
# Outputs: output/genomic_tf_results.csv, _cv_folds.csv, _preds_<S>.csv
# ==============================================================================
source("scripts/common.R")
source("scripts/models_tf.R")

prep <- prepare_omics("genomic")[["genomic"]]

hidden <- list(S0 = list(c(16, 16), c(10, 10)),
               SA = list(c(8, 8),   c(5, 5)),
               SB = list(c(16, 16), c(10, 10)),
               SC = list(c(16, 16), c(10, 10)))
run_all_scenarios(prep, "TF", lapply(hidden, tf_grid), tf_fit_predict, run_tag = "genomic_tf")
