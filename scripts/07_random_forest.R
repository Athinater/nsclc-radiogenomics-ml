# ==============================================================================
# 07_random_forest.R -- Random Forest, all three pipelines
# Same splits, folds, genes and feature selection as scripts 01-06.
# Outputs: output/<pipeline>_rf_results.csv, _cv_folds.csv, _preds_<S>.csv
#          (pipeline = radiomic, radiogenomic, genomic)
# ==============================================================================
source("scripts/common.R")
source("scripts/models_rf.R")

rf_grids <- setNames(rep(list(rf_grid), length(scenarios_all)), scenarios_all)

prep_rad <- prepare_radiomic()
run_all_scenarios(prep_rad, "RF", rf_grids, rf_fit_predict, run_tag = "radiomic_rf")

preps <- prepare_omics(c("radiogenomic", "genomic"))
run_all_scenarios(preps$radiogenomic, "RF", rf_grids, rf_fit_predict, run_tag = "radiogenomic_rf")
run_all_scenarios(preps$genomic,      "RF", rf_grids, rf_fit_predict, run_tag = "genomic_rf")
