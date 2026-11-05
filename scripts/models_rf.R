# ==============================================================================
# models_rf.R -- Random Forest classifier (used by script 07)
# randomForest, 500 trees, class weights 1:2 (ADC:SCC), seed 123 before every forest.
# Grid: mtry {floor(sqrt(p)), floor(p/3)} x maxnodes {10, 20, 30} (p = number of features).
# ==============================================================================
suppressPackageStartupMessages(library(randomForest))

rf_ntree <- 500

rf_grid <- function(n_features) {
  mtry <- unique(c(max(1, floor(sqrt(n_features))), max(1, floor(n_features / 3))))
  cfgs <- list()
  for (mt in mtry) for (mn in c(10, 20, 30)) {
    cfgs[[length(cfgs) + 1]] <- list(mtry = mt, maxnodes = mn,
                                     label = sprintf("mtry=%d;maxnodes=%d", mt, mn))
  }
  cfgs
}

rf_fit_predict <- function(X_a, y_a, X_b, cfg) {
  set.seed(seed_value)
  m <- randomForest(x = as.data.frame(X_a), y = factor(y_a, levels = c(0, 1)),
                    ntree = rf_ntree, mtry = cfg$mtry, maxnodes = cfg$maxnodes,
                    classwt = c("0" = 1, "1" = 2))
  predict(m, newdata = as.data.frame(X_b), type = "prob")[, "1"]
}
