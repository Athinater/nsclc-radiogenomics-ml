# ==============================================================================
# models_h2o.R -- H2O deep learning classifier (used by scripts 01, 03, 05)
# Two hidden layers, RectifierWithDropout activation (H2O default hidden dropout 0.5),
# 40 epochs without early stopping, inverse-frequency class weights (observation weights,
# the same weights as in the TensorFlow models; balance_classes = FALSE),
# reproducible = TRUE (single-threaded, deterministic) and seed 123 for every model.
# Grid: hidden-layer sizes (per pipeline and scenario) x input dropout ratio {0.0, 0.1}.
# ==============================================================================
suppressPackageStartupMessages(library(h2o))

h2o_epochs <- 40

h2o_start <- function() {
  h2o.init(strict_version_check = FALSE)
  h2o.no_progress()
}
h2o_cleanup <- function() invisible(h2o.removeAll())
software_note <- function() paste("H2O cluster version:", tryCatch(h2o.getVersion(), error = function(e) NA))

# hidden_list: list of hidden-layer sizes, e.g. list(c(100, 100), c(70, 70))
h2o_grid <- function(hidden_list) {
  force(hidden_list)
  function(n_features) {
    cfgs <- list()
    for (h in hidden_list) for (idr in c(0.0, 0.1)) {
      cfgs[[length(cfgs) + 1]] <- list(hidden = h, input_dropout = idr,
                                       label = sprintf("hidden=%s;input_dropout=%.1f", paste(h, collapse = "x"), idr))
    }
    cfgs
  }
}

h2o_fit_once <- function(X_a, y_a, X_b, cfg) {
  lab <- factor(ifelse(y_a == 1, "SCC", "ADC"), levels = c("ADC", "SCC"))
  n   <- length(y_a)
  w   <- ifelse(y_a == 1, n / (2 * sum(y_a == 1)), n / (2 * sum(y_a == 0)))   # inverse-frequency class weights (as in TF)
  tr  <- as.h2o(data.frame(X_a, Histology = lab, class_weight = w, check.names = FALSE))
  nw  <- as.h2o(data.frame(X_b, check.names = FALSE))
  m <- h2o.deeplearning(
    x = colnames(X_a), y = "Histology", training_frame = tr, weights_column = "class_weight",
    activation = "RectifierWithDropout", hidden = cfg$hidden, input_dropout_ratio = cfg$input_dropout,
    epochs = h2o_epochs, balance_classes = FALSE,
    stopping_rounds = 0, overwrite_with_best_model = FALSE,   # train exactly 40 epochs
    reproducible = TRUE, seed = seed_value)
  p <- as.data.frame(h2o.predict(m, nw))[["SCC"]]
  h2o.rm(list(m, tr, nw))
  p
}

# A model is refitted (up to 3 attempts) if the call to the local H2O server fails, e.g. after a
# connection timeout. With reproducible = TRUE a repeated fit gives the same model.
h2o_fit_predict <- function(X_a, y_a, X_b, cfg, attempts = 3) {
  for (i in seq_len(attempts)) {
    res <- tryCatch(h2o_fit_once(X_a, y_a, X_b, cfg), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (i < attempts) {
      cat("    ", cfg$label, ": attempt ", i, " failed (", conditionMessage(res), "), retrying\n", sep = "")
      Sys.sleep(10)
    }
  }
  stop(res)
}
