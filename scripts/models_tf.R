# ==============================================================================
# models_tf.R -- TensorFlow/Keras classifier (used by scripts 02, 04, 06)
# Dense(ReLU, L1/L2 = 1e-4) -> batch normalisation -> dropout -> Dense(ReLU, L1/L2 = 1e-4)
# -> dropout -> sigmoid; Adam (learning rate 0.001), binary cross-entropy, batch size 8,
# 40 epochs, inverse-frequency class weights. The same architecture is used in all three
# pipelines. Grid: hidden-layer sizes (per pipeline and scenario) x dropout {0.1, 0.2}.
# Reproducibility: op determinism is enabled once, and keras3::set_random_seed(123) is
# called before every model is built.
# ==============================================================================
suppressPackageStartupMessages({
  library(keras3)
  library(tensorflow)
})
tf$config$experimental$enable_op_determinism()

software_note <- function() c(
  paste("TensorFlow version:", tryCatch(as.character(tf$version$VERSION), error = function(e) NA)),
  paste("Python version:", tryCatch(as.character(reticulate::py_config()$version), error = function(e) NA)))

tf_epochs     <- 40
tf_batch_size <- 8

tf_grid <- function(hidden_list) {
  force(hidden_list)
  function(n_features) {
    cfgs <- list()
    for (h in hidden_list) for (dr in c(0.1, 0.2)) {
      cfgs[[length(cfgs) + 1]] <- list(units = h, dropout = dr,
                                       label = sprintf("units=%s;dropout=%.1f", paste(h, collapse = "x"), dr))
    }
    cfgs
  }
}

tf_class_weights <- function(y) {
  n <- length(y)
  list("0" = n / (2 * sum(y == 0)), "1" = n / (2 * sum(y == 1)))
}

tf_build_model <- function(n_features, units, dropout, l1 = 1e-4, l2 = 1e-4) {
  model <- keras_model_sequential(input_shape = c(n_features)) %>%
    layer_dense(units = units[1], activation = "relu",
                kernel_regularizer = regularizer_l1_l2(l1 = l1, l2 = l2)) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = units[2], activation = "relu",
                kernel_regularizer = regularizer_l1_l2(l1 = l1, l2 = l2)) %>%
    layer_dropout(rate = dropout) %>%
    layer_dense(units = 1, activation = "sigmoid")
  model %>% compile(optimizer = optimizer_adam(learning_rate = 0.001),
                    loss = "binary_crossentropy")
  model
}

tf_fit_predict <- function(X_a, y_a, X_b, cfg) {
  clear_session()
  set_random_seed(seed_value)
  m <- tf_build_model(ncol(X_a), cfg$units, cfg$dropout)
  m %>% fit(as.matrix(X_a), as.numeric(y_a), epochs = tf_epochs, batch_size = tf_batch_size,
            class_weight = tf_class_weights(y_a), verbose = 0)
  as.vector(predict(m, as.matrix(X_b), verbose = 0))
}
