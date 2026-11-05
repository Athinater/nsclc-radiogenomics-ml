# ==============================================================================
# analysis_helpers.R -- shared helpers for the post-hoc analyses (scripts 08-12)
# These scripts read the files written by scripts 01-07 (output/<run_tag>_results.csv and
# output/<run_tag>_preds_<scenario>.csv). They do not retrain any model.
# ==============================================================================
source("scripts/common.R")   # paths, class coding, roc_adc(), evaluate_test(), threshold

pipelines      <- c("radiomic", "radiogenomic", "genomic")
pipeline_names <- c(radiomic = "Radiomic", radiogenomic = "Radiogenomic", genomic = "Genomic")
classifiers    <- c(H2O = "h2o", TF = "tf", RF = "rf")      # label -> run_tag suffix

# ---- Event definition -------------------------------------------------------------
# The positive class is ADC, as in Tables 2 and 3. AUC, the Brier score and the calibration
# slope do not depend on this choice; the calibration intercept changes sign and the decision
# curves change with it. Set to "SCC" to repeat the analyses with SCC as the event.
event_class <- "ADC"

read_preds <- function(file) {
  d <- read.csv(out_path(file), stringsAsFactors = FALSE)
  need <- c("patient_id", "true_label", "prob_SCC", "prob_ADC")
  if (!all(need %in% names(d))) stop(file, " lacks columns: ", paste(setdiff(need, names(d)), collapse = ", "))
  if (any(abs(d$prob_ADC - (1 - d$prob_SCC)) > 1e-9)) stop("prob_ADC is not 1 - prob_SCC in ", file)
  d$patient_id <- norm_id(d$patient_id)
  if (anyDuplicated(d$patient_id)) stop("Duplicated patient IDs in ", file)
  d
}

# id, y (1 = event class), p (predicted probability of the event class)
read_event_preds <- function(file, event = event_class) {
  d <- read_preds(file)
  if (event == "ADC") {
    data.frame(id = d$patient_id, y = as.integer(d$true_label == 0), p = d$prob_ADC)
  } else {
    data.frame(id = d$patient_id, y = as.integer(d$true_label == 1), p = d$prob_SCC)
  }
}

# Two models matched by patient ID (columns y, p_1, p_2)
pair_by_id <- function(a, b) {
  m <- merge(a, b, by = "id", suffixes = c("_1", "_2"))
  if (any(m$y_1 != m$y_2)) stop("True labels differ between the two prediction files for the same patient IDs.")
  data.frame(id = m$id, y = m$y_1, p_1 = m$p_1, p_2 = m$p_2)
}

# Primary model per pipeline, written by script 08
load_primary_models <- function() {
  f <- out_path("primary_models.csv")
  if (!file.exists(f)) stop("Run scripts/08_tables_and_primary_models.R first (", f, " not found).")
  pm <- read.csv(f, stringsAsFactors = FALSE)
  pm[match(intersect(pipelines, pm$Pipeline), pm$Pipeline), ]
}

model_label <- function(pipeline, classifier, scenario)
  sprintf("%s (%s, %s)", pipeline_names[[pipeline]], classifier, scenario)

# Colour-blind-safe colours (Okabe-Ito) per pipeline
pipeline_colours <- c(radiomic = "#E69F00", radiogenomic = "#0072B2", genomic = "#009E73")
