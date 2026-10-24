# ==============================================================================
# analysis_helpers.R -- shared helpers for the post-hoc analyses (scripts 08-12)
# These scripts use the test-set prediction files written by scripts 01-07
# (output/<run_tag>_preds_<scenario>.csv). They do not retrain any model.
# ==============================================================================
output_dir <- "output"                       # same location as in common.R
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- function(f) file.path(output_dir, f)

# ---- Models analysed ------------------------------------------------------------
# One prediction file per model. The primary models are those reported in Table 2:
#   radiomic:      TensorFlow, scenario SA
#   radiogenomic:  TensorFlow, scenario SB / SC (same test AUC in Table 2; both are analysed)
#   genomic:       H2O deep learning, scenario SC
analysis_models <- list(
  "Radiomic (TF, SA)"     = "radiomic_tf_preds_SA.csv",
  "Radiogenomic (TF, SB)" = "radiogenomic_tf_preds_SB.csv",
  "Radiogenomic (TF, SC)" = "radiogenomic_tf_preds_SC.csv",
  "Genomic (H2O, SC)"     = "genomic_h2o_preds_SC.csv"
)

# ---- Event definition -------------------------------------------------------------
# The positive class is ADC, as in Tables 2 and 3 of the manuscript. Brier score and the
# calibration slope do not depend on this choice; the calibration intercept changes sign
# and decision-curve results change with it. Set to "SCC" to repeat the analyses with SCC
# as the event.
event_class <- "ADC"

# Returns the test-set labels and predicted probabilities for one prediction file:
#   y = 1 for the event class, p = predicted probability of the event class.
read_event_preds <- function(file, event = event_class) {
  d <- read.csv(out_path(file))
  stopifnot(all(c("true_label", "prob_SCC", "prob_ADC") %in% names(d)))
  if (event == "ADC") {
    list(y = as.integer(d$true_label == 0), p = d$prob_ADC)
  } else {
    list(y = as.integer(d$true_label == 1), p = d$prob_SCC)
  }
}
