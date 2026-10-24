# ==============================================================================
# 09_calibration_brier.R
# Brier score and Cox recalibration (intercept and slope, Wald 95% CI) on the test set.
# Manuscript values for comparison: Brier 0.041 (radiogenomic), 0.068 (genomic);
# calibration slope 1.73 (radiogenomic), 0.71 (genomic); both 95% CI include 1.
# ==============================================================================
source("scripts/analysis_helpers.R")

eps <- 1e-6     # probabilities are clipped away from 0 and 1 before taking the logit

calibration_stats <- function(y, p) {
  brier <- mean((p - y)^2)
  lp    <- qlogis(pmin(pmax(p, eps), 1 - eps))
  # Cox recalibration: logit P(y = 1) = a + b * logit(p); ideal a = 0, b = 1
  fit <- suppressWarnings(glm(y ~ lp, family = binomial))
  co  <- summary(fit)$coefficients
  ci  <- suppressWarnings(confint.default(fit))
  # calibration-in-the-large: intercept with the slope fixed at 1
  fit0 <- suppressWarnings(glm(y ~ 1 + offset(lp), family = binomial))
  ci0  <- suppressWarnings(confint.default(fit0))
  data.frame(
    Brier           = round(brier, 3),
    Cox_intercept   = round(co[1, 1], 2),
    Intercept_CI    = sprintf("%.2f to %.2f", ci[1, 1], ci[1, 2]),
    Cox_slope       = round(co[2, 1], 2),
    Slope_CI        = sprintf("%.2f to %.2f", ci[2, 1], ci[2, 2]),
    CITL_intercept  = round(coef(fit0)[1], 2),
    CITL_CI         = sprintf("%.2f to %.2f", ci0[1, 1], ci0[1, 2]),
    n = length(y), n_events = sum(y == 1)
  )
}

rows <- lapply(names(analysis_models), function(nm) {
  d <- read_event_preds(analysis_models[[nm]])
  cbind(Model = nm, calibration_stats(d$y, d$p))
})
out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, out_path("calibration_brier.csv"), row.names = FALSE)
cat("Saved: ", out_path("calibration_brier.csv"), "\n")
cat("Event class:", event_class, "(the sign of the intercept depends on this choice).\n")
cat("With few events the Wald intervals are wide and, for near-perfect separation, unstable.\n")
