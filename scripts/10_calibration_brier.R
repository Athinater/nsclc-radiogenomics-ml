# ==============================================================================
# 10_calibration_brier.R
# Brier score and Cox recalibration (intercept and slope, Wald 95% CI) of the primary
# models (output/primary_models.csv, script 08) on the test set. Event = ADC.
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
    Brier          = round(brier, 3),
    Cox_intercept  = round(co[1, 1], 2),
    Intercept_CI   = sprintf("%.2f to %.2f", ci[1, 1], ci[1, 2]),
    Cox_slope      = round(co[2, 1], 2),
    Slope_CI       = sprintf("%.2f to %.2f", ci[2, 1], ci[2, 2]),
    CITL_intercept = round(coef(fit0)[1], 2),
    CITL_CI        = sprintf("%.2f to %.2f", ci0[1, 1], ci0[1, 2]),
    Separation     = fit$converged == FALSE || any(abs(co[, 1]) > 20),
    n = length(y), n_events = sum(y == 1)
  )
}

pm <- load_primary_models()
out <- do.call(rbind, lapply(seq_len(nrow(pm)), function(i) {
  d <- read_event_preds(pm$Preds_file[i])
  cbind(Model = pm$Label[i], calibration_stats(d$y, d$p))
}))
print(out, row.names = FALSE)
write.csv(out, out_path("calibration_brier.csv"), row.names = FALSE)
cat("Saved:", out_path("calibration_brier.csv"), "\n")
cat("Event class:", event_class, "(the sign of the intercept depends on this choice).\n")
cat("With few events the Wald intervals are wide; 'Separation = TRUE' means the recalibration\n",
    "model did not converge or has extreme coefficients (near-perfect separation), so its\n",
    "intercept/slope should not be interpreted.\n")
