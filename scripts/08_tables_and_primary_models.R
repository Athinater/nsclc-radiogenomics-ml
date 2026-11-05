# ==============================================================================
# 08_tables_and_primary_models.R
# Builds the result tables of the manuscript from the outputs of scripts 01-07 and
# selects the primary model of each pipeline. All test-set metrics are RECOMPUTED from
# the prediction files (ADC positive, fixed threshold 0.5) and checked against the
# results files, so the tables, the figures and the post-hoc analyses use the same numbers.
#
# Selection rules (pre-specified; the test set is not used):
#   Table 2 (one row per pipeline x scenario): the classifier with the highest CV-AUC.
#   Primary model per pipeline: the pipeline/scenario/classifier with the highest CV-AUC.
#   Ties: fewer features, then H2O < TF < RF.
#
# Outputs (output/):
#   table_S1_all_models.csv / .tex   all pipeline x scenario x classifier combinations
#   table2_best_per_scenario.csv     Table 2 data;  table2_rows.tex, table3_rows.tex (LaTeX rows)
#   primary_models.csv               used by scripts 09-12
#   feature_selection_frequency.csv  selection frequency of every feature (final models and CV folds)
# ==============================================================================
source("scripts/analysis_helpers.R")

scenarios_used <- scenarios_to_run      # set in common.R (e.g. drop "SC" there)

# ---- Collect and recompute ------------------------------------------------------------
rows <- list(); fold_rows <- list()
for (p in pipelines) for (clf in names(classifiers)) {
  tag <- paste0(p, "_", classifiers[[clf]])
  f   <- out_path(paste0(tag, "_results.csv"))
  if (!file.exists(f)) { warning("Missing results file (script not run?): ", f); next }
  r <- read.csv(f, stringsAsFactors = FALSE)
  r <- r[r$Scenario %in% scenarios_used, , drop = FALSE]
  r$n_test <- NA_integer_
  for (i in seq_len(nrow(r))) {
    d  <- read_preds(r$Preds_file[i])
    ev <- evaluate_test(d$true_label, d$prob_SCC)
    diffs <- c(AUC = ev$AUC - r$Test_AUC[i], TP = ev$TP - r$TP[i], FN = ev$FN - r$FN[i],
               TN = ev$TN - r$TN[i], FP = ev$FP - r$FP[i])
    if (any(abs(diffs) > 1e-6))
      warning(tag, " ", r$Scenario[i], ": results file and prediction file disagree (",
              paste(names(diffs)[abs(diffs) > 1e-6], collapse = ", "), "). Values from the prediction file are used.")
    for (nm in c("CI_Low", "CI_High", "Sensitivity", "Specificity", "MCC", "F1", "TP", "FN", "TN", "FP"))
      r[[nm]][i] <- ev[[nm]]
    r$Test_AUC[i] <- ev$AUC
    r$Overfitting_Gap[i] <- r$CV_AUC[i] - ev$AUC
    r$n_test[i] <- nrow(d)
  }
  rows[[tag]] <- r
  ff <- out_path(paste0(tag, "_cv_folds.csv"))
  if (file.exists(ff)) fold_rows[[tag]] <- read.csv(ff, stringsAsFactors = FALSE)
}
if (length(rows) == 0) stop("No results files found in ", output_dir)
all_res <- do.call(rbind, rows)
all_res$Pipeline_order   <- match(all_res$Pipeline, pipelines)
all_res$Classifier_order <- match(all_res$Classifier, names(classifiers))
all_res$Scenario_order   <- match(all_res$Scenario, scenarios_all)
all_res$Label <- mapply(model_label, all_res$Pipeline, all_res$Classifier, all_res$Scenario)
all_res <- all_res[order(all_res$Pipeline_order, all_res$Scenario_order, all_res$Classifier_order), ]

n_test <- unique(all_res$n_test)
if (length(n_test) != 1) warning("Test-set size differs between models: ", paste(n_test, collapse = ", "))

# ---- Selection by CV-AUC -------------------------------------------------------------------
rank_cv <- function(d) d[order(-d$CV_AUC, d$Features, d$Classifier_order), , drop = FALSE]
best_per_scenario <- do.call(rbind, lapply(split(all_res, list(all_res$Pipeline, all_res$Scenario), drop = TRUE),
                                           function(d) rank_cv(d)[1, ]))
best_per_scenario <- best_per_scenario[order(best_per_scenario$Pipeline_order, best_per_scenario$Scenario_order), ]
primary <- do.call(rbind, lapply(split(all_res, all_res$Pipeline), function(d) rank_cv(d)[1, ]))
primary <- primary[order(primary$Pipeline_order), ]
best_per_scenario$Primary <- paste(best_per_scenario$Pipeline, best_per_scenario$Scenario, best_per_scenario$Classifier) %in%
                             paste(primary$Pipeline, primary$Scenario, primary$Classifier)
all_res$Primary <- paste(all_res$Pipeline, all_res$Scenario, all_res$Classifier) %in%
                   paste(primary$Pipeline, primary$Scenario, primary$Classifier)

# ---- Write CSV tables ------------------------------------------------------------------------
keep <- c("Pipeline", "Scenario", "Classifier", "Features", "Valid_Folds", "CV_AUC", "CV_AUC_SD",
          "Test_AUC", "CI_Low", "CI_High", "Overfitting_Gap", "Sensitivity", "Specificity", "MCC", "F1",
          "TP", "FN", "TN", "FP", "Threshold", "Selected_config", "Primary", "Preds_file", "Feature_list")
write.csv(all_res[, keep], out_path("table_S1_all_models.csv"), row.names = FALSE)
write.csv(best_per_scenario[, keep], out_path("table2_best_per_scenario.csv"), row.names = FALSE)
write.csv(data.frame(Pipeline = primary$Pipeline, Classifier = primary$Classifier, Scenario = primary$Scenario,
                     Label = primary$Label, Run_tag = primary$Run_tag, Preds_file = primary$Preds_file,
                     Features = primary$Features, CV_AUC = primary$CV_AUC, Test_AUC = primary$Test_AUC,
                     CI_Low = primary$CI_Low, CI_High = primary$CI_High),
          out_path("primary_models.csv"), row.names = FALSE)

# ---- LaTeX rows (paste between \midrule and \bottomrule) ---------------------------------------
f3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
scen_tex <- function(d, mark) paste0(d$Scenario, ifelse(d$Primary, paste0("\\textsuperscript{", mark, "}"), ""))
auc_tex  <- function(d) sprintf("%s [%s--%s]", f3(d$Test_AUC), f3(d$CI_Low), f3(d$CI_High))
t2 <- best_per_scenario
table2_tex <- sprintf("%s & %s & %s & %d & %s & %s & %s & %s & %s & %s\\\\",
                      pipeline_names[t2$Pipeline], scen_tex(t2, "e"), t2$Classifier, t2$Features,
                      f3(t2$CV_AUC), auc_tex(t2), f3(t2$Sensitivity), f3(t2$Specificity), f3(t2$MCC), f3(t2$F1))
writeLines(c("% Table 2 rows. Columns: Pipeline & Scenario & Classifier & Features & CV-AUC & Test AUC (95% CI) & Sens. & Spec. & MCC & F1",
             "% Classifier per scenario and primary model (superscript e) selected by the highest CV-AUC.",
             table2_tex), out_path("table2_rows.tex"))
table3_tex <- sprintf("%s & %s & %d & %d & %d & %d\\\\",
                      pipeline_names[t2$Pipeline], scen_tex(t2, "c"), t2$TP, t2$FN, t2$TN, t2$FP)
writeLines(c("% Table 3 rows (same models as Table 2). Columns: Pipeline & Scenario & TP (ADC) & FN (ADC) & TN (SCC) & FP (SCC)",
             table3_tex), out_path("table3_rows.tex"))
s1 <- all_res
tableS1_tex <- sprintf("%s & %s & %s & %d & %s & %s & %s & %s & %s & %s & %d/%d/%d/%d\\\\",
                       pipeline_names[s1$Pipeline], s1$Scenario, s1$Classifier, s1$Features, f3(s1$CV_AUC),
                       auc_tex(s1), f3(s1$Sensitivity), f3(s1$Specificity), f3(s1$MCC), f3(s1$F1),
                       s1$TP, s1$FN, s1$TN, s1$FP)
writeLines(c("% Table S1 rows. Columns: Pipeline & Scenario & Classifier & Features & CV-AUC & Test AUC (95% CI) & Sens. & Spec. & MCC & F1 & TP/FN/TN/FP",
             tableS1_tex), out_path("table_S1_all_models.tex"))

# ---- Feature-selection frequency ---------------------------------------------------------------
# Final models: one feature list per pipeline x selection scenario (SA/SB/SC); the lists are the
# same for the three classifiers because feature selection does not depend on the classifier.
# CV folds: one list per pipeline x scenario x fold.
is_radiomic <- function(f) grepl(radiomic_pattern, f)
sel_scen <- setdiff(scenarios_used, "S0")
count_lists <- function(lists, combos_pipeline) {
  feats <- sort(unique(unlist(lists)))
  data.frame(Feature = feats,
             Count = vapply(feats, function(f) sum(vapply(lists, function(l) f %in% l, logical(1))), numeric(1)),
             Eligible = vapply(feats, function(f) {
               ok <- if (is_radiomic(f)) c("radiomic", "radiogenomic") else c("radiogenomic", "genomic")
               sum(combos_pipeline %in% ok)
             }, numeric(1)), row.names = NULL)
}
fin <- unique(all_res[all_res$Scenario %in% sel_scen, c("Pipeline", "Scenario", "Classifier", "Feature_list")])
fin_lists <- list(); fin_pipe <- character(0)
for (key in unique(paste(fin$Pipeline, fin$Scenario))) {
  fl <- unique(fin$Feature_list[paste(fin$Pipeline, fin$Scenario) == key])
  if (length(fl) > 1) warning("Final feature lists differ between classifiers for ", key, "; the union is counted.")
  fin_lists[[key]] <- unique(unlist(strsplit(fl, ";")))
  fin_pipe <- c(fin_pipe, sub(" .*", "", key))
}
freq <- count_lists(fin_lists, fin_pipe)
names(freq)[2:3] <- c("Final_models", "Final_eligible")
if (length(fold_rows) > 0) {
  fr <- do.call(rbind, fold_rows)
  fr <- unique(fr[fr$Scenario %in% sel_scen, c("Pipeline", "Scenario", "Fold", "Feature_list")])
  fr <- fr[!duplicated(paste(fr$Pipeline, fr$Scenario, fr$Fold)), ]
  ff <- count_lists(lapply(fr$Feature_list, function(x) strsplit(x, ";")[[1]]), fr$Pipeline)
  names(ff)[2:3] <- c("CV_folds", "CV_eligible")
  freq <- merge(freq, ff, by = "Feature", all = TRUE)
}
freq$Type <- ifelse(is_radiomic(freq$Feature), "radiomic", "gene")
freq[is.na(freq)] <- 0
ord_cols <- intersect(c("Final_models", "CV_folds"), names(freq))
freq <- freq[do.call(order, c(list(freq$Type), lapply(ord_cols, function(cn) -freq[[cn]]), list(freq$Feature))), ]
write.csv(freq, out_path("feature_selection_frequency.csv"), row.names = FALSE)

# ---- Console summary --------------------------------------------------------------------------
cat("\nTable 2 (classifier per scenario = highest CV-AUC):\n")
print(data.frame(Pipeline = t2$Pipeline, Scen = t2$Scenario, Clf = t2$Classifier, Feat = t2$Features,
                 CV_AUC = round(t2$CV_AUC, 3), Test_AUC = sprintf("%s [%s-%s]", f3(t2$Test_AUC), f3(t2$CI_Low), f3(t2$CI_High)),
                 Sens = round(t2$Sensitivity, 3), Spec = round(t2$Specificity, 3), MCC = round(t2$MCC, 3),
                 F1 = round(t2$F1, 3), Primary = t2$Primary), row.names = FALSE)
cat("\nPrimary models (highest CV-AUC per pipeline):\n")
print(primary[, c("Label", "Features", "CV_AUC", "Test_AUC", "CI_Low", "CI_High")], row.names = FALSE, digits = 3)
cat("\nSaved tables to", output_dir, "\n")
