# Radiomic, genomic and radiogenomic machine-learning models for NSCLC subtype classification

Code accompanying the manuscript:

> Terzi A, Dikaios N, Protonotarios NE, Athanasiadis E, Kastis GA.
> *Evaluation of Radiomic, Genomic, and Radiogenomic Machine Learning Models for Non-Small Cell Lung Cancer Subtype Classification.*

The code classifies lung adenocarcinoma (ADC) versus squamous cell carcinoma (SCC) from CT-derived radiomic features, tumour RNA-Seq gene expression, and both combined, using H2O deep learning, TensorFlow/Keras and Random Forest classifiers.

## Data

**No patient data are distributed in this repository.** The analysis uses two public datasets:

| Data | Source | Licence / terms |
|---|---|---|
| CT images, tumour segmentations, clinical data | The Cancer Imaging Archive (TCIA), collection *NSCLC Radiogenomics*, DOI [10.7937/K9/TCIA.2017.7hs46erv](https://doi.org/10.7937/K9/TCIA.2017.7hs46erv) | CC BY 3.0; TCIA Data Usage Policy |
| RNA-Seq raw counts | NCBI Gene Expression Omnibus, accession [GSE103584](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103584) | see GEO |

Please cite:

- Bakr S, et al. *Data for NSCLC Radiogenomics (Version 4) [Data set].* The Cancer Imaging Archive, 2017. doi:10.7937/K9/TCIA.2017.7hs46erv
- Bakr S, et al. A radiogenomic dataset of non-small cell lung cancer. *Sci Data* 2018;5:180202. doi:10.1038/sdata.2018.202

### Radiomic features

Radiomic features (n = 142: first-order, GLCM, GLSZM, GLRLM and NGTDM) were extracted from the TCIA CT images with the provided tumour segmentations using LIFEx v7.6.7 (400 grey levels, 1 mm^3 isotropic resampling, intensity rescaling to [-1000, 3000] HU). The resulting feature tables are **not** included here. To run the radiomic and radiogenomic scripts, regenerate the features from the public images and segmentations as described in the manuscript (Methods).

### Expected files in `data/`

| File | Content |
|---|---|
| `GSE103584_raw_counts_GRCh38.p13_NCBI.tsv` | GEO raw count matrix (first column: Entrez GeneID; one column per sample) |
| `clinicaldatawithout.csv` | Radiomic features + `Histology` (0 = ADC, 1 = SCC), one row per patient (semicolon-separated, patient ID as first column) |
| `clinicaldata.csv` | Same structure, restricted to patients with matched RNA-Seq |
| `onlyclinical.csv` | Clinical metadata incl. `Histology` ("Adenocarcinoma" / "Squamous cell carcinoma"), used as DESeq2 `colData` |

File names are set at the top of `scripts/common.R`. The final cohort has 113 patients (87 ADC, 26 SCC); selection criteria are given in the manuscript.

## Repository structure

```
scripts/
  common.R                         settings, data preparation, feature selection, nested-CV runner
  models_h2o.R / models_tf.R / models_rf.R   classifier definitions and hyperparameter grids
  01_radiomic_h2o.R                radiomic pipeline, H2O deep learning
  02_radiomic_tf.R                 radiomic pipeline, TensorFlow/Keras
  03_radiogenomic_h2o.R            radiogenomic pipeline, H2O deep learning
  04_radiogenomic_tf.R             radiogenomic pipeline, TensorFlow/Keras
  05_genomic_h2o.R                 genomic pipeline, H2O deep learning
  06_genomic_tf.R                  genomic pipeline, TensorFlow/Keras
  07_random_forest.R               Random Forest, all three pipelines
  analysis_helpers.R               helpers for scripts 08-12
  08_tables_and_primary_models.R   Tables 2, 3 and S1, primary model per pipeline, feature-selection frequency
  09_delong_comparison.R           paired DeLong tests between the primary models
  10_calibration_brier.R           Brier score, Cox recalibration intercept/slope
  11_decision_curve.R              decision curve analysis (net benefit) + figure
  12_roc_figure.R                  ROC curves of the primary models + figure
  13_go_enrichment.R               GO enrichment of the selected genes
  14_radiomic_gene_correlations.R  Spearman correlations, radiomic features vs genes
  15_correlation_heatmap.R         heatmap of the correlations (Figure 5)
  16_figure3_roc_heatmap.R         Figure 3: ROC curves per pipeline + heatmap of the threshold metrics
data/                              input files (not distributed; see above)
output/                            created by the scripts
```

## How to run

Run every script from the repository root (the folder containing `scripts/` and `data/`), in this order:

```
Rscript scripts/01_radiomic_h2o.R        # 01-07 are independent of each other
...
Rscript scripts/07_random_forest.R
Rscript scripts/08_tables_and_primary_models.R
Rscript scripts/09_delong_comparison.R   # 09-12 need the output of 08
Rscript scripts/10_calibration_brier.R
Rscript scripts/11_decision_curve.R
Rscript scripts/12_roc_figure.R
Rscript scripts/13_go_enrichment.R
Rscript scripts/14_radiomic_gene_correlations.R   # uses the feature frequencies written by 08
Rscript scripts/15_correlation_heatmap.R          # needs the output of 14
Rscript scripts/16_figure3_roc_heatmap.R          # needs the output of 08 (package patchwork)
```

Each model script writes to `output/`: `<pipeline>_<classifier>_results.csv` (one row per scenario: CV-AUC, test AUC with 95% CI, sensitivity, specificity, MCC, F1, confusion-matrix counts, selected hyperparameters, selected features), `<pipeline>_<classifier>_cv_folds.csv` (per-fold AUC, features and hyperparameters), `<pipeline>_<classifier>_preds_<scenario>.csv` (test-set predictions with patient ID, true label, `prob_SCC`, `prob_ADC`) and `<pipeline>_<classifier>_sessionInfo.txt`. All numbers in the tables and figures are computed from these files.

To drop a feature-selection scenario from all analyses, edit `scenarios_to_run` in `scripts/common.R`.

## Methods implemented

- **Split:** stratified 70/30 split (training n = 80, test n = 33; test set 26 ADC, 7 SCC), seed 123. The test set is used once, for the final evaluation. The test-patient IDs of each pipeline are written to `output/<pipeline>_test_ids.txt`.
- **Preprocessing:** near-zero-variance filter (label-free) and z-scaling with training-set statistics; inside each cross-validation fold the z-scaling is repeated with the fold-training statistics. RNA-Seq counts: genes with >= 10 counts in >= 2 samples (label-free filter on all samples), Entrez ID to symbol mapping with `org.Hs.eg.db`.
- **Genes:** DESeq2 size factors, dispersions and dispersion trend estimated on the training samples; variance-stabilising transformation with that trend (`varianceStabilizingTransformation(blind = FALSE)`); the selected genes are the union of the top-10 absolute loadings of PC1 and PC2 of the training VST matrix (19 genes in the final run; listed in `output/selected_genes.csv`). Test (and CV-validation) samples get their own size factors and are transformed with the training dispersion trend (frozen VST); their labels are not used. Inside each cross-validation fold the genes are re-selected from the fold-training samples only.
- **Feature-selection scenarios** (training data only; the same features for every classifier): S0 all features; SA LASSO (`cv.glmnet`, alpha = 1, 5-fold, lambda.min; fallback top-5 features by absolute correlation); SB SA + mRMR (`mRMRe`, max(#SA, 5) features); SC SB + the features with the largest absolute loading on the principal components covering >= 90% of the variance (at most 10). In the radiogenomic pipeline the 5 radiomic features most correlated with the class are always added to SA (and therefore to SB and SC).
- **Nested cross-validation:** stratified 5-fold outer CV on the training set. In every outer fold, feature selection, hyperparameter selection (stratified inner 20% hold-out, highest inner AUC) and a refit of the selected configuration on the whole fold-training set are performed; the fold-validation patients are used only to compute the CV-AUC. The final model repeats the same procedure on all 80 training patients (selection on an inner 20% hold-out, then refit on all 80) and is evaluated once on the test set.
- **Classifiers:** H2O deep learning (two hidden layers, RectifierWithDropout, grid over layer sizes and input dropout {0, 0.1}, 40 epochs without early stopping, `reproducible = TRUE`, inverse-frequency class weights as observation weights); TensorFlow/Keras (two ReLU dense layers, batch normalisation after the first layer, L1/L2 = 1e-4, dropout grid {0.1, 0.2}, Adam lr 0.001, batch size 8, 40 epochs, inverse-frequency class weights, op determinism enabled and `set_random_seed(123)` before every model); Random Forest (`randomForest`, 500 trees, grid mtry {sqrt(p), p/3} x maxnodes {10, 20, 30}, class weights 1:2 for ADC:SCC). The same hyperparameter grids as in the manuscript are used.
- **Class coding and metrics:** histology is coded 0 = ADC, 1 = SCC; models output P(SCC). **ADC is the positive class** for AUC, sensitivity, specificity, MCC, F1 and the post-hoc analyses. All classifiers use the same fixed threshold (P(SCC) >= 0.5 → SCC). AUC and its 95% CI (DeLong) are computed with `pROC` with a fixed curve direction.
- **Model selection:** Table 2 shows, for each pipeline and scenario, the classifier with the highest CV-AUC; the primary model of each pipeline is the pipeline/scenario/classifier with the highest CV-AUC (ties: fewer features, then H2O, TF, RF). The test set is not used for these choices. All combinations are given in Table S1 (`output/table_S1_all_models.csv`).
- **Post-hoc analyses (primary models):** paired DeLong test matched by patient ID (radiogenomic vs genomic as the primary comparison), Brier score and Cox recalibration, decision curve analysis (event = ADC, thresholds 2-80%).
- **Seeds:** 123 for the split, the folds, the inner hold-outs and every model; the LASSO cross-validation uses seed 123 for the final model and 123 + k inside outer fold k.

## Software

The scripts were last run in the following environment:

| Software | Version |
|---|---|
| R | 4.4.1 (Windows, ucrt) |
| DESeq2 | 1.46.0 |
| glmnet | 4.1.10 |
| mRMRe | 2.1.2.2 |
| h2o (R package and server) | 3.44.0.3 (Java 1.8.0_421) |
| keras3 | 1.5.1 |
| tensorflow (R package) | 2.20.0.9000 (TensorFlow 2.20.0, Python 3.12) |
| randomForest | 4.7.1.2 |
| caret | 7.0.1 |
| pROC | 1.19.0.1 |
| clusterProfiler | 4.14.6 |
| org.Hs.eg.db | 3.20.0 |
| AnnotationDbi | 1.68.0 |
| data.table | 1.17.8 |
| LIFEx (feature extraction) | 7.6.7 |

Neural-network results can differ slightly across hardware and library versions even with fixed seeds. Each model script writes its `sessionInfo()` to `output/`.

## Licence

Code: GNU General Public License v3.0 (GPL-3.0), Copyright (C) 2026 Athina Terzi -- see `LICENSE`. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. The licence of the code does not apply to the TCIA and GEO data, which are governed by their own terms (see above).

## Citation

If you use this code, please cite the manuscript above and the data sources.
