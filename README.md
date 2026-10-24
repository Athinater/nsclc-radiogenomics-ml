# Radiomic, genomic and radiogenomic machine-learning models for NSCLC subtype classification

Code accompanying the manuscript:

> Terzi A, Dikaios N, Protonotarios NE, Athanasiadis E, Kastis GA.
> *Evaluation of Radiomic, Genomic, and Radiogenomic Machine Learning Models for Non-Small Cell Lung Cancer Subtype Classification.*
> Journal of Medical Imaging. [citation / DOI to be added]

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
| `clinicaldatawithout.csv` | Radiomic features + `Histology`, one row per patient (semicolon-separated, patient ID as first column) |
| `clinicaldata.csv` | Same structure, restricted to patients with matched RNA-Seq |
| `onlyclinical.csv` | Clinical metadata incl. `Histology` ("Adenocarcinoma" / "Squamous cell carcinoma"), used as DESeq2 `colData` |

File names are set at the top of `scripts/common.R`. The final cohort has 113 patients (87 ADC, 26 SCC); selection criteria are given in the manuscript.

The genomic scripts (`05`, `06`) need only the GEO counts and `onlyclinical.csv`/`clinicaldata.csv`.

## Repository structure

```
scripts/
  common.R                      paths, helper functions, shared gene-selection code
  01_radiomic_h2o.R             radiomic pipeline, H2O deep learning
  02_radiomic_tf.R              radiomic pipeline, TensorFlow/Keras
  03_radiogenomic_h2o.R         radiogenomic pipeline, H2O deep learning
  04_radiogenomic_tf.R          radiogenomic pipeline, TensorFlow/Keras
  05_genomic_h2o.R              genomic pipeline, H2O deep learning
  06_genomic_tf.R               genomic pipeline, TensorFlow/Keras
  07_random_forest.R            Random Forest, all three pipelines
  analysis_helpers.R            helpers for scripts 08-10
  08_delong_comparison.R        paired DeLong test, radiogenomic vs genomic
  09_calibration_brier.R        Brier score, Cox recalibration intercept/slope
  10_decision_curve.R           decision curve analysis (net benefit)
  11_go_enrichment.R            GO enrichment of the 18 selected genes
  12_radiomic_gene_correlations.R   Spearman correlations, radiomic features vs genes
data/                           input files (not distributed; see above)
output/                         created by the scripts
```

## How to run

Run every script from the repository root (the folder containing `scripts/` and `data/`):

```
Rscript scripts/01_radiomic_h2o.R
```

Scripts `01`-`07` are independent training/evaluation scripts. Each writes, to `output/`, a results table (`<run>_results.csv`, with AUC, sensitivity, specificity, MCC, F1 and confusion-matrix counts) and test-set predictions (`<run>_preds_<scenario>.csv`). Scripts `08`-`10` read these prediction files, so run `01`-`07` first. Scripts `11` and `12` need the data in `data/` but not the predictions. The scripts write tables only; no figures are produced.

## Methods implemented

- **Split:** stratified 70/30 split (training n = 80, test n = 33; test set 26 ADC, 7 SCC), seed 123. The test set is used once, for final evaluation.
- **Scenarios** (all feature selection uses training data only): S0 all features; SA LASSO (`cv.glmnet`, alpha = 1, 5-fold, lambda.min); SB LASSO + mRMR; SC SB + principal-component loading features covering >= 90% of the variance (capped at 10 components).
- **Genes:** 18 genes from the variance-stabilised training counts (DESeq2, `vst(blind = FALSE)`): top-10 absolute loadings of PC1 and PC2. Inside each cross-validation fold the genes are re-selected from the fold-training samples only. Validation/test counts are transformed with the training dispersions and their own size factors.
- **Evaluation:** nested 5-fold cross-validation with an inner stratified 20% hold-out for hyperparameter selection, then a final refit on the full training set and one evaluation on the test set.
- **Classifiers:** H2O deep learning (two hidden layers, grid search, 40 epochs); TensorFlow/Keras (two ReLU dense layers with batch normalisation, L1/L2 regularisation, Adam, batch size 8, 40 epochs, inverse-frequency class weights); Random Forest (`randomForest`, 500 trees, `maxnodes` grid, class weights 1:2 for ADC:SCC). H2O models do not use class weights (`balance_classes = FALSE`).
- **Class coding and metrics:** histology is coded 0 = ADC, 1 = SCC. Sensitivity, specificity, MCC and F1 use **ADC as the positive class**. TensorFlow and Random Forest predictions use a fixed 0.5 probability threshold; H2O uses its default predicted label. AUC and its confidence interval are computed with `pROC`.
- **Seeds:** 123 throughout (`set.seed(123)`, `tf$random$set_seed(123)`, `seed = 123` in H2O).

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

Neural-network results can differ slightly across hardware and library versions even with fixed seeds. Please record `sessionInfo()` when reproducing.

## Licence

Code: GNU General Public License v3.0 (GPL-3.0), Copyright (C) 2026 Athina Terzi -- see `LICENSE`. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. The licence of the code does not apply to the TCIA and GEO data, which are governed by their own terms (see above).

## Citation

If you use this code, please cite the manuscript above and the data sources.
