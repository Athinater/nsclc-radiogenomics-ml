# ==============================================================================
# 11_go_enrichment.R
# Gene Ontology enrichment of the 18 PCA-selected genes, separately for the
# SCC-associated and ADC-associated gene sets.
# clusterProfiler::enrichGO, org.Hs.eg.db, hypergeometric test, Benjamini-Hochberg
# adjusted p < 0.20 (relaxed because the gene sets are small). Writes tables only.
# ==============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
source("scripts/common.R")

# Gene sets: the 18 genes of the manuscript, split by subtype association
# (first 10 = SCC-associated, last 8 = ADC-associated; see manuscript_genes in common.R)
scc_genes <- manuscript_genes[1:10]
adc_genes <- manuscript_genes[11:18]

ontology   <- "BP"     # Biological Process
p_cutoff   <- 0.20     # BH-adjusted p-value cutoff
# Background (universe): clusterProfiler default = all genes with a GO annotation in org.Hs.eg.db.
# To use the expressed genes as background instead, pass universe = rownames(load_counts(counts_file)).

run_go <- function(genes, label) {
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = ontology,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = p_cutoff,
                  qvalueCutoff  = p_cutoff,
                  readable      = FALSE)
  df <- as.data.frame(ego)
  if (nrow(df) == 0) { cat(label, ": no term with adjusted p <", p_cutoff, "\n"); return(df) }
  df <- df[order(df$p.adjust), ]
  df$Group <- label
  df
}

res_scc <- run_go(scc_genes, "SCC")
res_adc <- run_go(adc_genes, "ADC")

for (r in list(res_scc, res_adc)) if (nrow(r) > 0) {
  cat("\n---", r$Group[1], "top terms ---\n")
  print(head(r[, c("Description", "Count", "p.adjust")], 10), row.names = FALSE)
}

all_res <- rbind(res_scc, res_adc)
write.csv(all_res, out_path("go_enrichment_all_terms.csv"), row.names = FALSE)
write.csv(rbind(head(res_scc, 10), head(res_adc, 10)), out_path("go_enrichment_top10.csv"), row.names = FALSE)
cat("\nSaved GO tables to", output_dir, "\n")
cat("Gene sets used: SCC (n =", length(scc_genes), "):", paste(scc_genes, collapse = ", "),
    "\n                ADC (n =", length(adc_genes), "):", paste(adc_genes, collapse = ", "), "\n")
