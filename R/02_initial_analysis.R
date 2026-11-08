# ============================================================
# GSE289159
# 02 - Initial exploratory analysis
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(ggplot2)
library(DESeq2)

dir.create(
  "results/initial_analysis/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/initial_analysis/tables",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Load count matrix
# ------------------------------------------------------------

counts <- read.delim(
  "data/raw/GSE289159_GeneCounts_allSamples.txt.gz",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# First five columns contain genomic annotations
count_matrix <- counts[, -(1:5)]

stopifnot(
  all(sapply(count_matrix, is.numeric))
)


# ------------------------------------------------------------
# 3. Load sample metadata
# ------------------------------------------------------------

sample_info <- read.csv(
  "results/sample_description/tables/sample_info.csv",
  stringsAsFactors = FALSE
)

sample_info$stage <- factor(
  sample_info$stage,
  levels = c("hiPSC", "Prog", "Neur")
)

sample_info$genotype <- factor(
  sample_info$genotype,
  levels = c("CTRL", "KMT2D_MUT")
)

sample_info$dataset <- factor(
  sample_info$dataset,
  levels = c("main", "validation_DS2")
)

stopifnot(
  nrow(sample_info) == ncol(count_matrix),
  identical(sample_info$sample, colnames(count_matrix))
)


# ------------------------------------------------------------
# 4. Filter lowly expressed genes
# ------------------------------------------------------------

# Keep genes with at least 10 counts
# in at least 3 samples

keep <- rowSums(count_matrix >= 10) >= 3

count_matrix_filtered <- count_matrix[keep, ]

filter_summary <- data.frame(
  genes_before = nrow(count_matrix),
  genes_after = nrow(count_matrix_filtered),
  genes_removed = sum(!keep)
)

print(filter_summary)

write.csv(
  count_matrix_filtered,
  "results/initial_analysis/tables/count_matrix_filtered.csv"
)

write.csv(
  filter_summary,
  "results/initial_analysis/tables/filtering_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Build DESeq2 object
# ------------------------------------------------------------

# Initial model:
# stage + genotype + stage:genotype

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = sample_info,
  design = ~ stage * genotype
)


# ------------------------------------------------------------
# 6. Normalize and transform
# ------------------------------------------------------------

dds <- estimateSizeFactors(dds)

cat("\nSize factor summary:\n")
print(summary(sizeFactors(dds)))

vsd <- vst(
  dds,
  blind = TRUE
)


# ------------------------------------------------------------
# 7. PCA
# ------------------------------------------------------------

pca <- plotPCA(
  vsd,
  intgroup = c(
    "stage",
    "genotype",
    "dataset"
  )
)

pca

ggsave(
  "results/initial_analysis/figures/pca_initial.png",
  plot = pca,
  width = 8,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 8. PCA with explicit coordinates
# ------------------------------------------------------------

vsd_mat <- assay(vsd)

gene_variances <- apply(
  vsd_mat,
  1,
  var
)

top_genes <- names(
  sort(
    gene_variances,
    decreasing = TRUE
  )
)[1:500]

pca_res <- prcomp(
  t(vsd_mat[top_genes, ]),
  center = TRUE,
  scale. = FALSE
)

variance_explained <-
  pca_res$sdev^2 /
  sum(pca_res$sdev^2)

pca_data <- as.data.frame(
  pca_res$x[, 1:3]
)

pca_data$sample <- rownames(pca_data)

pca_data$stage <- sample_info$stage
pca_data$genotype <- sample_info$genotype
pca_data$dataset <- sample_info$dataset
pca_data$replicate <- sample_info$replicate

stopifnot(
  !anyNA(pca_data),
  identical(
    pca_data$sample,
    rownames(pca_res$x)
  )
)

write.csv(
  pca_data,
  "results/initial_analysis/tables/pca_coordinates.csv",
  row.names = TRUE
)


# ------------------------------------------------------------
# 9. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)


# ------------------------------------------------------------
# 10. Inspect model coefficients
# ------------------------------------------------------------

print(resultsNames(dds))


# ------------------------------------------------------------
# 11. Differential expression by stage
# ------------------------------------------------------------

res_hiPSC <- results(
  dds,
  name = "genotype_KMT2D_MUT_vs_CTRL",
  alpha = 0.05
)

res_Prog <- results(
  dds,
  contrast = list(
    c(
      "genotype_KMT2D_MUT_vs_CTRL",
      "stageProg.genotypeKMT2D_MUT"
    )
  ),
  alpha = 0.05
)

res_Neur <- results(
  dds,
  contrast = list(
    c(
      "genotype_KMT2D_MUT_vs_CTRL",
      "stageNeur.genotypeKMT2D_MUT"
    )
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 12. Differential expression summary
# ------------------------------------------------------------

de_summary <- data.frame(
  stage = c(
    "hiPSC",
    "Prog",
    "Neur"
  ),
  significant_FDR_05 = c(
    sum(res_hiPSC$padj < 0.05, na.rm = TRUE),
    sum(res_Prog$padj < 0.05, na.rm = TRUE),
    sum(res_Neur$padj < 0.05, na.rm = TRUE)
  )
)

print(de_summary)

write.csv(
  de_summary,
  "results/initial_analysis/tables/DE_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Save full DESeq2 results
# ------------------------------------------------------------

write.csv(
  as.data.frame(res_hiPSC),
  "results/initial_analysis/tables/DE_hiPSC.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(res_Prog),
  "results/initial_analysis/tables/DE_Prog.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(res_Neur),
  "results/initial_analysis/tables/DE_Neur.csv",
  row.names = TRUE
)


# ------------------------------------------------------------
# 14. Interaction analysis
# ------------------------------------------------------------

res_interaction_Prog <- results(
  dds,
  name = "stageProg.genotypeKMT2D_MUT",
  alpha = 0.05
)

res_interaction_Neur <- results(
  dds,
  name = "stageNeur.genotypeKMT2D_MUT",
  alpha = 0.05
)

write.csv(
  as.data.frame(res_interaction_Prog),
  "results/initial_analysis/tables/interaction_Prog.csv",
  row.names = TRUE
)

write.csv(
  as.data.frame(res_interaction_Neur),
  "results/initial_analysis/tables/interaction_Neur.csv",
  row.names = TRUE
)


# ------------------------------------------------------------
# 15. Initial MA plots
# ------------------------------------------------------------

png(
  "results/initial_analysis/figures/MA_hiPSC.png",
  width = 1000,
  height = 800,
  res = 150
)

plotMA(
  res_hiPSC,
  alpha = 0.05,
  main = "Initial DE analysis - hiPSC"
)

dev.off()


png(
  "results/initial_analysis/figures/MA_Prog.png",
  width = 1000,
  height = 800,
  res = 150
)

plotMA(
  res_Prog,
  alpha = 0.05,
  main = "Initial DE analysis - Prog"
)

dev.off()


png(
  "results/initial_analysis/figures/MA_Neur.png",
  width = 1000,
  height = 800,
  res = 150
)

plotMA(
  res_Neur,
  alpha = 0.05,
  main = "Initial DE analysis - Neur"
)

dev.off()


# ------------------------------------------------------------
# 16. Initial DEG filtering
# ------------------------------------------------------------

sig_hiPSC <- res_hiPSC[
  !is.na(res_hiPSC$padj) &
    res_hiPSC$padj < 0.05 &
    abs(res_hiPSC$log2FoldChange) >= 1 &
    res_hiPSC$baseMean >= 10,
]

sig_Prog <- res_Prog[
  !is.na(res_Prog$padj) &
    res_Prog$padj < 0.05 &
    abs(res_Prog$log2FoldChange) >= 1 &
    res_Prog$baseMean >= 10,
]

sig_Neur <- res_Neur[
  !is.na(res_Neur$padj) &
    res_Neur$padj < 0.05 &
    abs(res_Neur$log2FoldChange) >= 1 &
    res_Neur$baseMean >= 10,
]

initial_deg_summary <- data.frame(
  stage = c(
    "hiPSC",
    "Prog",
    "Neur"
  ),
  DEG = c(
    nrow(sig_hiPSC),
    nrow(sig_Prog),
    nrow(sig_Neur)
  )
)

print(initial_deg_summary)

write.csv(
  initial_deg_summary,
  "results/initial_analysis/tables/DEG_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Initial shared DEG analysis
# ------------------------------------------------------------

genes_hiPSC <- rownames(sig_hiPSC)
genes_Prog <- rownames(sig_Prog)
genes_Neur <- rownames(sig_Neur)

common_all <- Reduce(
  intersect,
  list(
    genes_hiPSC,
    genes_Prog,
    genes_Neur
  )
)

cat(
  "\nDEGs shared across all three stages:",
  length(common_all),
  "\n"
)

write.csv(
  data.frame(gene = common_all),
  "results/initial_analysis/tables/common_DEGs.csv",
  row.names = FALSE
)

