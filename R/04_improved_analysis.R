# ============================================================
# GSE289159
# 04 - Final analysis after sex QC
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

dir.create(
  "results/improved_analysis/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/improved_analysis/tables",
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

count_matrix <- counts[, -(1:5), drop = FALSE]

if (!all(vapply(count_matrix, is.numeric, logical(1))))
  stop("Some columns are not numeric")
if (anyDuplicated(rownames(count_matrix)))
  stop("Gene identifiers are not unique")


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

print(sample_info)

stopifnot(
  nrow(sample_info) == ncol(count_matrix),
  identical(
    sample_info$sample,
    colnames(count_matrix)
  )
)


# ------------------------------------------------------------
# 4. Load sex QC
# ------------------------------------------------------------

sex_qc <- read.csv(
  "results/sex_qc/tables/sex_qc.csv",
  stringsAsFactors = FALSE
)

print(sex_qc)

sample_info$sex_NCBI <- sex_qc$sex_NCBI[
  match(
    sample_info$sample,
    sex_qc$sample
  )
]

sample_info$sex_inferred <- sex_qc$sex_inferred[
  match(
    sample_info$sample,
    sex_qc$sample
  )
]

stopifnot(
  !anyNA(sample_info$sex_NCBI),
  !anyNA(sample_info$sex_inferred)
)


# ------------------------------------------------------------
# 5. Restrict analysis to main dataset
# ------------------------------------------------------------

main <- sample_info$dataset == "main"

sample_info_main <- droplevels(sample_info[main, , drop = FALSE])

count_matrix_main <- count_matrix[, main, drop = FALSE]

rownames(sample_info_main) <- sample_info_main$sample

stopifnot(identical(sample_info_main$sample, colnames(count_matrix_main)))

write.csv(
  sample_info_main, 
  "results/improved_analysis/tables/sample_metadata_main.csv",
  row.names = FALSE
)

sample_summary <- as.data.frame(
  table(
    stage = sample_info_main$stage, 
    genotype = sample_info_main$genotype)
)

names(sample_summary)[3] <- "n_samples"

write.csv(
  sample_summary, 
  "results/improved_analysis/tables/sample_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Identify sex-chromosome genes
# ------------------------------------------------------------


chr_field <- gsub(
  "\\s+",
  "", 
  as.character(counts$Chr)
)

is_sex_chr <- grepl(
  "(^|;)(chr)?[XY](;|$)",
  chr_field,
  ignore.case = TRUE, 
  perl = TRUE
)

count_autosomal <- count_matrix_main[
  !is_sex_chr, 
  ,
   drop = FALSE
]

# ------------------------------------------------------------
# 7. Filter lowly expressed genes 
# ------------------------------------------------------------

keep <- rowSums(
  count_autosomal >= 10
) >= 3

count_final <- count_autosomal[
  keep, 
  ,
  drop = FALSE
]

filter_summary <- data.frame(
  step = c("Genes in input", "After excluding X/Y", "After low-count filtering"), 
  n_genes = c(nrow(count_matrix_main), nrow(count_autosomal), nrow(count_final))
)

write.csv(
  filter_summary, 
  "results/improved_analysis/tables/filtering_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Build final DESeq2 object
# ------------------------------------------------------------

# Sex is intentionally NOT included in the design.
#
# In the main dataset, genotype and sex are confounded:
# CTRL       -> male
# KMT2D_MUT  -> female
#
# Therefore sex cannot be estimated independently from genotype.


dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_final)),
  colData = sample_info_main,
  design = ~ stage * genotype
)

dds <- DESeq(dds)

write.csv(
  data.frame(
    sample = names(sizeFactors(dds)), 
    size_factor = unname(sizeFactors(dds))),    
  "results/improved_analysis/tables/size_factors.csv", 
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. VST transformation
# ------------------------------------------------------------

vsd <- vst(
  dds, 
  blind = TRUE
)

# ------------------------------------------------------------
# 10. PCA
# ------------------------------------------------------------

vsd_mat <- assay(vsd)

gene_var <- apply(
  vsd_mat, 
  1, 
  var
)

top_n <- min(500, length(gene_var))

top_genes <- names(
  sort(
    gene_var, 
    decreasing = TRUE
  )
)[seq_len(top_n)]

pca <- prcomp(
  t(
    vsd_mat[top_genes, , drop = FALSE]
  ), 
  center = TRUE
)

var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- as.data.frame(
  pca$x[, 1:2, drop = FALSE]
)

pca_data$sample <- rownames(pca_data)


meta_idx <- match(pca_data$sample, sample_info_main$sample)

pca_data$stage <- sample_info_main$stage[meta_idx]
pca_data$genotype <- sample_info_main$genotype[meta_idx]
pca_data$sex_inferred <- sample_info_main$sex_inferred[meta_idx]



stopifnot(
  !anyNA(
    pca_data[, c("stage", "genotype")]
  )
)


# ------------------------------------------------------------
# 11. Save PCA data
# ------------------------------------------------------------

write.csv(
  pca_data, 
  "results/improved_analysis/tables/PCA_coordinates.csv", 
  row.names = FALSE
)


write.csv(
  data.frame(
    component = paste0("PC", seq_along(var_pct)),
    percent_variance = var_pct
  ),
  "results/improved_analysis/tables/PCA_variance.csv", 
  row.names = FALSE
)


p_pca <- ggplot(
  pca_data, 
  aes(
    PC1, 
    PC2, 
    color = genotype, 
    shape = stage
  )
) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCA of main RNA-seq dataset",
    x = sprintf(
      "PC1 (%.1f%%)", 
      var_pct[1]
    ),
       
    y = sprintf(
      "PC2 (%.1f%%)", 
      var_pct[2]
    ),
    color = "Genotype", 
    shape = "Stade"
  ) +
  theme_minimal(base_size = 12)



ggsave(
  "results/improved_analysis/figures/PCA_PC1_PC2.png",
  p_pca,
  width = 8, 
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 12. Gene annotation for differential expression analysis
# ------------------------------------------------------------


gene_annot <- data.frame(
  ensembl_id = rownames(counts),
  Chr = as.character(counts$Chr),
  Length = as.numeric(counts$Length),
  stringsAsFactors = FALSE
)

gene_annot$ensembl_clean <- sub(
  "\\..*$",
  "", 
  gene_annot$ensembl_id
)

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = unique(gene_annot$ensembl_clean),
  keytype = "ENSEMBL", 
  column = "SYMBOL",
  multiVals = "first"
)


gene_annot$gene_symbol <- unname(
  symbol_map[gene_annot$ensembl_clean]
)

annotate_results <- function(res, annotation) {
  df <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  idx <- match(df$ensembl_id, annotation$ensembl_id)
  df$gene_symbol <- annotation$gene_symbol[idx]
  df$Chr <- annotation$Chr[idx]
  df$Length <- annotation$Length[idx]
  df
}


coef_names <- resultsNames(dds)

required <- c(
  "genotype_KMT2D_MUT_vs_CTRL",
  "stageProg.genotypeKMT2D_MUT",
  "stageNeur.genotypeKMT2D_MUT"
)


if (!all(required %in% coef_names))
  stop("Unexpected coefficient names : ", paste(coef_names, collapse = ", "))


# ------------------------------------------------------------
# 13. RNA-seq differential expression analysis
# ------------------------------------------------------------

de_results <- list(
  hiPSC = results(
    dds,
    name = "genotype_KMT2D_MUT_vs_CTRL", 
    alpha = 0.05
  ),
  Prog = results(
    dds, 
    contrast = list(
      c(
        "genotype_KMT2D_MUT_vs_CTRL",
        "stageProg.genotypeKMT2D_MUT"
      )
    ),
    alpha = 0.05
  ),
  Neur = results(
    dds, 
    contrast = list(
      c(
        "genotype_KMT2D_MUT_vs_CTRL",
        "stageNeur.genotypeKMT2D_MUT"
      )
    ), 
    alpha = 0.05
  )
)


# ------------------------------------------------------------
# 14. Result summary
# ------------------------------------------------------------

de_summary <- data.frame(
  stage = names(de_results),
  n_tested = integer(length(de_results)),
  n_FDR_05 = integer(length(de_results)),
  n_FDR_05_absLFC1 = integer(length(de_results))
)

print(de_summary)

for (i in seq_along(de_results)) {
  df <- as.data.frame(de_results[[i]])
  de_summary$n_tested[i] <- sum(!is.na(df$pvalue))
  de_summary$n_FDR_05[i] <- sum(df$padj < 0.05, na.rm = TRUE)
  de_summary$n_FDR_05_absLFC1[i] <- sum(
    df$padj < 0.05 & abs(df$log2FoldChange) >= 1 & df$baseMean >= 10,
    na.rm = TRUE
  )
}


write.csv(
  de_summary,
  "results/improved_analysis/tables/DE_summary.csv",
  row.names = FALSE
)


plot_volcano <- function(df, stage_name) {
  df$significance <- ifelse(
    !is.na(df$padj) & df$padj < 0.05,
    "FDR < 0.05", 
    "Not significant"
  )

  df$neg_log10_p <- -log10(pmax(df$pvalue, 1e-300))

  ggplot(
    df, 
    aes(
      log2FoldChange, 
      neg_log10_p, 
      color = significance
    )
  ) +
    geom_point(
      alpha = 0.65, 
      size = 1.2, 
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = c(-1, 1), 
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05), 
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "Not significant" = "grey70",
        "FDR < 0.05" = "firebrick"
      )
    ) +
    labs(
      title = paste("Differential expression", stage_name),
      subtitle = "KMT2D_MUT vs CTRL",
      x = "log2 fold change", 
      y = "-log10(p-value nominale)", 
      color = NULL
    ) +
    theme_minimal(base_size = 12)
}

fig_dir <- file.path("results/improved_analysis", "figures")
tab_dir <- file.path("results/improved_analysis", "tables")

for (stage_name in names(de_results)) {
  df <- annotate_results(de_results[[stage_name]], gene_annot)

  df <- df[order(df$padj, na.last = TRUE), , drop = FALSE]

  write.csv(
    df, 
    file.path(tab_dir, paste0("DE_", stage_name, ".csv")),
    row.names = FALSE
  )

  ggsave(
    file.path(fig_dir, paste0("volcano_", stage_name, ".png")),
    plot_volcano(df, stage_name), 
    width = 8,
    height = 6, 
    dpi = 300
  )
}

