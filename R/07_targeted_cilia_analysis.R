# ============================================================
# GSE289159
# 07 - Analysis of genes associated with the GO terms of Prog_down
# ============================================================

# ----------------------------------------------------------
# 1. Setup
# ----------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)

go_file <- "results/functional_enrichment/tables/GO_BP_Prog_down.csv"
counts_file <- "results/initial_analysis/tables/count_matrix_filtered.csv"
metadata_file <- "results/sample_description/tables/sample_info.csv"
de_dir <- "results/improved_analysis/tables"

out_dir <- "results/targeted_cilia_analysis"
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# 2. Retrieve the gene symbols from the GO file
# ----------------------------------------------------------

go_results <- read.csv(go_file, stringsAsFactors = FALSE)

gene_symbols <- unlist(
  strsplit(go_results$geneID, split = "/"),
  use.names = FALSE
)

gene_symbols <- sort(unique(
  gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
))

write.table(
  data.frame(gene_symbol = gene_symbols),
  file.path(table_dir, "cilia_gene_symbols.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Distinct symbols in GO terms :", length(gene_symbols), "\n")

# ----------------------------------------------------------
# 3. Read the differential results of the three stages
# ----------------------------------------------------------

de_hiPSC <- read.csv(
  file.path(de_dir, "DE_hiPSC.csv"),
  stringsAsFactors = FALSE
)

de_Prog <- read.csv(
  file.path(de_dir, "DE_Prog.csv"),
  stringsAsFactors = FALSE
)

de_Neur <- read.csv(
  file.path(de_dir, "DE_Neur.csv"),
  stringsAsFactors = FALSE
)

extract_targets <- function(de_table, stage, symbols) {
  selected <- de_table[
    !is.na(de_table$gene_symbol) &
      de_table$gene_symbol %in% symbols,
  ]
  selected$stage <- stage
  selected
}

target_hiPSC <- extract_targets(de_hiPSC, "hiPSC", gene_symbols)
target_Prog <- extract_targets(de_Prog, "Prog", gene_symbols)
target_Neur <- extract_targets(de_Neur, "Neur", gene_symbols)

target_de <- rbind(target_hiPSC, target_Prog, target_Neur)

write.csv(
  target_de,
  file.path(table_dir, "cilia_genes_DE_all_stages.csv"),
  row.names = FALSE
)

de_summary <- do.call(
  rbind,
  lapply(
    list(hiPSC = target_hiPSC, Prog = target_Prog, Neur = target_Neur),
    function(x) {
      data.frame(
        n_genes_found = length(unique(x$gene_symbol)),
        n_negative_log2FC = sum(x$log2FoldChange < 0, na.rm = TRUE),
        n_positive_log2FC = sum(x$log2FoldChange > 0, na.rm = TRUE),
        n_FDR_below_005 = sum(x$padj < 0.05, na.rm = TRUE)
      )
    }
  )
)

de_summary$stage <- rownames(de_summary)
rownames(de_summary) <- NULL

write.csv(
  de_summary,
  file.path(table_dir, "DE_summary_by_stage.csv"),
  row.names = FALSE
)

print(de_summary)

# ----------------------------------------------------------
# 4. Graph of log2 fold changes by stage
# ----------------------------------------------------------

p_fc <- ggplot(
  target_de,
  aes(x = stage, y = log2FoldChange, group = gene_symbol)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.25, color = "grey55") +
  geom_point(aes(shape = padj < 0.05), size = 2) +
  facet_wrap(~ gene_symbol, scales = "free_y") +
  labs(
    title = "Genes associated with the axoneme: variation in expression",
    subtitle = "KMT2D_MUT vs CTRL",
    x = "Stade",
    y = "log2 fold change",
    shape = "FDR < 0.05"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(figure_dir, "cilia_genes_log2FC_by_stage.png"),
  p_fc,
  width = 14,
  height = max(7, 0.3 * length(unique(target_de$gene_symbol))),
  limitsize = FALSE
)

# ----------------------------------------------------------
# 5. Read filtered accounts and metadata
# ----------------------------------------------------------

counts_table <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_ids <- counts_table[[1]]

count_matrix <- as.matrix(counts_table[, -1, drop = FALSE])
rownames(count_matrix) <- gene_ids
storage.mode(count_matrix) <- "integer"

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

metadata$sample_for_counts <- metadata$sample

missing_samples <- setdiff(
  metadata$sample_for_counts,
  colnames(count_matrix)
)

if (length(missing_samples) > 0) {
  stop(
    "Samples missing from count matrix : ",
    paste(missing_samples, collapse = ", ")
  )
}

metadata <- metadata[
  match(colnames(count_matrix), metadata$sample_for_counts),
]

stopifnot(
  all(metadata$sample == colnames(count_matrix))
)


# ----------------------------------------------------------
# 6. Normalize counts for visualization
# ----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(row.names = colnames(count_matrix)),
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)
log_normalized_counts <- log2(normalized_counts + 1)

# ----------------------------------------------------------
# 7. Associate symbols with Ensemble identifiers
# ----------------------------------------------------------

annotation <- unique(rbind(
  de_hiPSC[, c("ensembl_id", "gene_symbol")],
  de_Prog[, c("ensembl_id", "gene_symbol")],
  de_Neur[, c("ensembl_id", "gene_symbol")]
))

annotation <- annotation[
  !is.na(annotation$gene_symbol) &
    annotation$gene_symbol %in% gene_symbols,
]

annotation <- annotation[
  annotation$ensembl_id %in% rownames(log_normalized_counts),
]

write.csv(
  annotation,
  file.path(table_dir, "cilia_gene_id_mapping.csv"),
  row.names = FALSE
)

target_ids <- unique(annotation$ensembl_id)

if (length(target_ids) == 0) {
  stop("No Ensembl identifier corresponding to the GO genes was found.")
}

target_expression <- log_normalized_counts[
  target_ids,
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 8. Heatmap : expression in each sample
# ----------------------------------------------------------

expression_scaled <- t(scale(t(target_expression)))
expression_scaled[is.na(expression_scaled)] <- 0

sample_order <- order(
  metadata$stage,
  metadata$dataset,
  metadata$genotype,
  metadata$replicate
)

expression_scaled <- expression_scaled[, sample_order, drop = FALSE]

sample_annotation <- data.frame(
  stage = metadata$stage[sample_order],
  genotype = metadata$genotype[sample_order],
  dataset = metadata$dataset[sample_order]
)

rownames(sample_annotation) <- colnames(expression_scaled)

gene_labels <- annotation$gene_symbol[
  match(rownames(expression_scaled), annotation$ensembl_id)
]

duplicated_labels <- duplicated(gene_labels) |
  duplicated(gene_labels, fromLast = TRUE)

gene_labels[duplicated_labels] <- paste0(
  gene_labels[duplicated_labels],
  " [",
  rownames(expression_scaled)[duplicated_labels],
  "]"
)

rownames(expression_scaled) <- gene_labels

png(
  file.path(figure_dir, "cilia_genes_expression_heatmap.png"),
  width = 1800,
  height = max(1000, 28 * nrow(expression_scaled)),
  res = 150
)

pheatmap(
  expression_scaled,
  annotation_col = sample_annotation,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 7,
  fontsize_col = 7,
  main = "Expression of genes associated with the axoneme"
)

dev.off()

# ----------------------------------------------------------
# 9. Export the normalized values per sample
# ----------------------------------------------------------

expression_export <- data.frame(
  ensembl_id = rownames(target_expression),
  target_expression,
  check.names = FALSE
)

expression_export <- merge(
  annotation,
  expression_export,
  by = "ensembl_id",
  all.x = TRUE
)

write.csv(
  expression_export,
  file.path(table_dir, "cilia_genes_normalized_expression.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 10. Prog descriptive summary: main cohort and validation
# ----------------------------------------------------------

expression_long <- as.data.frame(
  as.table(target_expression),
  stringsAsFactors = FALSE
)

names(expression_long) <- c(
  "ensembl_id",
  "sample_for_counts",
  "normalized_expression"
)


expression_long <- merge(
  expression_long,
  metadata[, c("sample_for_counts", "stage", "genotype", "dataset")],
  by = "sample_for_counts",
  all.x = TRUE
)


prog_expression <- expression_long[
  expression_long$stage == "Prog",
]

prog_summary <- aggregate(
  normalized_expression ~ ensembl_id + genotype + dataset,
  data = prog_expression,
  FUN = mean
)

prog_summary <- merge(
  annotation,
  prog_summary,
  by = "ensembl_id",
  all.y = TRUE
)

write.csv(
  prog_summary,
  file.path(table_dir, "Prog_main_vs_validation_expression.csv"),
  row.names = FALSE
)

