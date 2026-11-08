# ============================================================
# GSE289159
# 8. Targeted screening of genes associated with ciliary motility
# ============================================================

# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create(
  "results/concluding_analysis/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/concluding_analysis/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

counts <- read.csv(
  "results/initial_analysis/tables/count_matrix_filtered.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  "results/sample_description/tables/sample_info.csv",
  stringsAsFactors = FALSE
)

go_results <- read.csv(
  "results/functional_enrichment/tables/GO_BP_Prog_down.csv",
  stringsAsFactors = FALSE
)


de_files <- c(
  "results/improved_analysis/tables/DE_hiPSC.csv",
  "results/improved_analysis/tables/DE_Prog.csv",
  "results/improved_analysis/tables/DE_Neur.csv"
)

de_annotation <- bind_rows(
  lapply(de_files, read.csv)
) %>%
  select(ensembl_id, gene_symbol) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  distinct()

# ------------------------------------------------------------
# 2. Identify the 28 genes
# ------------------------------------------------------------

target_genes <- unique(
  unlist(strsplit(go_results$geneID, "/"))
)

cat(
  "Number of targeted genes :",
  length(target_genes),
  "\n"
)

target_annotation <- de_annotation %>%
  filter(gene_symbol %in% target_genes) %>%
  distinct(ensembl_id, gene_symbol)

# ------------------------------------------------------------
# 3. Construct the count matrix
# ------------------------------------------------------------

rownames(counts) <- counts[[1]]

count_matrix <- as.matrix(
  counts[, -1, drop = FALSE]
)

storage.mode(count_matrix) <- "integer"

metadata$sample_for_counts <- metadata$sample

missing_samples <- setdiff(
  colnames(count_matrix),
  metadata$sample_for_counts
)

if (length(missing_samples) > 0) {
  stop(
    "Samples missing from metadata : ",
    paste(missing_samples, collapse = ", ")
  )
}

metadata <- metadata[
  match(colnames(count_matrix), metadata$sample_for_counts),
]

stopifnot(
  all(metadata$sample == colnames(count_matrix))
)

rownames(metadata) <- metadata$sample

# ------------------------------------------------------------
# 4. Normalize counts
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

# ------------------------------------------------------------
# 5. Prepare the expression data for the targeted genes
# ------------------------------------------------------------

target_expression <- normalized_counts[
  rownames(normalized_counts) %in%
    target_annotation$ensembl_id,
  ,
  drop = FALSE
]

expression_df <- as.data.frame(target_expression) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  pivot_longer(
    cols = -ensembl_id,
    names_to = "sample",
    values_to = "normalized_count"
  ) %>%
  left_join(
    target_annotation,
    by = "ensembl_id"
  ) %>%
  left_join(
    metadata %>%
      select(sample, stage, genotype, dataset),
    by = "sample"
  ) %>%
  mutate(
    log2_expression = log2(normalized_count + 1),
    stage = factor(stage, levels = c("hiPSC", "Prog", "Neur")),
    genotype = factor(
      genotype,
      levels = c("CTRL", "KMT2D_MUT")
    )
  )


cat(
  "Number of symbols found in data :",
  n_distinct(expression_df$gene_symbol),
  "\n"
)

# ------------------------------------------------------------
# 6. Visualize the 28 genes
# ------------------------------------------------------------

p_all_stages <- ggplot(
  expression_df,
  aes(
    x = stage,
    y = log2_expression,
    colour = genotype,
    shape = dataset
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0
    ),
    size = 2,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ gene_symbol,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    title = "Expression of genes associated with ciliary motility",
    subtitle = "Normalized counts, log2 scale",
    x = "Differentiation stage",
    y = "log2(normalized counts + 1)",
    colour = "Genotype",
    shape = "Dataset"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  "results/concluding_analysis/figures/cilia_genes_all_stages_datasets.png",
  p_all_stages,
  width = 14,
  height = 18,
  dpi = 300
)


# ------------------------------------------------------------
# 7. Extract only Prog samples from the main dataset
# ------------------------------------------------------------

prog_main <- expression_df %>%
  filter(
    stage == "Prog",
    dataset == "main"
  )

print(
  prog_main %>%
    distinct(sample, genotype) %>%
    count(genotype)
)

# ------------------------------------------------------------
# 8. Summarize the expression by gene and genotype
# ------------------------------------------------------------

prog_main_summary <- prog_main %>%
  group_by(gene_symbol, genotype) %>%
  summarise(
    n_samples = n(),
    mean_expression = mean(normalized_count),
    median_expression = median(normalized_count),
    sd_expression = sd(normalized_count),
    mean_log2_expression = mean(log2_expression),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = genotype,
    values_from = c(
      n_samples,
      mean_expression,
      median_expression,
      sd_expression,
      mean_log2_expression
    )
  ) %>%
  mutate(
    descriptive_log2FC = log2(
      (mean_expression_KMT2D_MUT + 1) /
      (mean_expression_CTRL + 1)
    )
  ) %>%
  arrange(descriptive_log2FC)

write.csv(
  prog_main_summary,
  "results/concluding_analysis/tables/Prog_main_cilia_descriptive_comparison.csv",
  row.names = FALSE
)

print(prog_main_summary)

# ------------------------------------------------------------
# 9. View the Prog main comparison
# ------------------------------------------------------------

p_prog_main <- ggplot(
  prog_main,
  aes(
    x = genotype,
    y = log2_expression,
    colour = genotype
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0
    ),
    size = 2.2
  ) +
  stat_summary(
    aes(group = genotype),
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    colour = "black"
  ) +
  facet_wrap(
    ~ gene_symbol,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    title = "Expression of genes associated with ciliary motility",
    subtitle = "Prog — main dataset only",
    x = "Genotype",
    y = "log2(normalized counts + 1)",
    colour = "Genotype"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  "results/concluding_analysis/figures/cilia_genes_Prog_main_only.png",
  p_prog_main,
  width = 14,
  height = 18,
  dpi = 300
)

# ----------------------------------------------------------
# 10. Article comparison
# ----------------------------------------------------------


article_genes <- tribble(
  ~gene_symbol, ~article_direction,
  "FEZ1",   "Down",
  "VANGL2", "Down",
  "CDK10",  "Down",
  "DNALI1", "Down",
  "AJAP1",  "Up",
  "IRX2",   "Up",
  "TBX1",   "Up",
  "EN2",    "Up"
)


de_prog <- read.csv(
  "results/improved_analysis/tables/DE_Prog.csv"
)


prog_summary <- read.csv(
  "results/concluding_analysis/tables/Prog_main_cilia_descriptive_comparison.csv"
)


article_comparison <- article_genes %>%
  left_join(
    de_prog %>%
      select(
        gene_symbol,
        log2FoldChange,
        padj
      ),
    by = "gene_symbol"
  ) %>%
  left_join(
    prog_summary %>%
      select(
        gene_symbol,
        descriptive_log2FC
      ),
    by = "gene_symbol"
  ) %>%
  mutate(
    our_direction = case_when(
      is.na(log2FoldChange) ~ "Not found",
      log2FoldChange < 0 ~ "Down",
      log2FoldChange > 0 ~ "Up",
      TRUE ~ "Unchanged"
    ),
    direction_concordance =
      article_direction == our_direction,
    significant_FDR05 =
      !is.na(padj) & padj < 0.05
  ) %>%
  select(
    gene_symbol,
    article_direction,
    our_direction,
    direction_concordance,
    log2FoldChange,
    padj,
    significant_FDR05,
    descriptive_log2FC
  )

print(article_comparison, n = Inf)
