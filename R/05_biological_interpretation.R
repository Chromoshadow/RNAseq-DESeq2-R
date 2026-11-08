# ============================================================
# GSE289159
# 05 - Biological interpretation of the results
# ============================================================

# Execution: run from the project root, after script 04.

# -------------------------------------------------------------------------
# 1. Setup
# -------------------------------------------------------------------------

library(ggplot2)

input_dir <- "results/improved_analysis/tables"
output_dir <- "results/biological_interpretation"
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

stages <- c("hiPSC", "Prog", "Neur")
fdr_threshold <- 0.05
lfc_threshold <- 1
minimum_base_mean <- 10
top_n <- 20


# -------------------------------------------------------------------------
# 2. Load results from 04 script
# -------------------------------------------------------------------------

de_results <- setNames(vector("list", length(stages)), stages)

for (stage in stages) {
  input_file <- file.path(input_dir, paste0("DE_", stage, ".csv"))
  
  if (!file.exists(input_file)) {
    stop(
      "Missing file : ", input_file, ". First run 04_improved_analysis.R."
    )
  }
  
  de_results[[stage]] <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

required_columns <- c(
  "ensembl_id", 
  "gene_symbol", 
  "baseMean",
  "log2FoldChange", 
  "pvalue", 
  "padj"
)

for (stage in stages) {
  missing_columns <- setdiff(required_columns, names(de_results[[stage]]))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in DE_", stage, ".csv : ",
      paste(missing_columns, collapse = ", ")
    )
  }
}


# -------------------------------------------------------------------------
# 3. Significant genes summary
# -------------------------------------------------------------------------

direction_summary <- data.frame(
  stage = stages,
  tested = integer(length(stages)),
  significant = integer(length(stages)),
  higher_in_mutant = integer(length(stages)),
  lower_in_mutant = integer(length(stages)),
  significant_and_large_effect = integer(length(stages))
)

for (i in seq_along(stages)) {
  stage <- stages[i]
  result <- de_results[[stage]]
  
  significant <- !is.na(result$padj) & result$padj < fdr_threshold
  
  direction_summary$tested[i] <- sum(!is.na(result$pvalue))
  direction_summary$significant[i] <- sum(significant)
  direction_summary$higher_in_mutant[i] <- sum(
    significant & result$log2FoldChange > 0,
    na.rm = TRUE
  )
  
  direction_summary$lower_in_mutant[i] <- sum(
    significant & result$log2FoldChange < 0,
    na.rm = TRUE
  )
  
  direction_summary$significant_and_large_effect[i] <- sum(
    significant &
      abs(result$log2FoldChange) >= lfc_threshold &
      result$baseMean >= minimum_base_mean,
    na.rm = TRUE
  )
}

print(direction_summary)

write.csv(
  direction_summary,
  file.path(tables_dir, "DE_direction_summary.csv"),
  row.names = FALSE
)

direction_plot_data <- rbind(
  data.frame(
    stage = direction_summary$stage,
    direction = "Higher in KMT2D_MUT",
    n_genes = direction_summary$higher_in_mutant
  ),
  data.frame(
    stage = direction_summary$stage,
    direction = "Lower in KMT2D_MUT",
    n_genes = direction_summary$lower_in_mutant
  )
)

p_direction <- ggplot(
  direction_plot_data,
  aes(
    x = stage, 
    y = n_genes, 
    fill = direction
  )
) +
  geom_col(position = "dodge") +
  labs(
    title = "Differentially expressed genes by stage",
    subtitle = paste0("FDR < ", fdr_threshold),
    x = "Stage",
    y = "Gene numbers",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(figures_dir, "DE_direction_by_stage.png"),
  p_direction,
  width = 8,
  height = 5,
  dpi = 300
)


# -------------------------------------------------------------------------
# 4. Most significant genes examination
# -------------------------------------------------------------------------

top_genes_by_stage <- list()

for (stage in stages) {
  result <- de_results[[stage]]
  
  result <- result[
    !is.na(result$padj) & !is.na(result$log2FoldChange),
    ,
    drop = FALSE
  ]
  
  result <- result[
    order(result$padj), 
    , 
    drop = FALSE
  ]
  
  result$direction <- ifelse(
    result$log2FoldChange > 0,
    "Higher in KMT2D_MUT",
    ifelse(
      result$log2FoldChange < 0,
      "Lower in KMT2D_MUT",
      "No change"
    )
  )
  
  top_table <- head(
    result[, c(
      "ensembl_id", 
      "gene_symbol", 
      "baseMean",
      "log2FoldChange", 
      "pvalue", 
      "padj", 
      "direction"
    )],
    top_n
  )
  
  top_genes_by_stage[[stage]] <- top_table
  
  write.csv(
    top_table,
    file.path(tables_dir, paste0("top_", top_n, "_genes_", stage, ".csv")),
    row.names = FALSE
  )
}

top_genes_combined <- do.call(
  rbind,
  lapply(names(top_genes_by_stage), function(stage) {
    data.frame(
      stage = stage,
      top_genes_by_stage[[stage]],
      row.names = NULL,
      check.names = FALSE
    )
  })
)

print(top_genes_combined)

write.csv(
  top_genes_combined,
  file.path(tables_dir, "top_genes_all_stages.csv"),
  row.names = FALSE
)


# -------------------------------------------------------------------------
# 5. Compare significant genes between stages
# -------------------------------------------------------------------------

significant_gene_sets <- lapply(de_results, function(result) {
  result$ensembl_id[
    !is.na(result$padj) & result$padj < fdr_threshold
  ]
})

all_gene_ids <- sort(unique(unlist(significant_gene_sets)))

gene_overlap <- data.frame(
  ensembl_id = all_gene_ids,
  stringsAsFactors = FALSE
)

for (stage in stages) {
  gene_overlap[[stage]] <- all_gene_ids %in% significant_gene_sets[[stage]]
}

gene_overlap$n_stages_significant <- rowSums(
  gene_overlap[, stages, drop = FALSE]
)

gene_overlap$pattern <- apply(
  gene_overlap[, stages, drop = FALSE],
  1,
  function(x) paste(stages[as.logical(x)], collapse = " + ")
)

annotation_lookup <- unique(do.call(
  rbind,
  lapply(de_results, function(result) {
    result[, c("ensembl_id", "gene_symbol"), drop = FALSE]
  })
))

annotation_index <- match(
  gene_overlap$ensembl_id,
  annotation_lookup$ensembl_id
)

gene_overlap$gene_symbol <- annotation_lookup$gene_symbol[annotation_index]

write.csv(
  gene_overlap,
  file.path(tables_dir, "significant_gene_overlap.csv"),
  row.names = FALSE
)

overlap_summary <- as.data.frame(
  table(gene_overlap$pattern),
  stringsAsFactors = FALSE
)

names(overlap_summary) <- c("stage_pattern", "n_genes")

overlap_summary <- overlap_summary[
  order(overlap_summary$n_genes, decreasing = TRUE),
  ,
  drop = FALSE
]

print(overlap_summary)

write.csv(
  overlap_summary,
  file.path(tables_dir, "significant_gene_overlap_summary.csv"),
  row.names = FALSE
)

p_overlap <- ggplot(
  overlap_summary,
  aes(
    x = reorder(stage_pattern, n_genes), 
    y = n_genes
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Share of significant genes between stages",
    subtitle = paste0("FDR < ", fdr_threshold),
    x = "Stages at which the gene is significant",
    y = "Gene number"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figures_dir, "significant_gene_overlap.png"),
  p_overlap,
  width = 9,
  height = 6,
  dpi = 300
)


# -------------------------------------------------------------------------
# 6. Export the lists for future enrichment analysis
# -------------------------------------------------------------------------

for (stage in stages) {
  result <- de_results[[stage]]
  
  tested_genes <- result$ensembl_id[!is.na(result$pvalue)]
  significant <- !is.na(result$padj) & result$padj < fdr_threshold
  
  up_genes <- result$ensembl_id[
    significant & result$log2FoldChange > 0
  ]
  
  down_genes <- result$ensembl_id[
    significant & result$log2FoldChange < 0
  ]
  
  write.table(
    unique(tested_genes),
    file.path(tables_dir, paste0("tested_genes_", stage, ".txt")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  
  write.table(
    unique(up_genes),
    file.path(tables_dir, paste0("significant_up_", stage, ".txt")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  
  write.table(
    unique(down_genes),
    file.path(tables_dir, paste0("significant_down_", stage, ".txt")),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

