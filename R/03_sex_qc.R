# ============================================================
# GSE289159
# 03 - Sex annotation quality control
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

dir.create(
  "results/sex_qc/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/sex_qc/tables",
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

stopifnot(
  nrow(sample_info) == ncol(count_matrix),
  identical(
    sample_info$sample,
    colnames(count_matrix)
  )
)


# ------------------------------------------------------------
# 4. Sex reported by NCBI
# ------------------------------------------------------------

# These annotations correspond to the manually collected NCBI
# metadata used during the initial analysis.
#
# IMPORTANT:
# This is named sex_NCBI rather than sex because some annotations
# conflict with the published supplementary table.

sample_info$sex_NCBI <- c(
  "M", "M", "M",       
  "F", "F", "F",       
  "M", "M", "M",      
  "F", "F", "F",      
  "M", "M", "M",    
  "F", "F", "F",      
  "M", "M",           
  "F", "F", "F"        
)

stopifnot(
  length(sample_info$sex_NCBI) == nrow(sample_info)
)


# ------------------------------------------------------------
# 5. Prepare Ensembl IDs
# ------------------------------------------------------------

# Remove Ensembl version suffix:
# ENSG00000129824.11 -> ENSG00000129824

gene_ids_clean <- sub(
  "\\..*$",
  "",
  rownames(count_matrix)
)


# ------------------------------------------------------------
# 6. Define chromosome-Y markers
# ------------------------------------------------------------

y_markers <- c(
  "RPS4Y1",
  "DDX3Y",
  "KDM5D",
  "USP9Y",
  "EIF1AY",
  "UTY",
  "ZFY"
)

y_ids <- mapIds(
  org.Hs.eg.db,
  keys = y_markers,
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
)

print(y_ids)


# ------------------------------------------------------------
# 7. Extract chromosome-Y counts
# ------------------------------------------------------------

y_rows <- gene_ids_clean %in% unname(y_ids)

stopifnot(
  sum(y_rows) == length(y_markers)
)

y_counts <- count_matrix[y_rows, ]

# Replace Ensembl IDs by gene symbols for readability
y_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_clean[y_rows],
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

rownames(y_counts) <- unname(y_symbols)


# ------------------------------------------------------------
# 8. Calculate Y-chromosome expression metrics
# ------------------------------------------------------------

library_sizes <- colSums(count_matrix)

y_score <- colSums(y_counts)

y_detected <- colSums(
  y_counts > 0
)

y_cpm <- y_score /
  library_sizes *
  1e6

sex_qc <- data.frame(
  sample = colnames(count_matrix),
  sex_NCBI = sample_info$sex_NCBI,
  Y_counts = as.numeric(y_score),
  Y_genes_detected = as.numeric(y_detected),
  Y_counts_per_million = as.numeric(y_cpm)
)


# ------------------------------------------------------------
# 9. XIST expression
# ------------------------------------------------------------

xist_id <- mapIds(
  org.Hs.eg.db,
  keys = "XIST",
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "first"
)

xist_rows <- gene_ids_clean %in% unname(xist_id)

stopifnot(
  sum(xist_rows) == 1
)

xist_counts <- as.numeric(
  count_matrix[xist_rows, ]
)

xist_cpm <- xist_counts /
  library_sizes *
  1e6

sex_qc$XIST_counts <- xist_counts
sex_qc$XIST_CPM <- xist_cpm


# ------------------------------------------------------------
# 10. Infer sex from RNA-seq evidence
# ------------------------------------------------------------

# This is an exploratory classification, not a biological
# ground truth determination.

sex_qc$sex_inferred <- ifelse(
  sex_qc$Y_genes_detected >= 5 &
    sex_qc$Y_counts_per_million >= 100,
  "likely_M",
  
  ifelse(
    sex_qc$Y_genes_detected == 0 &
      sex_qc$XIST_CPM > 0,
    "likely_F",
    
    "ambiguous"
  )
)


# ------------------------------------------------------------
# 11. Display QC table
# ------------------------------------------------------------

sex_qc <- sex_qc[
  order(
    sex_qc$Y_counts_per_million,
    decreasing = TRUE
  ),
]

print(sex_qc)


# ------------------------------------------------------------
# 12. Summary tables
# ------------------------------------------------------------

cat("\nSex reported by NCBI:\n")
print(
  table(sample_info$sex_NCBI)
)

cat("\nRNA-seq inferred sex:\n")
print(
  table(
    sex_qc$sex_inferred
  )
)

cat("\nNCBI annotation vs RNA-seq evidence:\n")
print(
  table(
    sex_qc$sex_NCBI,
    sex_qc$sex_inferred
  )
)


# ------------------------------------------------------------
# 13. Y-chromosome marker matrix
# ------------------------------------------------------------

y_marker_counts <- data.frame(
  gene_symbol = rownames(y_counts),
  y_counts,
  check.names = FALSE
)

write.csv(
  y_marker_counts,
  "results/sex_qc/tables/Y_marker_counts.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Sex QC table
# ------------------------------------------------------------

write.csv(
  sex_qc,
  "results/sex_qc/tables/sex_qc.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Y expression plot
# ------------------------------------------------------------

p_y <- ggplot(
  sex_qc,
  aes(
    x = reorder(sample, Y_counts_per_million),
    y = Y_counts_per_million,
    fill = sex_NCBI
  )
) +
  geom_col() +
  labs(
    title = "Chromosome-Y expression by sample",
    x = "Sample",
    y = "Y-marker counts per million",
    fill = "Sex reported by NCBI"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(
  "results/sex_qc/figures/Y_expression.png",
  plot = p_y,
  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 16. XIST expression plot
# ------------------------------------------------------------

p_xist <- ggplot(
  sex_qc,
  aes(
    x = reorder(sample, XIST_CPM),
    y = XIST_CPM,
    fill = sex_NCBI
  )
) +
  geom_col() +
  labs(
    title = "XIST expression by sample",
    x = "Sample",
    y = "XIST counts per million",
    fill = "Sex reported by NCBI"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave(
  "results/sex_qc/figures/XIST_expression.png",
  plot = p_xist,
  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 17. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  "results/sex_qc/tables/sessionInfo.txt"
)

