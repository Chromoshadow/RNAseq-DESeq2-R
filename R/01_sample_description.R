# ============================================================
# GSE289159
# 01 - Sample description and basic QC
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(ggplot2)

dir.create(
  "results/sample_description/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/sample_description/tables",
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

# Basic dimensions
cat(
  "Count matrix dimensions:",
  nrow(counts), "features x",
  ncol(counts), "columns\n"
)

# Column names
cat("\nColumn names:\n")
print(colnames(counts))

# First rows
cat("\nFirst rows:\n")
print(head(counts[, 1:min(10, ncol(counts))]))


# ------------------------------------------------------------
# 3. Identify annotation and count columns
# ------------------------------------------------------------

annotation_columns <- c(
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

stopifnot(
  identical(
    colnames(counts)[1:5],
    annotation_columns
  )
)

# First five columns contain genomic annotations
count_matrix <- counts[, -(1:5)]

cat(
  "\nCount matrix:",
  nrow(count_matrix), "genes x",
  ncol(count_matrix), "samples\n"
)

# Check that all count columns are numeric
stopifnot(
  all(sapply(count_matrix, is.numeric))
)


# ------------------------------------------------------------
# 4. Describe sample structure
# ------------------------------------------------------------

sample_info <- data.frame(
  sample = colnames(count_matrix)
)

sample_info$stage <- c(
  rep("hiPSC", 6),
  rep("Prog", 12),
  rep("Neur", 5)
)

sample_info$genotype <- c(
  rep("CTRL", 3),
  rep("KMT2D_MUT", 3),
  
  rep("CTRL", 3),
  rep("KMT2D_MUT", 3),
  
  rep("CTRL", 3),
  rep("KMT2D_MUT", 3),
  
  rep("CTRL", 2),
  rep("KMT2D_MUT", 3)
)

sample_info$dataset <- c(
  rep("main", 12),
  rep("validation_DS2", 6),
  rep("main", 5)
)

sample_info$replicate <- c(
  1, 2, 3, 1, 2, 3,
  1, 2, 3, 1, 2, 3,
  1, 2, 3, 1, 2, 3,
  1, 2, 1, 2, 3
)

# Convert categorical variables to factors
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


# ------------------------------------------------------------
# 5. Validate sample metadata
# ------------------------------------------------------------

stopifnot(
  nrow(sample_info) == ncol(count_matrix)
)

stopifnot(
  identical(
    sample_info$sample,
    colnames(count_matrix)
  )
)

stopifnot(
  !anyNA(sample_info)
)

cat("\nSample metadata:\n")
print(sample_info)

cat("\nSamples per stage:\n")
print(table(sample_info$stage))

cat("\nSamples per genotype:\n")
print(table(sample_info$genotype))

cat("\nSamples per dataset:\n")
print(table(sample_info$dataset))

cat("\nSamples per stage and genotype:\n")
print(
  table(
    sample_info$stage,
    sample_info$genotype,
    sample_info$dataset
  )
)


# ------------------------------------------------------------
# 6. Basic library QC
# ------------------------------------------------------------

library_sizes <- colSums(count_matrix)

detected_genes <- colSums(
  count_matrix > 0
)

qc <- data.frame(
  sample = colnames(count_matrix),
  library_size = as.numeric(library_sizes),
  detected_genes = as.numeric(detected_genes)
)

qc$library_size_millions <-
  qc$library_size / 1e6

# Add metadata to QC table
qc$stage <- sample_info$stage
qc$genotype <- sample_info$genotype
qc$dataset <- sample_info$dataset
qc$replicate <- sample_info$replicate


# ------------------------------------------------------------
# 7. Library size summary
# ------------------------------------------------------------

cat("\nLibrary size summary:\n")
print(
  summary(
    qc$library_size_millions
  )
)

cat("\nDetected genes summary:\n")
print(
  summary(
    qc$detected_genes
  )
)


# ------------------------------------------------------------
# 8. Library size plot
# ------------------------------------------------------------

p_library <- ggplot(
  qc,
  aes(
    x = sample,
    y = library_size_millions
  )
) +
  geom_col() +
  labs(
    title = "Library sizes",
    x = "Sample",
    y = "Total counts (millions)"
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
  "results/sample_description/figures/library_sizes.png",
  plot = p_library,
  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 9. Detected genes plot
# ------------------------------------------------------------

p_detected <- ggplot(
  qc,
  aes(
    x = sample,
    y = detected_genes
  )
) +
  geom_col() +
  labs(
    title = "Detected genes per sample",
    x = "Sample",
    y = "Number of detected genes"
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
  "results/sample_description/figures/detected_genes.png",
  plot = p_detected,
  width = 10,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 10. Library complexity
# ------------------------------------------------------------

model <- lm(
  detected_genes ~ library_size_millions,
  data = qc
)

r <- cor(
  qc$library_size_millions,
  qc$detected_genes,
  method = "pearson"
)

cor_test <- cor.test(
  qc$library_size_millions,
  qc$detected_genes,
  method = "pearson"
)

slope <- coef(model)[2]


# ------------------------------------------------------------
# 11. Library complexity plot
# ------------------------------------------------------------

p_complexity <- ggplot(
  qc,
  aes(
    x = library_size_millions,
    y = detected_genes
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = sprintf(
      "r = %.3f\np = %.2e\nslope = %.1f genes / million counts",
      r,
      cor_test$p.value,
      slope
    ),
    hjust = 1.1,
    vjust = -0.5
  ) +
  labs(
    title = "Library size and gene detection",
    x = "Library size (million counts)",
    y = "Detected genes"
  ) +
  theme_classic()

ggsave(
  "results/sample_description/figures/library_size_vs_detected_genes.png",
  plot = p_complexity,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 12. Save QC tables
# ------------------------------------------------------------

write.csv(
  sample_info,
  "results/sample_description/tables/sample_info.csv",
  row.names = FALSE
)

write.csv(
  qc,
  "results/sample_description/tables/qc_metrics.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  metric = c(
    "Number of genes/features",
    "Number of samples",
    "Minimum library size (million counts)",
    "Maximum library size (million counts)",
    "Minimum detected genes",
    "Maximum detected genes",
    "Pearson correlation: library size vs detected genes",
    "Correlation p-value"
  ),
  value = c(
    nrow(count_matrix),
    ncol(count_matrix),
    min(qc$library_size_millions),
    max(qc$library_size_millions),
    min(qc$detected_genes),
    max(qc$detected_genes),
    r,
    cor_test$p.value
  )
)

write.csv(
  qc_summary,
  "results/sample_description/tables/qc_summary.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  "results/sample_description/tables/sessionInfo.txt"
)

