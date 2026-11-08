# ============================================================
# GSE289159
# 06 - Biological interpretation of the results
# ============================================================

# -------------------------------------------------------------------------
# 1. Setup
# -------------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

input_directory <- file.path(
  "results",
  "biological_interpretation",
  "tables"
)

output_directory <- file.path(
  "results",
  "functional_enrichment"
)

table_directory <- file.path(output_directory, "tables")
figure_directory <- file.path(output_directory, "figures")

dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)


# -------------------------------------------------------------------------
# 2. Define stages and read gene lists
# -------------------------------------------------------------------------

stages <- c("hiPSC", "Prog", "Neur")

read_gene_list <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("Required input file not found: ", file_path)
  }
  
  gene_ids <- readLines(file_path, warn = FALSE)
  gene_ids <- trimws(gene_ids)
  gene_ids <- gene_ids[gene_ids != ""]
  gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
  
  unique(gene_ids)
}

gene_lists <- list()

input_summary <- data.frame(
  stage = character(),
  list_type = character(),
  number_of_ids = integer(),
  stringsAsFactors = FALSE
)

for (stage in stages) {
  
  tested_file <- file.path(
    input_directory,
    paste0("tested_genes_", stage, ".txt")
  )
  
  up_file <- file.path(
    input_directory,
    paste0("significant_up_", stage, ".txt")
  )
  
  down_file <- file.path(
    input_directory,
    paste0("significant_down_", stage, ".txt")
  )
  
  tested_ids <- read_gene_list(tested_file)
  up_ids <- read_gene_list(up_file)
  down_ids <- read_gene_list(down_file)
  
  gene_lists[[stage]] <- list(
    tested = tested_ids,
    up = up_ids,
    down = down_ids
  )
  
  input_summary <- rbind(
    input_summary,
    data.frame(
      stage = stage,
      list_type = c("tested", "up", "down"),
      number_of_ids = c(
        length(tested_ids),
        length(up_ids),
        length(down_ids)
      ),
      stringsAsFactors = FALSE
    )
  )
  
  up_outside_universe <- setdiff(up_ids, tested_ids)
  down_outside_universe <- setdiff(down_ids, tested_ids)
  
  if (length(up_outside_universe) > 0) {
    warning(
      stage, ": ", length(up_outside_universe),
      " up-regulated ID(s) are absent from the tested-gene universe."
    )
  }
  
  if (length(down_outside_universe) > 0) {
    warning(
      stage, ": ", length(down_outside_universe),
      " down-regulated ID(s) are absent from the tested-gene universe."
    )
  }
}

write.csv(
  input_summary,
  file.path(table_directory, "input_gene_list_summary.csv"),
  row.names = FALSE
)

print(input_summary)


# -------------------------------------------------------------------------
# 3. Set Gene Ontology enrichment parameters
# -------------------------------------------------------------------------

go_ontology <- "BP"
p_value_cutoff <- 0.05
q_value_cutoff <- 0.20
minimum_gene_set_size <- 10
maximum_gene_set_size <- 500


# -------------------------------------------------------------------------
# 4. Define a function to run enrichment for one list
# -------------------------------------------------------------------------

run_go_enrichment <- function(
    gene_ids,
    universe_ids,
    stage_name,
    direction_name
) {
  
  if (length(gene_ids) == 0) {
    message(
      "Skipping ", stage_name, " / ", direction_name,
      ": the input gene list is empty."
    )
    return(NULL)
  }
  
  if (length(universe_ids) == 0) {
    warning(
      "Skipping ", stage_name, " / ", direction_name,
      ": the tested-gene universe is empty."
    )
    return(NULL)
  }
  
  gene_ids <- intersect(gene_ids, universe_ids)
  
  if (length(gene_ids) == 0) {
    message(
      "Skipping ", stage_name, " / ", direction_name,
      ": no input genes remain inside the universe."
    )
    return(NULL)
  }
  
  message(
    "\nRunning GO enrichment for ", stage_name, " / ", direction_name,
    " (", length(gene_ids), " input IDs; ",
    length(universe_ids), " universe IDs)."
  )
  
  result <- tryCatch(
    {
      clusterProfiler::enrichGO(
        gene = gene_ids,
        universe = universe_ids,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = go_ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = p_value_cutoff,
        qvalueCutoff = q_value_cutoff,
        minGSSize = minimum_gene_set_size,
        maxGSSize = maximum_gene_set_size,
        readable = TRUE
      )
    },
    error = function(error_condition) {
      warning(
        "GO enrichment failed for ", stage_name, " / ",
        direction_name, ": ", conditionMessage(error_condition)
      )
      return(NULL)
    }
  )
  
  result
}


# -------------------------------------------------------------------------
# 5. Run six enrichment analyses: 3 stages x 2 directions
# -------------------------------------------------------------------------

enrichment_results <- list()

for (stage in stages) {
  
  current_universe <- gene_lists[[stage]]$tested
  
  for (direction in c("up", "down")) {
    
    current_genes <- gene_lists[[stage]][[direction]]
    result_name <- paste(stage, direction, sep = "_")
    
    enrichment_results[[result_name]] <- run_go_enrichment(
      gene_ids = current_genes,
      universe_ids = current_universe,
      stage_name = stage,
      direction_name = direction
    )
  }
}


# -------------------------------------------------------------------------
# 6. Export the complete enrichment tables
# -------------------------------------------------------------------------


for (result_name in names(enrichment_results)) {
  
  result <- enrichment_results[[result_name]]
  
  output_file <- file.path(
    table_directory,
    paste0("GO_BP_", result_name, ".csv")
  )
  
  if (is.null(result)) {
    
    write.csv(data.frame(), output_file, row.names = FALSE)
    
  } else {
    
    result_table <- as.data.frame(result)
    
    write.csv(
      result_table,
      output_file,
      row.names = FALSE
    )
    
    message(
      result_name, ": ",
      nrow(result_table),
      " GO term(s) returned."
    )
  }
}


# -------------------------------------------------------------------------
# 7. Create dot plots for analyses with returned terms
# -------------------------------------------------------------------------

for (result_name in names(enrichment_results)) {
  
  result <- enrichment_results[[result_name]]
  
  if (is.null(result)) {
    next
  }
  
  result_table <- as.data.frame(result)
  
  if (nrow(result_table) == 0) {
    message(
      "No GO terms to plot for ", result_name,
      "; skipping the figure."
    )
    next
  }
  
  number_of_terms_to_plot <- min(20, nrow(result_table))
  
  enrichment_plot <- enrichplot::dotplot(
    result,
    showCategory = number_of_terms_to_plot,
    x = "GeneRatio",
    color = "p.adjust",
    font.size = 10
  ) +
    ggtitle(
      paste(
        "GO Biological Process enrichment:",
        gsub("_", " / ", result_name)
      )
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  plot_file <- file.path(
    figure_directory,
    paste0("GO_BP_", result_name, ".png")
  )
  
  ggsave(
    filename = plot_file,
    plot = enrichment_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}


# -------------------------------------------------------------------------
# 8. Summarize the number of returned GO terms
# -------------------------------------------------------------------------

enrichment_summary <- data.frame(
  stage_direction = names(enrichment_results),
  number_of_terms = vapply(
    enrichment_results,
    function(result) {
      if (is.null(result)) {
        return(0L)
      }
      
      nrow(as.data.frame(result))
    },
    FUN.VALUE = integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  enrichment_summary,
  file.path(table_directory, "GO_BP_enrichment_summary.csv"),
  row.names = FALSE
)

print(enrichment_summary)


