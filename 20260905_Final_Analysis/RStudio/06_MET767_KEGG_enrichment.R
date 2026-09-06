# KEGG enrichment of MET767-associated RNA genes.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c(
  "data.table", "ggplot2", "AnnotationDbi", "org.Mm.eg.db", "clusterProfiler",
  "ragg", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

working_directory <- getwd()
if (basename(working_directory) == "RStudio") {
  analysis_root <- dirname(working_directory)
} else if (dir.exists(file.path(working_directory, "Tables"))) {
  analysis_root <- working_directory
} else {
  stop("Open OPEN_IN_RSTUDIO.Rproj")
}

table_directory <- file.path(analysis_root, "Tables")
figure_directory <- file.path(analysis_root, "Figure")
gene_file <- file.path(table_directory, "05_MET767_correlated_genes.csv")
universe_file <- file.path(table_directory, "01_RNA_key_gene_condition_means.csv")
if (any(!file.exists(c(gene_file, universe_file)))) stop("Run Code 5 before Code 6")

selected_genes <- fread(gene_file)
universe_genes <- fread(universe_file)[, .(ensembl_gene_id)]
gene_map <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = unique(universe_genes$ensembl_gene_id),
  columns = "ENTREZID",
  keytype = "ENSEMBL"
)
gene_map <- unique(as.data.table(gene_map)[!is.na(ENTREZID)])
setnames(gene_map, "ENSEMBL", "ensembl_gene_id")

mapped_selected <- merge(
  selected_genes, gene_map,
  by = "ensembl_gene_id", all.x = TRUE, allow.cartesian = TRUE
)
universe_ids <- unique(gene_map$ENTREZID)

run_KEGG <- function(direction) {
  gene_ids <- unique(mapped_selected[association == direction & !is.na(ENTREZID), ENTREZID])
  if (length(gene_ids) < 5L) return(data.table())
  result <- clusterProfiler::enrichKEGG(
    gene = gene_ids,
    organism = "mmu",
    keyType = "ncbi-geneid",
    universe = universe_ids,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
  )
  output <- as.data.table(as.data.frame(result))
  if (nrow(output)) output[, association := direction]
  output
}

KEGG_results <- rbindlist(list(
  run_KEGG("Positive"),
  run_KEGG("Negative")
), fill = TRUE)

if (nrow(KEGG_results)) {
  KEGG_results[, Description := sub(" - Mus musculus.*$", "", Description)]
  KEGG_results[, GeneRatio_numeric := vapply(
    strsplit(GeneRatio, "/", fixed = TRUE),
    function(x) as.numeric(x[1]) / as.numeric(x[2]),
    numeric(1)
  )]
  setorder(KEGG_results, association, p.adjust, -Count)
}

significant_KEGG <- KEGG_results[p.adjust < 0.05]
colour_limit <- max(-log10(KEGG_results$p.adjust), na.rm = TRUE)
fwrite(KEGG_results, file.path(table_directory, "06_MET767_KEGG_all_results.csv"), bom = TRUE)
fwrite(
  significant_KEGG,
  file.path(table_directory, "06_MET767_KEGG_significant.csv"),
  bom = TRUE
)

make_KEGG_plot <- function(direction) {
  direction_results <- KEGG_results[association == direction]
  direction_significant <- direction_results[p.adjust < 0.05]
  display_KEGG <- if (nrow(direction_significant)) {
    head(direction_significant, 10)
  } else {
    head(direction_results, 10)
  }
  if (!nrow(display_KEGG)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No KEGG pathways returned") +
        xlim(0, 1) + ylim(0, 1) + theme_void()
    )
  }
  display_KEGG[, pathway_label := factor(Description, levels = rev(Description))]
  subtitle_text <- if (nrow(direction_significant)) {
    paste(nrow(direction_significant), "pathways passed BH-FDR < 0.05")
  } else {
    "No pathways passed BH-FDR < 0.05"
  }
  ggplot(display_KEGG, aes(GeneRatio_numeric, pathway_label)) +
    geom_segment(
      aes(x = 0, xend = GeneRatio_numeric, yend = pathway_label),
      colour = "grey80", linewidth = 0.5
    ) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    scale_colour_gradient(
      low = "#4E79A7", high = "#B2182B",
      limits = c(0, colour_limit)
    ) +
    labs(
      x = "Gene ratio", y = NULL, size = "Genes",
      colour = expression(-log[10]("BH-FDR")),
      title = paste(direction, "MET767-associated genes"),
      subtitle = subtitle_text
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 11, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11, face = "bold")
    )
}

save_KEGG_plot <- function(direction) {
  plot_object <- make_KEGG_plot(direction)
  displayed_count <- min(10L, KEGG_results[association == direction, .N])
  plot_height <- max(3.8, 2.5 + 0.35 * displayed_count)
  filename_stub <- paste0("06_MET767_", tolower(direction), "_KEGG_enrichment")
  ggsave(
    file.path(figure_directory, paste0(filename_stub, ".png")),
    plot_object, width = 7.5, height = plot_height, dpi = 600,
    device = ragg::agg_png, bg = "white"
  )
  ggsave(
    file.path(figure_directory, paste0(filename_stub, ".pdf")),
    plot_object, width = 7.5, height = plot_height, device = cairo_pdf
  )
  ggsave(
    file.path(figure_directory, paste0(filename_stub, ".svg")),
    plot_object, width = 7.5, height = plot_height,
    device = svglite::svglite, bg = "white"
  )
}

save_KEGG_plot("Positive")
save_KEGG_plot("Negative")

summary_table <- data.table(
  metric = c(
    "KEGG_background", "positive_input_genes", "negative_input_genes",
    "positive_significant_pathways", "negative_significant_pathways"
  ),
  value = c(
    length(universe_ids),
    uniqueN(mapped_selected[association == "Positive" & !is.na(ENTREZID), ENTREZID]),
    uniqueN(mapped_selected[association == "Negative" & !is.na(ENTREZID), ENTREZID]),
    significant_KEGG[association == "Positive", .N],
    significant_KEGG[association == "Negative", .N]
  )
)
fwrite(summary_table, file.path(table_directory, "06_MET767_KEGG_summary.csv"), bom = TRUE)

cat("Code 6 completed successfully.\n")
cat("Positive significant pathways:", significant_KEGG[association == "Positive", .N], "\n")
cat("Negative significant pathways:", significant_KEGG[association == "Negative", .N], "\n")
