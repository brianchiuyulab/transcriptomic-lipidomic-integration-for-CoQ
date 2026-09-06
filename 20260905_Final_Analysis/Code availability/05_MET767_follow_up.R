# MET767-associated RNA genes.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c(
  "data.table", "ggplot2", "ComplexHeatmap", "circlize", "ragg", "svglite"
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
  library(ComplexHeatmap)
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
target_feature <- "MET767"
rho_cutoff <- 0.90
FDR_cutoff <- 0.05

correlation_file <- file.path(table_directory, "02_all_key_feature_correlations.csv")
RNA_mean_file <- file.path(table_directory, "01_RNA_key_gene_condition_means.csv")
lipid_mean_file <- file.path(table_directory, "01_lipid_key_feature_condition_means.csv")
matching_file <- file.path(table_directory, "01_matching_design.csv")
input_files <- c(correlation_file, RNA_mean_file, lipid_mean_file, matching_file)
if (any(!file.exists(input_files))) stop("Run Codes 1-4 before Code 5")

target_correlations <- fread(correlation_file)[met_id == target_feature]
target_correlations[, association := fcase(
  spearman_rho >= rho_cutoff & fdr_BH_key_feature_space < FDR_cutoff, "Positive",
  spearman_rho <= -rho_cutoff & fdr_BH_key_feature_space < FDR_cutoff, "Negative",
  default = NA_character_
)]

positive_genes <- target_correlations[association == "Positive"][order(-spearman_rho, symbol)]
negative_genes <- target_correlations[association == "Negative"][order(spearman_rho, symbol)]
positive_genes[, rank_within_direction := frank(-spearman_rho, ties.method = "min")]
negative_genes[, rank_within_direction := frank(spearman_rho, ties.method = "min")]
selected_genes <- rbindlist(list(positive_genes, negative_genes), use.names = TRUE)

selected_table <- selected_genes[, .(
  association,
  rank_within_direction,
  ensembl_gene_id,
  symbol,
  spearman_rho,
  p_value = p_value_unadjusted,
  BH_FDR = fdr_BH_key_feature_space,
  RNA_direction,
  RNA_average_post_log2FC,
  RNA_average_post_FDR = RNA_average_post_fdr
)]
fwrite(
  selected_table,
  file.path(table_directory, "05_MET767_correlated_genes.csv"),
  bom = TRUE
)

RNA_means <- fread(RNA_mean_file)
lipid_means <- fread(lipid_mean_file)
matching_design <- fread(matching_file)
condition_columns <- matching_design$condition

selected_RNA <- RNA_means[match(selected_genes$ensembl_gene_id, ensembl_gene_id)]
RNA_matrix <- as.matrix(selected_RNA[, ..condition_columns])
storage.mode(RNA_matrix) <- "numeric"
rownames(RNA_matrix) <- make.unique(fifelse(
  is.na(selected_RNA$symbol) | selected_RNA$symbol == "",
  selected_RNA$ensembl_gene_id,
  selected_RNA$symbol
))
RNA_z <- t(scale(t(RNA_matrix)))
RNA_z[!is.finite(RNA_z)] <- 0
RNA_z <- pmax(pmin(RNA_z, 2.5), -2.5)

MET767_values <- as.numeric(lipid_means[met_id == target_feature, ..condition_columns])
MET767_z <- as.numeric(scale(MET767_values))

passage_levels <- c("Early Passage", "Middle Passage", "Late Passage")
matching_design[, passage := factor(passage, levels = passage_levels)]
matching_design[, day := factor(day, levels = c("D0", "D1", "D6"))]
column_labels <- paste(
  sub(" Passage$", " P.", as.character(matching_design$passage)),
  matching_design$day,
  sep = " · "
)

passage_colours <- c(
  "Early Passage" = "#3775BA",
  "Middle Passage" = "#E69F00",
  "Late Passage" = "#C44E67"
)
day_colours <- c("D0" = "#F0F0F0", "D1" = "#969696", "D6" = "#252525")
association_colours <- c("Positive" = "#B2182B", "Negative" = "#2166AC")
expression_colour <- circlize::colorRamp2(
  c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")
)
MET767_colour <- circlize::colorRamp2(
  range(c(MET767_z, 0)), c("#2166AC", "#B2182B")
)

make_column_annotation <- function() {
  HeatmapAnnotation(
    Passage = matching_design$passage,
    Day = matching_design$day,
    `MET767 z-score` = MET767_z,
    col = list(
      Passage = passage_colours,
      Day = day_colours,
      `MET767 z-score` = MET767_colour
    ),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    simple_anno_size = grid::unit(3.5, "mm")
  )
}

association_factor <- factor(selected_genes$association, levels = c("Positive", "Negative"))
overview_heatmap <- Heatmap(
  RNA_z,
  name = "RNA z-score",
  col = expression_colour,
  top_annotation = make_column_annotation(),
  left_annotation = rowAnnotation(
    Association = association_factor,
    col = list(Association = association_colours),
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(3.5, "mm")
  ),
  row_split = association_factor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_labels = column_labels,
  column_names_rot = 45,
  column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  row_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  column_title = "MET767-associated RNA genes",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    at = c(-2.5, 0, 2.5),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  use_raster = TRUE,
  raster_quality = 3
)

draw_heatmap <- function(heatmap_object) {
  draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
}

ragg::agg_png(
  file.path(figure_directory, "05_MET767_correlated_gene_set_heatmap.png"),
  width = 7.5, height = 8.5, units = "in", res = 600, background = "white"
)
draw_heatmap(overview_heatmap)
dev.off()
grDevices::cairo_pdf(
  file.path(figure_directory, "05_MET767_correlated_gene_set_heatmap.pdf"),
  width = 7.5, height = 8.5
)
draw_heatmap(overview_heatmap)
dev.off()
svglite::svglite(
  file.path(figure_directory, "05_MET767_correlated_gene_set_heatmap.svg"),
  width = 7.5, height = 8.5, bg = "white"
)
draw_heatmap(overview_heatmap)
dev.off()

make_direction_heatmap <- function(direction) {
  rows <- which(selected_genes$association == direction)
  matrix_subset <- RNA_z[rows, , drop = FALSE]
  label_colours <- rep("black", nrow(matrix_subset))
  if (direction == "Positive") {
    label_colours[rownames(matrix_subset) == "Coq8a"] <- "#D73027"
  }
  Heatmap(
    matrix_subset,
    name = "RNA z-score",
    col = expression_colour,
    top_annotation = make_column_annotation(),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = grid::gpar(
      fontsize = if (nrow(matrix_subset) > 50L) 6.5 else 9,
      col = label_colours,
      fontface = "bold"
    ),
    column_labels = column_labels,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_title = paste0(
      direction, " RNA correlations with MET767 (n = ", nrow(matrix_subset), ")"
    ),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      at = c(-2.5, 0, 2.5),
      title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 10, fontface = "bold")
    ),
    use_raster = TRUE,
    raster_quality = 3
  )
}

save_direction_heatmap <- function(direction) {
  heatmap_object <- make_direction_heatmap(direction)
  gene_count <- sum(selected_genes$association == direction)
  plot_height <- max(5.5, 2.5 + 0.115 * gene_count)
  filename_stub <- paste0("05_MET767_", tolower(direction), "_genes_heatmap")

  ragg::agg_png(
    file.path(figure_directory, paste0(filename_stub, ".png")),
    width = 9.5, height = plot_height, units = "in", res = 600,
    background = "white"
  )
  draw_heatmap(heatmap_object)
  dev.off()
  svglite::svglite(
    file.path(figure_directory, paste0(filename_stub, ".svg")),
    width = 9.5, height = plot_height, bg = "white"
  )
  draw_heatmap(heatmap_object)
  dev.off()
  grDevices::cairo_pdf(
    file.path(figure_directory, paste0(filename_stub, ".pdf")),
    width = 9.5, height = plot_height
  )
  draw_heatmap(heatmap_object)
  dev.off()
}

save_direction_heatmap("Positive")
save_direction_heatmap("Negative")

Coq8a_result <- positive_genes[symbol == "Coq8a"]
all_ranked_genes <- copy(target_correlations)
all_ranked_genes[, rank_by_rho := frank(-spearman_rho, ties.method = "min")]
all_ranked_genes[, status := fcase(
  spearman_rho >= rho_cutoff & fdr_BH_key_feature_space < FDR_cutoff,
  "Significant positive",
  spearman_rho <= -rho_cutoff & fdr_BH_key_feature_space < FDR_cutoff,
  "Significant negative",
  default = "Not selected"
)]
setorder(all_ranked_genes, -spearman_rho, symbol)

Coq8a_position <- which(all_ranked_genes$symbol == "Coq8a")
Coq8a_rank <- all_ranked_genes[symbol == "Coq8a", rank_by_rho]
all_ranked_genes[, plot_order := seq_len(.N)]
all_ranked_genes[, status := factor(
  status,
  levels = c("Significant positive", "Not selected", "Significant negative")
)]
all_ranked_genes[, plot_FDR := pmin(-log10(pmax(
  fdr_BH_key_feature_space, .Machine$double.xmin
)), 12)]

fwrite(
  all_ranked_genes[, .(
    rank_by_rho, symbol, ensembl_gene_id, spearman_rho,
    p_value_unadjusted, fdr_BH_key_feature_space, status
  )],
  file.path(table_directory, "05_MET767_all_gene_correlations_ranked.csv"),
  bom = TRUE
)

Coq8a_plot <- all_ranked_genes[symbol == "Coq8a"]
significance_plot <- ggplot(
  all_ranked_genes,
  aes(spearman_rho, plot_FDR, colour = status)
) +
  geom_point(size = 1.9, alpha = 0.9) +
  geom_vline(
    xintercept = c(-rho_cutoff, rho_cutoff),
    linetype = 2, colour = "grey45", linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(FDR_cutoff),
    linetype = 2, colour = "grey45", linewidth = 0.4
  ) +
  geom_point(
    data = Coq8a_plot,
    shape = 21, size = 4.2, stroke = 0.9,
    fill = "#F2A900", colour = "black"
  ) +
  annotate(
    "segment",
    x = 0.68, y = 2.8,
    xend = Coq8a_plot$spearman_rho, yend = Coq8a_plot$plot_FDR,
    linewidth = 0.45, colour = "black"
  ) +
  annotate(
    "text",
    x = 0.67, y = 2.9,
    hjust = 1, vjust = 0,
    size = 4.2,
    fontface = "bold",
    label = sprintf(
      "Coq8a\nrho = %.3f; BH-FDR = %.4f\nrank %d/%d",
      Coq8a_plot$spearman_rho,
      Coq8a_plot$fdr_BH_key_feature_space,
      Coq8a_plot$rank_by_rho,
      nrow(all_ranked_genes)
    )
  ) +
  scale_colour_manual(
    values = c(
      "Significant positive" = "#B2182B",
      "Not selected" = "#CFCFCF",
      "Significant negative" = "#2166AC"
    ),
    breaks = c("Significant positive", "Not selected", "Significant negative"),
    name = "Status"
  ) +
  scale_x_continuous(
    breaks = c(-1, -0.9, -0.5, 0, 0.5, 0.9, 1),
    limits = c(-1.05, 1.05)
  ) +
  scale_y_continuous(
    breaks = c(0, -log10(FDR_cutoff), 5, 10, 12),
    labels = c("0", "1.30", "5", "10", "12"),
    limits = c(0, 12.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "RNA correlations with MET767 (putative CoQ9)",
    x = "Spearman rho",
    y = expression(-log[10]("BH-FDR")~"(capped at 12)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(figure_directory, "05_MET767_correlation_significance_plot.png"),
  significance_plot, width = 8, height = 5.8, dpi = 600,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  file.path(figure_directory, "05_MET767_correlation_significance_plot.pdf"),
  significance_plot, width = 8, height = 5.8, device = cairo_pdf
)
ggsave(
  file.path(figure_directory, "05_MET767_correlation_significance_plot.svg"),
  significance_plot, width = 8, height = 5.8, device = svglite::svglite,
  bg = "white"
)

ranked_plot <- ggplot(all_ranked_genes, aes(plot_order, spearman_rho)) +
  geom_line(colour = "grey75", linewidth = 0.45) +
  geom_point(aes(colour = status), size = 1.3, alpha = 0.85) +
  geom_hline(
    yintercept = c(-rho_cutoff, rho_cutoff),
    linetype = 2, colour = "grey45", linewidth = 0.4
  ) +
  geom_point(
    data = all_ranked_genes[symbol == "Coq8a"],
    shape = 21, size = 3.8, stroke = 0.8,
    fill = "#D73027", colour = "black"
  ) +
  annotate(
    "text",
    x = Coq8a_position + 35,
    y = Coq8a_result$spearman_rho - 0.08,
    hjust = 0,
    colour = "#D73027",
    fontface = "bold",
    size = 4.2,
    label = sprintf(
      "Coq8a\nrank %d/%d; rho = %.3f; BH-FDR = %.4f",
      Coq8a_rank, nrow(all_ranked_genes),
      Coq8a_result$spearman_rho, Coq8a_result$fdr_BH_key_feature_space
    )
  ) +
  scale_colour_manual(
    values = c(
      "Significant positive" = "#B2182B",
      "Not selected" = "#D3D3D3",
      "Significant negative" = "#2166AC"
    ),
    breaks = c("Significant positive", "Not selected", "Significant negative"),
    name = "Status"
  ) +
  scale_x_continuous(
    breaks = c(1, 500, 1000, 1498),
    limits = c(1, 1498),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(-1, -0.9, -0.5, 0, 0.5, 0.9, 1),
    limits = c(-1.05, 1.05)
  ) +
  labs(
    title = "RNA correlations with MET767 (putative CoQ9)",
    x = "Genes ordered by Spearman rho",
    y = "Spearman rho"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(figure_directory, "05_MET767_all_gene_correlation_rank_plot.png"),
  ranked_plot, width = 8, height = 5.2, dpi = 600,
  device = ragg::agg_png, bg = "white"
)
ggsave(
  file.path(figure_directory, "05_MET767_all_gene_correlation_rank_plot.pdf"),
  ranked_plot, width = 8, height = 5.2, device = cairo_pdf
)
ggsave(
  file.path(figure_directory, "05_MET767_all_gene_correlation_rank_plot.svg"),
  ranked_plot, width = 8, height = 5.2, device = svglite::svglite,
  bg = "white"
)

summary_table <- data.table(
  metric = c(
    "tested_genes", "positive_genes", "negative_genes", "selection_rule",
    "Coq8a_spearman_rho", "Coq8a_BH_FDR", "Coq8a_positive_rank"
  ),
  value = c(
    nrow(target_correlations), nrow(positive_genes), nrow(negative_genes),
    "|Spearman rho| >= 0.90; BH-FDR < 0.05",
    Coq8a_result$spearman_rho,
    Coq8a_result$fdr_BH_key_feature_space,
    Coq8a_result$rank_within_direction
  )
)
fwrite(summary_table, file.path(table_directory, "05_MET767_summary.csv"), bom = TRUE)

cat("Code 5 completed successfully.\n")
cat("Positive genes:", nrow(positive_genes), "\n")
cat("Negative genes:", nrow(negative_genes), "\n")
