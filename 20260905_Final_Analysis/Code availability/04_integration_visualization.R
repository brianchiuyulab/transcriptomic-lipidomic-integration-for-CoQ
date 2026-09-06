# Global transcript-lipid correlation heatmaps.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "ComplexHeatmap", "circlize")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
})

working_directory <- getwd()
if (basename(working_directory) %in% c("C2C12_Integration analysis", "RStudio")) {
  analysis_root <- dirname(working_directory)
} else if (
  dir.exists(file.path(working_directory, "Tables")) &&
    dir.exists(file.path(working_directory, "Objects"))
) {
  analysis_root <- working_directory
} else {
  stop("Open the OPEN_IN_RSTUDIO RStudio project")
}

table_directory <- file.path(analysis_root, "Tables")
figure_directory <- file.path(analysis_root, "Figure")
object_directory <- file.path(analysis_root, "Objects")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

correlation_object <- file.path(object_directory, "02_key_feature_correlation_matrices.rds")
high_pair_file <- file.path(table_directory, "02_high_correlation_pairs.csv")
class_file <- file.path(table_directory, "03_lipid_annotations.csv")
rna_file <- file.path(table_directory, "01_RNA_key_gene_condition_means.csv")
if (any(!file.exists(c(correlation_object, high_pair_file, class_file, rna_file)))) {
  stop("Run Codes 1-3 before visualization")
}

correlation <- readRDS(correlation_object)
rho <- correlation$rho
correlation_fdr <- correlation$fdr_BH_key_feature_space
high_pairs <- fread(high_pair_file)
feature_classes <- fread(class_file)
rna_features <- fread(rna_file)

feature_classes[, plot_class := fcase(
  annotation_level == "Putative QI annotation" & lipid_class == "Phosphatidylcholine", "PC",
  annotation_level == "Putative QI annotation" & lipid_class == "Phosphatidylethanolamine", "PE",
  annotation_level == "Putative QI annotation" & lipid_class == "Phosphatidylglycerol", "PG",
  annotation_level == "Putative QI annotation" & lipid_class == "Lysophosphatidylcholine", "LPC",
  annotation_level == "Putative QI annotation" & lipid_class == "Sphingomyelin", "SM",
  annotation_level == "Putative QI annotation" & lipid_class == "Free fatty acid", "FFA",
  annotation_level == "Putative LMSD candidate from QI-inferred neutral mass" & lipid_class == "Ubiquinone", "Ubiquinone",
  annotation_level == "Putative LMSD candidate from m/z-adduct search" & lipid_class == "Glycerolipid candidate", "Glycerolipid candidate",
  annotation_level == "Putative LMSD candidate from m/z-adduct search" & lipid_class == "Sphingolipid candidate", "Sphingolipid candidate",
  default = "Unassigned"
)]

class_colours <- c(
  "PC" = "#4E79A7", "PE" = "#F28E2B", "PG" = "#59A14F",
  "LPC" = "#76B7B2", "SM" = "#E15759", "FFA" = "#EDC948",
  "Ubiquinone" = "#6A3D9A", "Glycerolipid candidate" = "#9C755F",
  "Sphingolipid candidate" = "#1B9E77",
  "Unassigned" = "#D3D3D3"
)
class_order <- c(
  "FFA", "LPC", "PC", "PE", "PG", "SM",
  "Glycerolipid candidate", "Sphingolipid candidate", "Ubiquinone", "Unassigned"
)
direction_colours <- c(Higher = "#B2182B", Lower = "#2166AC")
rho_colours <- circlize::colorRamp2(
  c(-1, -0.90, 0, 0.90, 1),
  c("#053061", "#2166AC", "white", "#B2182B", "#67001F")
)

significant_strong <- abs(rho) >= 0.90 & correlation_fdr < 0.05
gene_degree <- rowSums(significant_strong)
feature_degree <- colSums(significant_strong)

gene_labels <- rna_features[, .(
  ensembl_gene_id, symbol, RNA_direction = trajectory_direction
)]
gene_direction <- setNames(
  c(higher = "Higher", lower = "Lower")[gene_labels$RNA_direction],
  gene_labels$ensembl_gene_id
)

all_feature_annotation <- feature_classes[match(colnames(rho), met_id)]
all_feature_annotation[, trajectory_direction := c(
  higher = "Higher", lower = "Lower"
)[trajectory_direction]]
all_feature_annotation[, plot_class := factor(plot_class, levels = class_order)]
all_display_matrix <- rho
all_display_matrix[!significant_strong] <- NA_real_
all_column_order <- order(
  all_feature_annotation$plot_class,
  -feature_degree[all_feature_annotation$met_id],
  all_feature_annotation$met_id
)
all_feature_annotation <- all_feature_annotation[all_column_order]
all_display_matrix <- all_display_matrix[, all_feature_annotation$met_id, drop = FALSE]

all_heatmap_object <- Heatmap(
  all_display_matrix,
  name = "Spearman\nrho",
  col = rho_colours,
  na_col = "white",
  top_annotation = HeatmapAnnotation(
    `Lipid class` = all_feature_annotation$plot_class,
    `Lipid direction` = all_feature_annotation$trajectory_direction,
    col = list(
      `Lipid class` = class_colours,
      `Lipid direction` = direction_colours
    ),
    simple_anno_size = grid::unit(3.5, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  left_annotation = rowAnnotation(
    `RNA direction` = gene_direction[rownames(all_display_matrix)],
    col = list(`RNA direction` = direction_colours),
    simple_anno_size = grid::unit(3.5, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  cluster_rows = hclust(dist(rho)),
  cluster_columns = FALSE,
  column_split = all_feature_annotation$plot_class,
  cluster_column_slices = FALSE,
  column_labels = all_feature_annotation$display_label,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 7, fontface = "bold"),
  column_names_rot = 90,
  column_title = "All selected RNA and lipid features grouped by lipid class",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    at = c(-1, -0.9, 0, 0.9, 1),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
)

pdf(
  file.path(figure_directory, "04_global_integration_heatmap_all_key_features.pdf"),
  width = 16, height = 12, useDingbats = FALSE
)
draw(all_heatmap_object, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()

annotated_feature_annotation <- feature_classes[
  reportable_putative == TRUE & plot_class != "Unassigned" & met_id %in% colnames(rho)
]
annotated_feature_ids <- annotated_feature_annotation$met_id
annotated_full_matrix <- rho[, annotated_feature_ids, drop = FALSE]
annotated_fdr_matrix <- correlation_fdr[, annotated_feature_ids, drop = FALSE]
annotated_display_matrix <- annotated_full_matrix
annotated_display_matrix[
  abs(annotated_display_matrix) < 0.90 | annotated_fdr_matrix >= 0.05
] <- NA_real_
annotated_feature_annotation <- feature_classes[match(
  colnames(annotated_display_matrix), met_id
)]
annotated_feature_annotation[, trajectory_direction := c(
  higher = "Higher", lower = "Lower"
)[trajectory_direction]]
annotated_class_order <- setdiff(class_order, "Unassigned")
annotated_feature_annotation[, plot_class := factor(
  plot_class, levels = annotated_class_order
)]
annotated_column_order <- order(
  annotated_feature_annotation$plot_class,
  -feature_degree[annotated_feature_annotation$met_id],
  annotated_feature_annotation$met_id
)
annotated_feature_annotation <- annotated_feature_annotation[annotated_column_order]
annotated_full_matrix <- annotated_full_matrix[, annotated_feature_annotation$met_id, drop = FALSE]
annotated_display_matrix <- annotated_display_matrix[, annotated_feature_annotation$met_id, drop = FALSE]

annotated_heatmap_object <- Heatmap(
  annotated_display_matrix,
  name = "Spearman\nrho",
  col = rho_colours,
  na_col = "white",
  top_annotation = HeatmapAnnotation(
    `Lipid class` = annotated_feature_annotation$plot_class,
    `Lipid direction` = annotated_feature_annotation$trajectory_direction,
    col = list(
      `Lipid class` = class_colours[annotated_class_order],
      `Lipid direction` = direction_colours
    ),
    simple_anno_size = grid::unit(3.5, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  left_annotation = rowAnnotation(
    `RNA direction` = gene_direction[rownames(annotated_display_matrix)],
    col = list(`RNA direction` = direction_colours),
    simple_anno_size = grid::unit(3.5, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  ),
  cluster_rows = hclust(dist(annotated_full_matrix)),
  cluster_columns = FALSE,
  column_split = annotated_feature_annotation$plot_class,
  cluster_column_slices = FALSE,
  column_labels = annotated_feature_annotation$display_label,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
  column_names_rot = 90,
  column_title = "RNA correlations with lipid-class assigned putative features",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    at = c(-1, -0.9, 0, 0.9, 1),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
)

pdf(
  file.path(figure_directory, "04_global_integration_heatmap_annotated_lipids_only.pdf"),
  width = 12, height = 12, useDingbats = FALSE
)
draw(annotated_heatmap_object, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()
png(
  file.path(figure_directory, "04_global_integration_heatmap_annotated_lipids_only.png"),
  width = 3600, height = 3600, res = 300, bg = "white"
)
draw(annotated_heatmap_object, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()
png(
  file.path(figure_directory, "04_global_integration_heatmap_all_key_features.png"),
  width = 4800, height = 3600, res = 300, bg = "white"
)
draw(all_heatmap_object, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()
visualization_summary <- data.table(
  metric = c(
    "tested_gene_lipid_pairs", "significant_high_correlation_pairs",
    "all_heatmap_genes", "all_heatmap_lipid_features",
    "annotated_heatmap_genes", "annotated_heatmap_lipid_features"
  ),
  value = c(
    length(rho), nrow(high_pairs),
    nrow(rho), ncol(rho),
    nrow(annotated_display_matrix), ncol(annotated_display_matrix)
  )
)
fwrite(visualization_summary, file.path(table_directory, "04_visualization_summary.csv"), bom = TRUE)

cat("Code 4 completed successfully.\n")
cat("All-feature heatmap:", nrow(rho), "genes x", ncol(rho), "lipid features\n")
cat(
  "Annotated-lipid heatmap:", nrow(annotated_display_matrix), "genes x",
  ncol(annotated_display_matrix), "lipid features\n"
)
