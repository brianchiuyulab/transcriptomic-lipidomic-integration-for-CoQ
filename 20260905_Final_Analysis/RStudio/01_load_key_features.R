# Load independently selected RNA and lipid features and create exact-day condition means.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

required_packages <- c("data.table", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(data.table))

find_script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) return(sub("^--file=", "", file_argument[1]))
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(path)) return(path)
  }
  ""
}

working_directory <- getwd()
script_path <- find_script_path()
if (basename(working_directory) %in% c("C2C12_Integration analysis", "RStudio")) {
  analysis_root <- dirname(working_directory)
} else if (
  dir.exists(file.path(working_directory, "Tables")) &&
    dir.exists(file.path(working_directory, "Objects"))
) {
  analysis_root <- working_directory
} else if (nzchar(script_path)) {
  analysis_root <- dirname(dirname(script_path))
} else {
  stop("Open the OPEN_IN_RSTUDIO RStudio project")
}

setwd(analysis_root)
table_directory <- file.path(analysis_root, "Tables")
object_directory <- file.path(analysis_root, "Objects")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(object_directory, recursive = TRUE, showWarnings = FALSE)

data_root <- dirname(dirname(analysis_root))
rna_root <- file.path(data_root, "RNA-seq", "20260801_Final_Analysis")
lipid_root <- file.path(data_root, "LC-MS", "20260905_Final_Analysis")

rna_selected_file <- file.path(
  rna_root, "Tables", "02_selected_consistent_direction_genes.csv"
)
rna_expression_file <- file.path(dirname(rna_root), "DESeq2_Chiu.xlsx")
lipid_selected_file <- file.path(
  lipid_root, "Tables", "02_selected_consistent_direction_features.csv"
)

input_files <- c(rna_selected_file, rna_expression_file, lipid_selected_file)
if (any(!file.exists(input_files))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = "; "))
}

rna_selected <- fread(rna_selected_file)
lipid_selected <- fread(lipid_selected_file)
rna_input <- as.data.table(readxl::read_excel(
  rna_expression_file, sheet = 1, .name_repair = "minimal"
))

passages <- c("Early", "Middle", "Late")
days <- c("D0", "D1", "D6")
condition_order <- unlist(lapply(passages, function(passage) {
  paste(passage, days, sep = "_")
}))

sample_columns <- names(rna_input)[grepl("^(P11|P22|P33)_D(0|1|6)_[1-3]$", names(rna_input))]
sample_pattern <- "^(P11|P22|P33)_D(0|1|6)_([1-3])$"
sample_parts <- regmatches(sample_columns, regexec(sample_pattern, sample_columns))
passage_lookup <- c(P11 = "Early", P22 = "Middle", P33 = "Late")
canonical_sample_columns <- sprintf(
  "%s_D%s_R%02d",
  passage_lookup[vapply(sample_parts, `[[`, character(1), 2)],
  vapply(sample_parts, `[[`, character(1), 3),
  as.integer(vapply(sample_parts, `[[`, character(1), 4))
)
rna_expression <- as.matrix(rna_input[, ..sample_columns])
storage.mode(rna_expression) <- "double"
rownames(rna_expression) <- rna_input$ensembl_gene_id
colnames(rna_expression) <- canonical_sample_columns
rna_log2 <- log2(rna_expression + 1)

if (!all(rna_selected$ensembl_gene_id %in% rownames(rna_log2))) {
  stop("Some selected RNA genes are missing from the normalized expression matrix")
}

rna_key_matrix <- rna_log2[rna_selected$ensembl_gene_id, , drop = FALSE]
rna_condition_matrix <- sapply(condition_order, function(condition_name) {
  columns <- sprintf("%s_R%02d", condition_name, 1:3)
  if (!all(columns %in% colnames(rna_key_matrix))) {
    stop("Missing RNA samples for condition ", condition_name)
  }
  rowMeans(rna_key_matrix[, columns, drop = FALSE])
})
rownames(rna_condition_matrix) <- rna_selected$ensembl_gene_id

if (!all(condition_order %in% names(lipid_selected))) {
  stop("The selected lipid table does not contain all exact-day condition means")
}
lipid_condition_matrix <- as.matrix(lipid_selected[, ..condition_order])
storage.mode(lipid_condition_matrix) <- "double"
rownames(lipid_condition_matrix) <- lipid_selected$met_id

rna_variable <- apply(rna_condition_matrix, 1, var) > 0
lipid_variable <- apply(lipid_condition_matrix, 1, var) > 0
rna_condition_matrix <- rna_condition_matrix[rna_variable, , drop = FALSE]
lipid_condition_matrix <- lipid_condition_matrix[lipid_variable, , drop = FALSE]

if (!identical(colnames(rna_condition_matrix), colnames(lipid_condition_matrix))) {
  stop("RNA and lipid condition columns are not aligned")
}

matching_design <- data.table(
  condition = condition_order,
  passage = rep(paste(passages, "Passage"), each = length(days)),
  day = rep(days, times = length(passages)),
  integration_unit = "passage-by-day condition mean",
  biological_sample_pairing = FALSE
)

rna_condition_table <- as.data.table(rna_condition_matrix, keep.rownames = "ensembl_gene_id")
rna_condition_table <- merge(
  rna_selected[, .(ensembl_gene_id, symbol, trajectory_direction, average_post_log2FC, average_post_fdr)],
  rna_condition_table,
  by = "ensembl_gene_id", sort = FALSE
)
lipid_condition_table <- as.data.table(lipid_condition_matrix, keep.rownames = "met_id")
lipid_condition_table <- merge(
  lipid_selected[, .(
    met_id, feature, mz, rt_min, accepted_compound_id, accepted_description,
    adducts, formula, trajectory_direction, average_post_log2FC, average_post_fdr
  )],
  lipid_condition_table,
  by = "met_id", sort = FALSE
)

input_summary <- data.table(
  metric = c(
    "integration_scope", "integration_unit", "common_days", "matched_conditions",
    "selected_RNA_genes", "variable_RNA_genes", "selected_lipid_features",
    "variable_lipid_features", "planned_gene_lipid_pairs", "RNA_scale", "lipid_scale"
  ),
  value = c(
    "independently selected key features only",
    "passage-by-day condition means; no individual sample pairing",
    paste(days, collapse = ", "), nrow(matching_design), nrow(rna_selected),
    nrow(rna_condition_matrix), nrow(lipid_selected), nrow(lipid_condition_matrix),
    nrow(rna_condition_matrix) * nrow(lipid_condition_matrix),
    "log2(DESeq2-normalized count + 1)",
    "log2 Progenesis QI normalized abundance after half-minimum imputation"
  )
)

fwrite(matching_design, file.path(table_directory, "01_matching_design.csv"), bom = TRUE)
fwrite(rna_condition_table, file.path(table_directory, "01_RNA_key_gene_condition_means.csv"), bom = TRUE)
fwrite(lipid_condition_table, file.path(table_directory, "01_lipid_key_feature_condition_means.csv"), bom = TRUE)
fwrite(input_summary, file.path(table_directory, "01_integration_input_summary.csv"), bom = TRUE)

saveRDS(
  list(
    matching_design = matching_design,
    rna_matrix = rna_condition_matrix,
    lipid_matrix = lipid_condition_matrix,
    rna_annotation = rna_selected,
    lipid_annotation = lipid_selected
  ),
  file.path(object_directory, "01_key_feature_condition_matrices.rds")
)

cat("Code 1 completed successfully.\n")
cat("Matched conditions:", nrow(matching_design), "\n")
cat("RNA key genes:", nrow(rna_condition_matrix), "\n")
cat("Lipid key features:", nrow(lipid_condition_matrix), "\n")
cat("Planned correlations:", nrow(rna_condition_matrix) * nrow(lipid_condition_matrix), "\n")
