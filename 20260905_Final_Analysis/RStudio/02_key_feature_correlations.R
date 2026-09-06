# Spearman correlations between independently selected RNA and lipid features.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table")
suppressPackageStartupMessages(library(data.table))

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
object_directory <- file.path(analysis_root, "Objects")
input_object <- file.path(object_directory, "01_key_feature_condition_matrices.rds")
if (!file.exists(input_object)) {
  code1 <- file.path("01_load_key_features.R")
  if (!file.exists(code1)) stop("Run Code 1 first")
  source(code1)
}

integration <- readRDS(input_object)
rna_matrix <- integration$rna_matrix
lipid_matrix <- integration$lipid_matrix
rna_annotation <- as.data.table(integration$rna_annotation)
lipid_annotation <- as.data.table(integration$lipid_annotation)

if (!identical(colnames(rna_matrix), colnames(lipid_matrix))) {
  stop("RNA and lipid condition columns are not aligned")
}

n_conditions <- ncol(rna_matrix)
rho <- cor(t(rna_matrix), t(lipid_matrix), method = "spearman")
denominator <- pmax(1 - rho^2, .Machine$double.eps)
statistic <- rho * sqrt((n_conditions - 2) / denominator)
p_value <- 2 * pt(-abs(statistic), df = n_conditions - 2)
fdr <- matrix(
  p.adjust(as.vector(p_value), method = "BH"),
  nrow = nrow(p_value), ncol = ncol(p_value), dimnames = dimnames(p_value)
)

positions <- which(matrix(TRUE, nrow(rho), ncol(rho)), arr.ind = TRUE)
correlation_table <- data.table(
  ensembl_gene_id = rownames(rho)[positions[, 1]],
  met_id = colnames(rho)[positions[, 2]],
  n_conditions = n_conditions,
  spearman_rho = rho[positions],
  p_value_unadjusted = p_value[positions],
  fdr_BH_key_feature_space = fdr[positions]
)

gene_fields <- rna_annotation[, .(
  ensembl_gene_id,
  symbol,
  RNA_direction = trajectory_direction,
  RNA_average_post_log2FC = average_post_log2FC,
  RNA_average_post_fdr = average_post_fdr
)]
lipid_fields <- lipid_annotation[, .(
  met_id,
  feature,
  mz,
  rt_min,
  accepted_compound_id,
  accepted_description,
  adducts,
  formula,
  lipid_direction = trajectory_direction,
  lipid_average_post_log2FC = average_post_log2FC,
  lipid_average_post_fdr = average_post_fdr
)]
correlation_table <- gene_fields[correlation_table, on = "ensembl_gene_id"]
correlation_table <- lipid_fields[correlation_table, on = "met_id"]
correlation_table[, absolute_spearman_rho := abs(spearman_rho)]
correlation_table[, direction_relationship := fifelse(
  RNA_direction == lipid_direction, "same passage direction", "opposite passage direction"
)]
setorder(correlation_table, -absolute_spearman_rho, fdr_BH_key_feature_space, met_id, symbol)

high_correlation <- correlation_table[
  absolute_spearman_rho >= 0.90 & fdr_BH_key_feature_space < 0.05
]

leave_one_out_summary <- function(gene_id, feature_id) {
  x <- as.numeric(rna_matrix[gene_id, ])
  y <- as.numeric(lipid_matrix[feature_id, ])
  estimates <- vapply(seq_along(x), function(i) {
    cor(x[-i], y[-i], method = "spearman")
  }, numeric(1))
  data.table(
    leave_one_out_min_rho = min(estimates),
    leave_one_out_max_rho = max(estimates),
    leave_one_out_min_absolute_rho = min(abs(estimates)),
    leave_one_out_sign_consistent = all(sign(estimates) == sign(cor(x, y, method = "spearman")))
  )
}

if (nrow(high_correlation)) {
  stability <- rbindlist(lapply(seq_len(nrow(high_correlation)), function(i) {
    leave_one_out_summary(
      high_correlation$ensembl_gene_id[i],
      high_correlation$met_id[i]
    )
  }))
  high_correlation <- cbind(high_correlation, stability)
}

feature_connectivity <- correlation_table[, .(
  tested_genes = .N,
  genes_abs_rho_ge_0_90 = sum(absolute_spearman_rho >= 0.90),
  positive_genes_abs_rho_ge_0_90 = sum(spearman_rho >= 0.90),
  negative_genes_abs_rho_ge_0_90 = sum(spearman_rho <= -0.90),
  maximum_absolute_rho = max(absolute_spearman_rho)
), by = .(met_id, feature, mz, rt_min, accepted_compound_id, accepted_description, lipid_direction)]
setorder(feature_connectivity, -genes_abs_rho_ge_0_90, -maximum_absolute_rho)

analysis_summary <- data.table(
  metric = c(
    "RNA_key_genes", "lipid_key_features", "matched_conditions", "tested_pairs",
    "correlation_unit", "individual_sample_pairing", "p_value_method", "FDR_scope",
    "rho_cutoff", "FDR_cutoff", "high_correlation_pairs", "genes_in_high_correlation_pairs",
    "lipid_features_in_high_correlation_pairs", "positive_high_correlation_pairs",
    "negative_high_correlation_pairs", "stable_sign_high_correlation_pairs"
  ),
  value = c(
    nrow(rna_matrix), nrow(lipid_matrix), n_conditions, nrow(correlation_table),
    "9 passage-by-day condition means", "FALSE",
    "two-sided t approximation of Spearman rho; df = 7",
    "BH adjustment across all tested key gene-lipid pairs",
    0.90, 0.05, nrow(high_correlation), uniqueN(high_correlation$ensembl_gene_id),
    uniqueN(high_correlation$met_id), sum(high_correlation$spearman_rho > 0),
    sum(high_correlation$spearman_rho < 0),
    sum(high_correlation$leave_one_out_sign_consistent, na.rm = TRUE)
  )
)

fwrite(correlation_table, file.path(table_directory, "02_all_key_feature_correlations.csv"), bom = TRUE)
fwrite(high_correlation, file.path(table_directory, "02_high_correlation_pairs.csv"), bom = TRUE)
fwrite(feature_connectivity, file.path(table_directory, "02_lipid_feature_correlation_connectivity.csv"), bom = TRUE)
fwrite(analysis_summary, file.path(table_directory, "02_correlation_summary.csv"), bom = TRUE)

saveRDS(
  list(rho = rho, p_value_unadjusted = p_value, fdr_BH_key_feature_space = fdr),
  file.path(object_directory, "02_key_feature_correlation_matrices.rds")
)

cat("Code 2 completed successfully.\n")
cat("Tested pairs:", nrow(correlation_table), "\n")
cat("High-correlation pairs:", nrow(high_correlation), "\n")
cat("Genes represented:", uniqueN(high_correlation$ensembl_gene_id), "\n")
cat("Lipid features represented:", uniqueN(high_correlation$met_id), "\n")
