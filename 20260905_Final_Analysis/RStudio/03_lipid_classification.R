# Lipid annotation and class assignment.

options(stringsAsFactors = FALSE, scipen = 999)

if (.Platform$OS.type == "windows" && !l10n_info()[["UTF-8"]]) {
  Sys.setlocale("LC_CTYPE", "Chinese (Traditional)_Taiwan.utf8")
}

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing R package: data.table")
suppressPackageStartupMessages(library(data.table))

working_directory <- getwd()
if (basename(working_directory) %in% c("C2C12_Integration analysis", "RStudio")) {
  analysis_root <- dirname(working_directory)
} else if (dir.exists(file.path(working_directory, "Tables"))) {
  analysis_root <- working_directory
} else {
  stop("Open the OPEN_IN_RSTUDIO RStudio project")
}

input_file <- file.path(analysis_root, "Tables", "01_lipid_key_feature_condition_means.csv")
data_root <- dirname(dirname(analysis_root))
annotation_file <- file.path(
  data_root, "LC-MS", "20260905_Final_Analysis", "Tables",
  "03_lipid_feature_annotation_audit.csv"
)
if (!file.exists(input_file)) stop("Run Code 1 first")
if (!file.exists(annotation_file)) stop("Run lipidomics Code 3 systematic annotation first")
features <- fread(input_file)
annotation <- fread(annotation_file)
features <- annotation[, .(
  met_id, neutral_mass, qi_concise_name, annotation_evidence, selected_putative_name,
  score, fragmentation_score,
  display_label, broad_class, lmsd_candidate_count,
  unique_candidate_names, eligible_candidate_count, candidate_names,
  eligible_candidate_names, candidate_categories,
  candidate_main_classes, candidate_sub_classes, candidate_broad_classes,
  has_neutral_mass_hit, neutral_mass_unique_names,
  unique_neutral_candidate, unique_mz_candidate,
  best_absolute_delta_ppm, identification_note
)][features, on = "met_id"]

features[, annotation_level := annotation_evidence]
features[, annotation_source := fcase(
  annotation_evidence == "Putative QI annotation", "Progenesis QI report",
  annotation_evidence == "Putative LMSD candidate from QI-inferred neutral mass", "LIPID MAPS LMSD search using QI-inferred neutral mass (5 ppm)",
  annotation_evidence == "Putative LMSD candidate from m/z-adduct search", "LIPID MAPS LMSD m/z-adduct search (5 ppm)",
  annotation_evidence == "Ambiguous LMSD candidates", "LIPID MAPS LMSD ambiguous mass matches (5 ppm)",
  annotation_evidence == "Excluded isotope-labelled LMSD match", "LIPID MAPS LMSD excluded isotope-labelled match",
  default = "No database assignment"
)]
features[, reportable_putative := annotation_evidence %in% c(
  "Putative QI annotation",
  "Putative LMSD candidate from QI-inferred neutral mass",
  "Putative LMSD candidate from m/z-adduct search"
)]
features[, display_label := fifelse(
  reportable_putative,
  paste0(met_id, " (putative ", selected_putative_name, ")"),
  met_id
)]
features[, broad_signal_category := fcase(
  broad_class != "Unassigned", broad_class,
  mz >= 500 & rt_min >= 2, "Unassigned high-mass lipid-like signal",
  mz >= 300, "Unassigned intermediate-mass signal",
  default = "Unassigned low-mass or fragment-like signal"
)]

lipid_annotations <- features[, .(
  met_id,
  feature,
  neutral_mass,
  mz,
  rt_min,
  compound_name = selected_putative_name,
  display_label,
  lipid_class = broad_class,
  annotation_level,
  annotation_source,
  mass_error_ppm = best_absolute_delta_ppm,
  QI_score = score,
  QI_fragmentation_score = fragmentation_score,
  reportable_putative,
  trajectory_direction,
  average_post_log2FC,
  average_post_fdr
)]

annotation_summary <- lipid_annotations[, .N, by = .(annotation_level, lipid_class)]
setorder(annotation_summary, annotation_level, lipid_class)

fwrite(lipid_annotations, file.path(
  analysis_root, "Tables", "03_lipid_annotations.csv"
), bom = TRUE)
fwrite(annotation_summary, file.path(
  analysis_root, "Tables", "03_lipid_annotation_summary.csv"
), bom = TRUE)

cat("Code 3 completed successfully.\n")
cat("Annotated features:", nrow(lipid_annotations), "\n")
