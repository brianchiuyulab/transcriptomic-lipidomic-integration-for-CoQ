scripts <- c(
  "01_load_key_features.R",
  "02_key_feature_correlations.R",
  "03_lipid_classification.R",
  "04_integration_visualization.R",
  "05_MET767_follow_up.R",
  "06_MET767_KEGG_enrichment.R"
)

script_directory <- getwd()
if (!all(file.exists(scripts))) stop("Open the integration RStudio project")
for (script in scripts) {
  setwd(script_directory)
  source(script, local = .GlobalEnv)
}

cat("Transcriptomic and lipidomic integration completed successfully.\n")
