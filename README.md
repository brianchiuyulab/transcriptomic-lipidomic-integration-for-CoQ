# Transcriptomic and lipidomic key-feature integration

Open `20260905_Final_Analysis/RStudio/OPEN_IN_RSTUDIO.Rproj` and run the numbered scripts in order. `20260905_Final_Analysis/Code availability` contains matching manuscript copies; `20260905_Final_Analysis/Tables` and `20260905_Final_Analysis/Figure` contain final outputs.

This project integrates the 1,498 RNA genes and 60 LC-MS features selected independently in the single-omics analyses.

The integration unit is the passage-by-day condition mean. Only exact shared days (D0, D1 and D6) are used across Early, Middle and Late Passage, giving nine matched conditions. Biological replicates are not paired across omics.

Passage groups are reported as Early Passage, Middle Passage and Late Passage throughout all processed objects, tables and statistical outputs. Dense figure labels use `Early P.`, `Middle P.` and `Late P.`. Axis, legend and annotation text is enlarged and bold for publication readability.

## Lipid annotation provenance

1. Step 1 reads the Progenesis QI Accepted Compound ID and Description. These names remain putative because no authentic-standard confirmation is available.
2. Step 2 searches all 60 selected MET features against LIPID MAPS LMSD at 5 ppm. Neutral mass is used when available; otherwise observed m/z and prespecified positive-mode adducts are searched. Every candidate is retained.
3. Step 3 records the final reporting decision. All single-candidate QI or LMSD annotations are reported as putative with their evidence source; ambiguous, excluded and unassigned results remain MET IDs.

The complete per-feature record is `Tables/03_lipid_annotation_provenance.csv`.

## Analysis steps

1. `01_load_key_features.R` loads the selected features and creates aligned condition-mean matrices.
2. `02_key_feature_correlations.R` calculates all 89,880 Spearman correlations, unadjusted P values and BH-FDR values within the prespecified key-feature search space. High correlation is defined as absolute rho at least 0.90 and BH-FDR below 0.05.
3. `03_lipid_classification.R` imports the systematic 60-feature lipid annotation audit. `03_lipid_annotation_provenance.csv` lists, for every MET feature, the Step 1 QI result, Step 2 LMSD result and Step 3 final annotation decision.
4. `04_integration_visualization.R` generates global correlation heatmaps for all 60 lipid features and for the 26 lipid-class assigned features. Lipid features are grouped by class, and RNA genes are displayed without targeted emphasis.
5. `05_MET767_follow_up.R` performs the prespecified follow-up after prioritizing MET767 as putative CoQ9. It reports the complete 1,498-gene correlation ranking, selects genes at absolute rho at least 0.90 and BH-FDR below 0.05, and generates the correlation significance plot, ranked correlation plot and expression heatmaps. Coq8a is highlighted only at this targeted stage.
6. `06_MET767_KEGG_enrichment.R` performs separate KEGG enrichment analyses for positively and negatively associated genes using the 1,498 selected RNA genes as the background.

Complete correlation tables are stored in `Tables`; compact matrices are stored in `Objects`; publication-oriented figures are stored in `Figure`. Each final figure is supplied as PDF, 600-dpi PNG and SVG. Use SVG for PowerPoint or vector manuscript layout, PNG for raster submission systems, and PDF for archival or editorial workflows. A raw-feature Circos plot is not generated because 3,899 edges cannot be interpreted reliably.
