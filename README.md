# Multi-Omics Alzheimer's Disease Figures Pipeline

R scripts for Figures 1-4 in a narrative review on multi-omic and vascular contributions to Alzheimer's disease. Each script starts after cohort-specific cleaning and produces the analysis and figure panels for one part of the manuscript.

## Cohorts

Figures 1-2 use ROSMAP, MSBB, and AMP-AD:
- ROSMAP: proteomics + metabolomics
- MSBB: proteomics + gene expression + ATAC-seq
- AMP-AD: proteomics + metabolomics

Figures 3-4 use DiCAD, Emory_Vascular, MC-CAA, and ROSMAP_Lipidomics_Emory. These cohorts don't overlap with ROSMAP/MSBB/AMP-AD above, no shared subjects.

Data comes from the AD Knowledge Portal / Synapse (M2OVE-AD Consortium and related studies). Raw data isn't included here, see Data access below.

## Scripts

figures1_2_pipeline_circos_heatmap.R
Cross-cohort multi-omic analysis. Loads and cleans each cohort's proteomics plus a second omics layer, maps specimens to individuals and averages replicates, drops features correlated with sex (|r| >= 0.6), and picks top-variance features per layer. Cross-omics correlation uses cis-pair Pearson for MSBB (paired protein/transcript) and a full Spearman grid with FDR for AMP-AD/ROSMAP. Produces the Figure 1 heatmaps and protein-overlap Venn diagram, and the Figure 2 circos plots, FDR summary table, and |r| violin plot.

Needs diablo3_model.RData (a pre-fit ROSMAP DIABLO object) and org.Hs.eg.db for the MSBB gene ID mapping.

figure3_vascular_risk_pipeline.R
Vascular risk factors vs AD-related outcomes in DiCAD, Emory_Vascular, and MC-CAA: forest plot of effects across cohorts, DiCAD amyloid PET by diabetes status, Emory hippocampal volume by hypertension, MC-CAA Aβ42 vs CAA severity, and MC-CAA transcriptomic PCA by Braak stage.

Expects each cohort's cleaned data frames already in the session (dicad_final, emory_check, mccaa_demographics, mccaa_biospecimen, mccaa_biochemical_final, mccaa_pca_data) from earlier cleaning steps.

figure4_HDL_LDL_pipeline.R
Lipid biology across clinical and tissue levels: DiCAD HDL/LDL by diabetes status adjusted for statin use, brain and plasma lipidomics volcano plots by Braak stage, a brain-plasma concordance scatter, and LPE 18:1 boxplots by tissue and Braak stage.

The lipid volcano screens (panels B-C) don't use a multiple-comparison correction, see the manuscript caption for why.

## Requirements

R 4.x with:

```r
install.packages(c("dplyr", "tidyverse", "broom", "ggplot2", "ggrepel",
                    "readxl", "gridExtra", "VennDiagram"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("circlize", "ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"))
```

## Usage

Update the path variables at the top of each script (base_msbb, base_ampad, base_rosmap, data_dir, lipid_dir, out_dir) for your local data. Figures 3-4 need the cleaned cohort data frames already loaded, those come from upstream cleaning scripts not included here. Run a script and it writes PNG panels to out_dir. The three Figure 1 heatmaps render as separate PNGs and get assembled into one figure outside R (Canva/Figma works fine).

## Known data quirks

MC-CAA's biochemical data has to be bridged through the biospecimen table (filtered to the ELISA assay) before joining on specimenID. Joining straight to individualID drops or mismatches rows silently.

AMP-AD protein IDs are a mix of SYMBOL|UniProt and bare UniProt accessions. The bare ones get dropped before cross-cohort comparison since they're unannotated, not missing.

AMP-AD's metabolomics sheet is stored transposed in the source Excel file (rows are metabolites, columns are specimen IDs) and needs to be flipped back before use.

The circos correlation cutoffs (0.12-0.35) were tuned per cohort just so the plots stay readable, they're not significance thresholds. The actual FDR-based significance numbers are in the Figure 2 panel D/E summary.

spear_ampad\$r and spear_rosmap\$r come out as matrices, not flat vectors, and need as.vector() before combining with MSBB's cis-pair vector. Skip that step and you get empty violins with no error.

## Data access

Raw data is not redistributed in this repository. Access requires registration with the AD Knowledge Portal (https://adknowledgeportal.synapse.org/) and, for restricted cohorts, an approved data use agreement.

## Acknowledgement

The data available in the AD Knowledge Portal would not be possible without the participation of research volunteers and the contribution of data by collaborating researchers.

## Citation

Citation for the manuscript to be added on publication :)
