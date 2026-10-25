# Multi-Omics Alzheimer's Disease Figures Pipeline

R scripts for the figures in a narrative review on multi-omic and vascular contributions to Alzheimer's disease. Each script starts after cohort-specific cleaning and produces the analysis and figure panels for one part of the manuscript.

**Figure order is still a work in progress!!** The figure order for the manuscript is being reworked to fit the newer findings and to show the most important results first. The script names below still follow the original draft order (Figures 1-4), so they won't match the final figure numbers. The final mapping of scripts to figures will be updated here once the order is settled.

## Cohorts

Figures 1-2 (as of now) use ROSMAP, MSBB, and AMP-AD:
- ROSMAP: proteomics + metabolomics
- MSBB: proteomics + gene expression + ATAC-seq
- AMP-AD: proteomics + metabolomics

Figures 3-4 (as of now) use DiCAD, Emory_Vascular, MC-CAA, and ROSMAP_Lipidomics_Emory. These cohorts don't overlap with ROSMAP/MSBB/AMP-AD above, no shared subjects.

The Metabolon lipid analysis uses ROSMAP brain Metabolon HD4 metabolomics (brain tissue per the Brain514 assay file name, region not confirmed from the study documentation). It is ROSMAP, so participants probably overlap with ROSMAP_Lipidomics_Emory and with the ROSMAP cohort in Figures 1-2. The overlap hasn't been counted.

Data comes from the AD Knowledge Portal / Synapse (M2OVE-AD Consortium and related studies). Raw data isn't included here, see Data access below.

## Scripts

### figures1_2_pipeline_circos_heatmap.R
Cross-cohort multi-omic analysis. Loads and cleans each cohort's proteomics plus a second omics layer, maps specimens to individuals and averages replicates, drops features correlated with sex (|r| >= 0.6), and picks top-variance features per layer. Cross-omics correlation uses cis-pair Pearson for MSBB (paired protein/transcript) and a full Spearman grid with FDR for AMP-AD/ROSMAP. Produces the Figure 1 heatmaps and protein-overlap Venn diagram, and the Figure 2 circos plots, FDR summary table, and |r| violin plot.

### figure3_vascular_risk_pipeline.R
Vascular risk factors vs AD-related outcomes in DiCAD, Emory_Vascular, and MC-CAA: forest plot of effects across cohorts, DiCAD amyloid PET by diabetes status, Emory hippocampal volume by hypertension, MC-CAA Aβ42 vs CAA severity, and MC-CAA transcriptomic PCA by Braak stage.

Expects each cohort's cleaned data frames already in the session (dicad_final, emory_check, mccaa_demographics, mccaa_biospecimen, mccaa_biochemical_final, mccaa_pca_data) from earlier cleaning steps.

### figure4_HDL_LDL_pipeline.R
Lipid biology across clinical and tissue levels: DiCAD HDL/LDL by diabetes status adjusted for statin use, brain and plasma lipidomics volcano plots by Braak stage, a brain-plasma concordance scatter, and LPE 18:1 boxplots by tissue and Braak stage.

Panel A (HDL/LDL by diabetes) is still used. Panels B-E are superseded by lipidomics_braak_adjusted.R below. The plasma input has 662 rows but only 430 people (technical replicates) and this script treated the rows as independent, so its plasma p-values are too small. It also used no multiple-comparison correction, while the new analysis uses BH FDR.

### lipidomics_braak_adjusted.R
Brain (DLPFC) and plasma lipidomics vs Braak stage (High = Braak >= 4), reading the rosmap_lipid_brain, rosmap_lipid_plasma and rosmap_lipid_clinical sheets of processed/M2OVE.xlsx. Averages plasma technical replicates to one row per person, drops deuterated internal standards, keeps species with at least 80% non-missing values, then fits each species on Braak group + age + sex with BH FDR. Also runs a batch-adjusted sensitivity model and summarises results by lipid class. Produces brain and plasma volcano plots and a class summary panel, plus direction and batch-by-group tables. Standalone, nothing needs to be loaded in the session first.

### metabolon_lipid_braak.R
Tests whether lipid features carry stronger Braak associations than other metabolites in ROSMAP brain Metabolon HD4 data (High = Braak >= 4, adjusted for age and sex, BH FDR). Joins features to the Metabolon data dictionary on CHEM_ID for SUPER_PATHWAY and SUB_PATHWAY. Includes a run-batch sensitivity model, a Braak coding sensitivity check (continuous, >= 3, >= 4, >= 5), and an exploratory, post hoc test of each lipid sub-pathway against the other lipids. The class-level comparison (lipid vs non-lipid) is null. Sub-pathway results are exploratory, and the hexosylceramide direction depends on the Braak cutoff. Standalone.

### cross_platform_lipid_classes.R
Compares median Braak effects by lipid class between the lipidomics results above and ROSMAP brain Metabolon HD4 (9 matched classes). Reads results/rosmap_lipidomics_braak_adjusted.csv and results/inputs/metabolon_braak_effects.csv, a copy of the feature-level table written by metabolon_lipid_braak.R. The class mapping is inexact (see comments in the script), and the Spearman correlation across classes is descriptive only. Run lipidomics_braak_adjusted.R and metabolon_lipid_braak.R first.

## Requirements

R 4.x with:

```r
install.packages(c("dplyr", "tidyverse", "broom", "ggplot2", "ggrepel",
                    "readxl", "gridExtra", "VennDiagram", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("circlize", "ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"))
```

## Usage

Update the path variables at the top of each script (base_msbb, base_ampad, base_rosmap, data_dir, lipid_dir, out_dir) for your local data. Figures 3-4 need the cleaned cohort data frames already loaded, those come from upstream cleaning scripts not included here. Run a script and it writes PNG panels to out_dir. The three Figure 1 heatmaps render as separate PNGs and get assembled into one figure outside R (Canva/Figma works fine).

The three newer scripts (lipidomics_braak_adjusted.R, metabolon_lipid_braak.R, cross_platform_lipid_classes.R) are standalone. Set the paths in the CONFIG block at the top of each, and they write figures to figures/ and tables to results/. 

## Known data quirks

MC-CAA's biochemical data has to be bridged through the biospecimen table (filtered to the ELISA assay) before joining on specimenID. Joining straight to individualID drops or mismatches rows silently.

AMP-AD protein IDs are a mix of SYMBOL|UniProt and bare UniProt accessions. The bare ones get dropped before cross-cohort comparison since they're unannotated, not missing.

AMP-AD's metabolomics sheet is stored transposed in the source Excel file (rows are metabolites, columns are specimen IDs) and needs to be flipped back before use.

The circos correlation cutoffs (0.12-0.35) were tuned per cohort just so the plots stay readable, they're not significance thresholds. The actual FDR-based significance numbers are in the Figure 2 panel D/E summary.

spear_ampad\$r and spear_rosmap\$r come out as matrices, not flat vectors, and need as.vector() before combining with MSBB's cis-pair vector. Skip that step and you get empty violins with no error.

ROSMAP plasma lipidomics has 662 rows but only 430 people. 232 people have two technical replicate rows (same batch, consecutive injections). Average per person before modelling, or p-values come out too small.

ROSMAP brain lipidomics batch is unevenly spread across Braak groups (the High share runs from about 36 to 69% by batch), so batch-adjusted models lose power. Both versions are reported.

Brain lipidomics has deuterated internal standards (for example ChE 18:1-d7) among the species columns. Drop them before any modelling.

Plasma lipidomics is negative ionization only, brain is positive and negative combined, so some classes (LPC for one) can't be compared across tissues.

Metabolon HD4 columns are bare CHEM_IDs. Join to the data dictionary on CHEM_ID to get SUPER_PATHWAY and SUB_PATHWAY. The file also has BATCH_209/305/400/402, SAMPLE_TYPE and MATRIX_DESIGNATION columns that are not metabolites and must be excluded. Each BATCH column belongs to one platform (LIB_ID in the dictionary).

The 80% completeness filter on Metabolon features keeps more lipids (69%) than non-lipids (55%), so class comparisons are not perfectly even.

## Data access

Raw data is not redistributed in this repository. Access requires registration with the AD Knowledge Portal (https://adknowledgeportal.synapse.org/) and, for restricted cohorts, an approved data use agreement with a project proposal and supervision in place.

## Acknowledgements

The data available in the AD Knowledge Portal would not be possible without the participation of research volunteers and the contribution of data by collaborating researchers.

## Citation

Citation for the manuscript to be added on publication :)

- Jaisal