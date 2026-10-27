# Multi-Omics Alzheimer's Disease Figures Pipeline

R scripts for the figures in a narrative review on multi-omic and vascular contributions to Alzheimer's disease. Each script starts after cohort-specific cleaning and produces the analysis and figure panels for one part of the manuscript.

**Figure order is still a work in progress!!** The figure order for the manuscript is being reworked to fit the newer findings and to show the most important results first. The script names below still follow the original draft order (Figures 1-4), so they won't match the final figure numbers. The final mapping of scripts to figures will be updated here once the order is settled.

## Cohorts

Planned mapping: Fig 1 lipidomics_braak_adjusted.R; Fig 2 metabolon_lipid_braak.R and cross_platform_lipid_classes.R; Fig 3 vascular_risk_figure.R; Fig 4 (context, may go to supplement) figures1_2_pipeline_circos_heatmap.R.

Figures 1-2 (as of now) use ROSMAP, MSBB, and AMP-AD:
- ROSMAP: proteomics + metabolomics
- MSBB: proteomics + gene expression + ATAC-seq
- AMP-AD: proteomics + metabolomics

Figures 3-4 (as of now) use DiCAD, Emory_Vascular, MC-CAA, and ROSMAP_Lipidomics_Emory. DiCAD, Emory_Vascular and MC-CAA don’t overlap with ROSMAP/MSBB/AMP-AD above, no shared subjects. ROSMAP_Lipidomics_Emory is ROSMAP, so it can overlap with the ROSMAP cohort in Figures 1-2 (not counted).

The Metabolon lipid analysis uses ROSMAP brain Metabolon HD4 metabolomics (dorsolateral prefrontal cortex, 514 donors per the AD Knowledge Portal 18.7 data release notes; the cleaned file has 512 people and 509 have Braak stage, age and sex). It is ROSMAP, so participants overlap with ROSMAP_Lipidomics_Emory: 262 of the 384 lipidomics participants (135 High, 127 Low Braak) are among the 509 analysed Metabolon participants. Overlap with the ROSMAP cohort in Figures 1-2 hasn’t been counted. Its metabolites appear to come from the same Metabolon HD4 data.

Data comes from the AD Knowledge Portal / Synapse (M2OVE-AD Consortium and related studies). Raw data isn't included here, see Data access below.

## Scripts

### figures1_2_pipeline_circos_heatmap.R
Cross-cohort multi-omic analysis. Loads and cleans each cohort's proteomics plus a second omics layer, maps specimens to individuals and averages replicates, drops features correlated with sex (|r| >= 0.6), and picks top-variance features per layer. Cross-omics correlation uses cis-pair Pearson for MSBB (paired protein/transcript) and a full Spearman grid with FDR for AMP-AD/ROSMAP. Produces the Figure 1 heatmaps and protein-overlap Venn diagram, and the Figure 2 circos plots, FDR summary table, and |r| violin plot.

### figure3_vascular_risk_pipeline.R
**Superseded by vascular_risk_figure.R below.** Vascular risk factors vs AD-related outcomes in DiCAD, Emory_Vascular, and MC-CAA: forest plot of effects across cohorts, DiCAD amyloid PET by diabetes status, Emory hippocampal volume by hypertension, MC-CAA Aβ42 vs CAA severity, and MC-CAA transcriptomic PCA by Braak stage.

Expects each cohort's cleaned data frames already in the session (dicad_final, emory_check, mccaa_demographics, mccaa_biospecimen, mccaa_biochemical_final, mccaa_pca_data) from earlier cleaning steps.

Kept only for the MC-CAA transcriptomic PCA by Braak stage (panel E), which is not rebuilt in the new script and goes to the supplement. The other panels used DiCAD repeat visits as independent rows (323 rows for 199 people), had no covariate adjustment, and included one extreme MC-CAA Ab42 value without reporting the estimate with and without it.

### figure4_HDL_LDL_pipeline.R
**Superseded.** Lipid biology across clinical and tissue levels: DiCAD HDL/LDL by diabetes status adjusted for statin use, brain and plasma lipidomics volcano plots by Braak stage, a brain-plasma concordance scatter, and LPE 18:1 boxplots by tissue and Braak stage.

Panel A (HDL/LDL by diabetes) is rebuilt per person in vascular_risk_figure.R (the old version used 354 lab rows for 200 people and counted any statin-like string, including nystatin). Panels B-E are superseded by lipidomics_braak_adjusted.R below. The plasma input has 662 rows but only 430 people (technical replicates) and this script treated the rows as independent, so its plasma p-values are too small. It also used no multiple-comparison correction, while the new analysis uses BH FDR.

### vascular_risk_figure.R
Standalone replacement for figure3 and for panel A of figure4. Reads processed/M2OVE.xlsx plus the raw DiCAD and MC-CAA CSVs, collapses everything to one row per person, and fits standardised effects (outcome in SD units per unit of exposure) with and without covariates: diabetes vs amyloid PET (DiCAD), hypertension vs hippocampal volume (Emory_Vascular), and CAA severity vs Ab42 TX fraction (MC-CAA). Covariates are age, sex and APOE e4 where available (plus intracranial volume for Emory). Any covariate that is mostly missing is dropped and the covariates actually used are saved with the estimates. For MC-CAA, all participants is the primary estimate, and the estimate excluding the single largest Ab42 value (individual 759), a log2 version and a Spearman test are saved alongside. Also fits DiCAD HDL and LDL by diabetes status (unadjusted, statin-adjusted, and statin + age + sex adjusted), with the statin flag excluding nystatin. Standalone, nothing needs to be loaded in the session first.

### lipidomics_braak_adjusted.R
Brain (DLPFC) and plasma lipidomics vs Braak stage (High = Braak >= 4), reading the rosmap_lipid_brain, rosmap_lipid_plasma and rosmap_lipid_clinical sheets of processed/M2OVE.xlsx. Averages plasma technical replicates to one row per person, drops deuterated internal standards, keeps species with at least 80% non-missing values, then fits each species on Braak group + age + sex with BH FDR. Also runs a batch-adjusted sensitivity model and summarises results by lipid class. Produces brain and plasma volcano plots and a class summary panel, plus direction and batch-by-group tables. Standalone, nothing needs to be loaded in the session first.

### metabolon_lipid_braak.R
Tests whether lipid features carry stronger Braak associations than other metabolites in ROSMAP brain Metabolon HD4 data (High = Braak >= 4, adjusted for age and sex, BH FDR). Joins features to the Metabolon data dictionary on CHEM_ID for SUPER_PATHWAY and SUB_PATHWAY. Includes a run-batch sensitivity model, a Braak coding sensitivity check (continuous, >= 3, >= 4, >= 5), and an exploratory, post hoc test of each lipid sub-pathway against the other lipids. The class-level comparison (lipid vs non-lipid) is null. Sub-pathway results are exploratory, and the hexosylceramide direction depends on the Braak cutoff. Standalone.

### cross_platform_lipid_classes.R
Compares median Braak effects by lipid class between the lipidomics results above and ROSMAP brain (DLPFC) Metabolon HD4 (9 matched classes). Reads results/rosmap_lipidomics_braak_adjusted.csv and the feature-level table written by metabolon_lipid_braak.R (path set in the CONFIG block), and saves a copy in results/inputs/. Same direction in 5 of 9 classes, Spearman rho 0.43 across classes (p = 0.25, descriptive only). The class mapping is inexact. Metabolon plasmalogens are restricted to GPE species to match lipidomics plasmenylPE (5 species after the 80% filter), Metabolon LPE is the lysophospholipid sub-pathway restricted to GPE, the lysoplasmalogens are a separate sub-pathway and are in no class, and lipidomics plasmanylPE has no Metabolon match. Run lipidomics_braak_adjusted.R and metabolon_lipid_braak.R first.

## Requirements

R R 4.5.1 with:

```r
install.packages(c("dplyr", "tidyverse", "broom", "ggplot2", "ggrepel",
                    "readxl", "gridExtra", "VennDiagram", "patchwork"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("circlize", "ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"))
```

## Usage

Update the path variables at the top of each script (base_msbb, base_ampad, base_rosmap, data_dir, lipid_dir, out_dir) for your local data. Figures 3-4 need the cleaned cohort data frames already loaded, those come from upstream cleaning scripts not included here. Run a script and it writes PNG panels to out_dir. The three Figure 1 heatmaps render as separate PNGs and get assembled into one figure outside R (Canva/Figma works fine).

The four newer scripts (vascular_risk_figure.R, lipidomics_braak_adjusted.R, metabolon_lipid_braak.R, cross_platform_lipid_classes.R) are standalone. Set the paths in the CONFIG block at the top of each, and they write figures to figures/ and tables to results/ holds the summary tables and sessionInfo files from the four newer scripts. They contain no person-level data, and figures are not committed.

## Known data quirks

MC-CAA's biochemical data has to be bridged through the biospecimen table (filtered to the ELISA assay) before joining on specimenID. Joining straight to individualID drops or mismatches rows silently.

MC-CAA's Ab42 column names contain literal newline characters, and the TX fraction is in pg/mg protein (the FA fraction is in ug/mg). Pick the column with grep("^Ab42_TX", names(...)) instead of retyping the name. One participant (individual 759, Ab42_TX about 15365, CAA severity 0) is an extreme value that moves the CAA severity vs Ab42 estimate from p = 0.24 to p = 0.006 when excluded, so both estimates are reported.

DiCAD has repeat records per person. The diabetes_cat table has 323 rows for 199 people, the amyloid PET file 317 scans for 200 people, and the laboratory file 354 rows for 200 people. Collapse to one row per person before modelling, or p-values come out too small. diabetes_cat does not change across a person's rows.

DiCAD diabetes_cat lives in the raw DiCAD_metabolomics CSV (skip = 1, Windows-1252 encoding), not in the M2OVE workbook. Its levels are NGT, Pre-diabetes, Treated diabetes and Untreated diabetes. Only the two treated and untreated categories are coded as diabetic, so the comparison group is NGT plus pre-diabetes.

DiCAD statin flags built by searching medication strings for "STATIN" also match nystatin (an antifungal). Exclude it. 34 of the 199 people have no medication record and are coded as not on a statin.

Emory_Vascular hxbp has stray "." missing codes (2 people). apoeGenotype is empty for every person in the emory_demographics sheet, so APOE can't be used as a covariate, and samplingAge is missing for 39 of 196 people.

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

Brain lipidomics has 92 species columns. Seven are deuterated internal standards (ChE, LPC, PC, PE, PS, SM and TG), leaving 85 species, all above the 80% completeness filter. Plasma has 117 species and no standards.

Age at death recorded as “90+” is set to 90 in the lipidomics and Metabolon scripts. Postmortem interval is missing for only 2 people and is not in the models.

The Metabolon dictionary lists 11 plasmalogens (6 GPE, 5 GPC) but only 9 reach the effects table. One GPC species (P-16:0/20:4) is 53% non-missing and fails the 80% filter, and one GPE species (P-18:0/18:2) is entirely missing. The 80% filter is computed in the 509 analysed people.

MSBB diagnosis groups in Figures 1-2 are CDR-based (0 = Control, 0.5 = MCI, 3 or more = AD), so people with CDR 1 or 2 are excluded. The Figure 1 top-100 protein sets for ROSMAP and AMP-AD are made unique with make.unique, so names ending in .1 or .2 cannot match across cohorts and overlap may be undercounted.

## Data access

Raw data is not redistributed in this repository. Access requires registration with the AD Knowledge Portal (https://adknowledgeportal.synapse.org/) and, for restricted cohorts, an approved data use agreement with a project proposal and supervision in place.

## Acknowledgements

The data available in the AD Knowledge Portal would not be possible without the participation of research volunteers and the contribution of data by collaborating researchers.

## Citation

Citation for the manuscript to be added on publication :)

Jaisal