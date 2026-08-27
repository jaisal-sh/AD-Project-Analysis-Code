# ============================================================================
# Multi-Omics Alzheimer's Disease Cohort Analysis -> code for Figures 1-2
# Cohorts: ROSMAP, MSBB, AMP-AD
# 
# Figure 1. Multi-omic feature landscapes and top-variance genes overlap across three diverse AD cohorts.
# Figure 1 - Panels A-C: Heatmaps (One per cohort)
# Figure 1 - Panel D: Upset and Venn Diagram for Top 100 Protein Overlap across Cohorts
# Figure 1 - Panels E-G: Demographics -> not included in this analysis code (no code used to create these)
#
# Figure 2. Cross-Omic Correlation Networks and Multi-Cohort Association Profiles in Alzheimer’s Disease.
# Figure 2 - Panels A-C: Circos Plots (One per cohort)
# Figure 2 - Panel D: FDR Significance Summary Table
# Figure 2 - Panel E: |r| Violin Plot
#
# Pipeline: load & clean each cohort's proteomics + second omics layer,
# match samples to individuals, filter sex-confounded features, select
# top-variance features, test cross-omics correlation (cohort-appropriate
# method), and generate circos / heatmap / summary figures.
#
# Author: Jaisal Sharma
# ============================================================================

library(dplyr)
library(readxl)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(gridExtra)
library(grid)

# ---- EDIT THESE PATHS FOR YOUR ENVIRONMENT --------------------------------
base_msbb   <- "PATH/TO/MSBB multi-omics/"
base_ampad  <- "PATH/TO/AMP-AD multi-omics/"
base_rosmap <- "PATH/TO/ROSMAP dataset/"
out_dir     <- "PATH/TO/outputs/"
# ----------------------------------------------------------------------------


## ============================================================================
## 0. HELPER FUNCTIONS (used across all three cohorts)
## ============================================================================

top_var <- function(mat, n) {
  # Select the n columns (features) with highest variance
  vars <- apply(mat, 2, var, na.rm = TRUE)
  mat[, order(vars, decreasing = TRUE)[1:min(n, ncol(mat))]]
}

impute_median <- function(mat) {
  # Replace non-finite values with the column median
  mat[!is.finite(mat)] <- NA
  col_med <- apply(mat, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) mat[is.na(mat[, j]), j] <- col_med[j]
  mat
}

dedupe_columns <- function(mat) {
  # Average columns sharing the same base name (e.g. "GFAP", "GFAP.1", "GFAP.2")
  base_names <- gsub("\\.\\d+$", "", colnames(mat))
  unique_names <- unique(base_names)
  dedup_mat <- sapply(unique_names, function(nm) {
    cols <- which(base_names == nm)
    if (length(cols) == 1) mat[, cols] else rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(dedup_mat) <- unique_names
  dedup_mat
}

build_individual_matrix <- function(expr_data, id_col = "individualID", exclude_cols) {
  # Collapse multiple specimens per individual (e.g. cell-fraction replicates) to a mean
  numeric_cols <- setdiff(colnames(expr_data), exclude_cols)
  agg <- aggregate(expr_data[, numeric_cols],
                    by = list(individualID = expr_data[[id_col]]),
                    FUN = mean, na.rm = TRUE)
  rownames(agg) <- agg$individualID
  agg$individualID <- NULL
  agg
}

sex_cor <- function(mat, sex_vec) {
  # Per-feature correlation with a binary sex vector, used to flag/drop sex-confounded features
  apply(mat, 2, function(x) suppressWarnings(cor(x, sex_vec, use = "pairwise.complete.obs")))
}

r_to_p <- function(r, n) {
  # Convert a Pearson/Spearman correlation to a two-sided p-value
  t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)
  2 * pt(-abs(t_stat), df = n - 2)
}

r_to_p_spearman_grid <- function(mat1, mat2) {
  # Full N x M Spearman correlation + p-value grid between two feature matrices
  n_row <- ncol(mat1); n_col <- ncol(mat2)
  r_mat <- matrix(NA, n_row, n_col, dimnames = list(colnames(mat1), colnames(mat2)))
  p_mat <- matrix(NA, n_row, n_col, dimnames = list(colnames(mat1), colnames(mat2)))
  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      test <- suppressWarnings(cor.test(mat1[, i], mat2[, j], method = "spearman"))
      r_mat[i, j] <- test$estimate
      p_mat[i, j] <- test$p.value
    }
  }
  list(r = r_mat, p = p_mat)
}

r_for_bonferroni_sig <- function(n, n_tests, alpha = 0.05) {
  # Back-calculate the |r| needed to survive Bonferroni correction, given n and test count
  target_p <- alpha / n_tests
  r_seq <- seq(0.001, 0.99, by = 0.001)
  p_seq <- r_to_p(r_seq, n)
  r_seq[which.min(abs(p_seq - target_p))]
}

clinical_col_names <- c(
  "individualID","specimenID","specimenIdSource","organ","tissue","BrodmannArea",
  "sampleStatus","tissueWeight","tissueVolume","nucleicAcidSource","cellType",
  "fastingState","assay","isPostMortem","samplingAge","visitNumber","exclude",
  "excludeReason","samplingAgeUnits","dataContributionGroup","cohort","species",
  "sex","race","isHispanic","ageDeath","PMI","apoeGenotype","apoe4Status",
  "amyCerad","amyAny","amyThal","amyA","Braak","Braak_numeric","bScore",
  "yearsEducation","causeDeath","mannerDeath","pH","brainWeight","diagnosis",
  "diagnosisCriteria","CDR","CDR_num","plaqueMean","group","mayoDx","reag",
  "ADoutcome","derivedOutcomeBasedOnMayoDx","clinicalMetadataSource",
  "individualID_AMPAD_1.0","Component","treatmentType","treatmentDose","Comments"
)

group_col_3 <- c("Control" = "#4575B4", "MCI" = "#E6820A", "AD" = "#D73027")
group_col_3_alt <- c("NCI" = "#4575B4", "MCI" = "#E6820A", "AD" = "#D73027")   # ROSMAP labels healthy group "NCI"
group_col_2 <- c("Control" = "#4575B4", "AD" = "#D73027")


## ============================================================================
## 1. MSBB — Proteomics + RNA-seq + ATAC-seq (epigenomics)
## ============================================================================

msbb_prot         <- read.csv(paste0(base_msbb, "processed/proteomics/prot_linked_full.csv"), stringsAsFactors = FALSE)
msbb_rna          <- read.csv(paste0(base_msbb, "processed/geneExpression/rna_linked_full.csv"), stringsAsFactors = FALSE)
msbb_atac         <- read.csv(paste0(base_msbb, "processed/epigenetics/atac_linked_full.csv"), stringsAsFactors = FALSE)
msbb_atac_meta    <- read.csv(paste0(base_msbb, "processed/epigenetics/atac_linked_metadata.csv"), stringsAsFactors = FALSE)
msbb_biospecimen  <- read.csv(paste0(base_msbb, "metadata/clinical/MSBB_biospecimen_metadata.csv"), stringsAsFactors = FALSE)
msbb_clinical     <- read.csv(paste0(base_msbb, "metadata/clinical/MSBB_clinical_clean.csv"), stringsAsFactors = FALSE)

# --- map specimenID -> individualID for prot & rna, aggregate to individual level ---
id_map   <- msbb_biospecimen[, c("specimenID", "individualID", "BrodmannArea", "assay")]
msbb_prot$individualID <- id_map$individualID[match(msbb_prot$specimenID, id_map$specimenID)]
msbb_rna$individualID  <- id_map$individualID[match(msbb_rna$specimenID, id_map$specimenID)]

rna_gene_cols <- setdiff(colnames(msbb_rna), c("specimenID", "individualID", "BrodmannArea"))
msbb_rna[rna_gene_cols] <- lapply(msbb_rna[rna_gene_cols], as.numeric)   # RNA cols read in as character - must convert

prot_by_individual <- build_individual_matrix(msbb_prot, exclude_cols = c("specimenID","individualID","BrodmannArea"))
rna_by_individual   <- build_individual_matrix(msbb_rna,  exclude_cols = c("specimenID","individualID","BrodmannArea"))

# --- ATAC: map + aggregate (rowsum-based, faster than aggregate() at 257k columns) ---
atac_id_map <- msbb_atac_meta[, c("specimenID", "individualID")]
msbb_atac$individualID <- atac_id_map$individualID[match(msbb_atac$specimenID, atac_id_map$specimenID)]
peak_cols <- setdiff(colnames(msbb_atac), c("specimenID", "individualID"))
msbb_atac[peak_cols] <- lapply(msbb_atac[peak_cols], as.numeric)
atac_by_individual <- build_individual_matrix(msbb_atac, exclude_cols = c("specimenID","individualID"))

# --- strip any leaked clinical/metadata columns from the omics matrices ---
X_msbb_prot_matched <- prot_by_individual[, !colnames(prot_by_individual) %in% clinical_col_names]
X_msbb_rna_matched  <- rna_by_individual[,  !colnames(rna_by_individual)  %in% clinical_col_names]
atac_by_individual   <- atac_by_individual[, !colnames(atac_by_individual) %in% clinical_col_names]

# --- 3-way sample intersection + diagnosis group (CDR-based: 0=Control, 0.5=MCI, >=3=AD) ---
common_all3 <- Reduce(intersect, list(rownames(X_msbb_prot_matched), rownames(X_msbb_rna_matched), rownames(atac_by_individual)))

msbb_clinical$CDR_num <- as.numeric(as.character(msbb_clinical$CDR))
msbb_dx <- msbb_clinical |>
  mutate(group = case_when(CDR_num == 0 ~ "Control", CDR_num == 0.5 ~ "MCI", CDR_num >= 3 ~ "AD", TRUE ~ NA_character_)) |>
  filter(!is.na(group)) |> select(individualID, group)

Y_df <- data.frame(individualID = common_all3) |> left_join(msbb_dx, by = "individualID")
keep3 <- !is.na(Y_df$group)
common_all3 <- common_all3[keep3]
X_prot_all3 <- as.matrix(X_msbb_prot_matched[common_all3, ])
X_rna_all3  <- as.matrix(X_msbb_rna_matched[common_all3, ])
X_atac_all3 <- as.matrix(atac_by_individual[common_all3, ])
Y_all3 <- factor(Y_df$group[keep3], levels = c("Control","MCI","AD"))

X_prot_clean <- impute_median(X_prot_all3)
X_rna_clean  <- impute_median(X_rna_all3)
X_atac_clean <- impute_median(X_atac_all3)

# --- sex-confound filter (drop features with |r| >= 0.6 vs sex) ---
msbb_sex_df <- data.frame(individualID = common_all3) |> left_join(msbb_clinical |> select(individualID, sex), by = "individualID")
sex_num <- ifelse(msbb_sex_df$sex == "male", 1, ifelse(msbb_sex_df$sex == "female", 0, NA))

X_prot_full_dedup <- dedupe_columns(X_prot_clean)   # collapse duplicate-named protein columns first
prot_sex_r <- sex_cor(X_prot_full_dedup, sex_num)
rna_sex_r  <- sex_cor(X_rna_clean, sex_num)
atac_sex_r <- sex_cor(X_atac_clean, sex_num)

prot_keep <- names(prot_sex_r)[abs(prot_sex_r) < 0.6 | is.na(prot_sex_r)]
rna_keep  <- names(rna_sex_r)[abs(rna_sex_r)   < 0.6 | is.na(rna_sex_r)]
atac_keep <- names(atac_sex_r)[abs(atac_sex_r) < 0.6 | is.na(atac_sex_r)]

# --- top-variance feature selection ---
X_prot_top <- top_var(X_prot_full_dedup[, prot_keep, drop = FALSE], 15)
X_rna_top  <- top_var(X_rna_clean[, rna_keep, drop = FALSE], 20)
X_atac_top <- top_var(X_atac_clean[, atac_keep, drop = FALSE], 10)

cat("MSBB: N =", length(common_all3), "| prot:", ncol(X_prot_top),
    "rna:", ncol(X_rna_top), "atac:", ncol(X_atac_top), "\n")


## ============================================================================
## 2. AMP-AD — Proteomics + Metabolomics
## ============================================================================

ampad_prot         <- read_excel(paste0(base_ampad, "processed/AMPAD_multiomics_processed.xlsx"), sheet = "Proteomics_Frontal")
ampad_met          <- read_excel(paste0(base_ampad, "processed/AMPAD_multiomics_processed.xlsx"), sheet = "Metabolomics")
ampad_link         <- read_excel(paste0(base_ampad, "processed/AMPAD_multiomics_processed.xlsx"), sheet = "Omic_ID_Link")
ampad_indiv_meta   <- read.csv(paste0(base_ampad, "metadata/clinical/AMP-AD_DiverseCohorts_individual_metadata_harmonized.csv"), stringsAsFactors = FALSE)
ampad_met_meta     <- read.csv(paste0(base_ampad, "processed/metabolomics/metab_linked_metadata.csv"), stringsAsFactors = FALSE)

ampad_link_clean <- ampad_link[!is.na(ampad_link$individualID), ]   # link sheet has ~1M padding rows

# --- NOTE: the raw Metabolomics sheet is transposed (rows=metabolite IDs, cols=CMTRX specimenIDs) ---
# Confirmed by checking whether colnames(ampad_met) match specimenID in the assay metadata.
met_mat <- as.matrix(ampad_met[, -1])
rownames(met_mat) <- ampad_met$specimenID
ampad_met_fixed <- as.data.frame(t(met_mat))
ampad_met_fixed$specimenID <- rownames(ampad_met_fixed)
ampad_met_fixed <- ampad_met_fixed[, c("specimenID", setdiff(colnames(ampad_met_fixed), "specimenID"))]
rownames(ampad_met_fixed) <- NULL

# --- diagnosis groups + specimen-to-individual bridges ---
ampad_met_groups <- ampad_met_meta |> filter(ADoutcome %in% c("AD","Control")) |>
  select(specimenID, individualID, ADoutcome, sex) |> rename(group = ADoutcome)

ampad_bridge_prot <- ampad_link_clean |> filter(!is.na(specimenID_frontal)) |>
  left_join(ampad_indiv_meta |> select(individualID, ADoutcome, sex), by = "individualID") |>
  filter(ADoutcome %in% c("AD","Control")) |> rename(group = ADoutcome)

matched_ampad_prot_spec <- intersect(ampad_prot$specimenID, ampad_bridge_prot$specimenID_frontal)

prot_individuals <- ampad_bridge_prot |> filter(specimenID_frontal %in% matched_ampad_prot_spec) |> distinct(individualID)
met_individuals  <- ampad_met_groups  |> filter(specimenID %in% ampad_met_fixed$specimenID) |> distinct(individualID)
common_ampad_individuals <- intersect(prot_individuals$individualID, met_individuals$individualID)

# --- aggregate specimen-level data to individual level (mean across duplicate specimens) ---
prot_spec_to_indiv <- ampad_bridge_prot |>
  filter(specimenID_frontal %in% matched_ampad_prot_spec, individualID %in% common_ampad_individuals) |>
  select(specimenID_frontal, individualID) |> distinct()
ampad_prot_matched <- ampad_prot |> filter(specimenID %in% prot_spec_to_indiv$specimenID_frontal)
ampad_prot_matched$individualID <- prot_spec_to_indiv$individualID[match(ampad_prot_matched$specimenID, prot_spec_to_indiv$specimenID_frontal)]
prot_gene_cols <- setdiff(colnames(ampad_prot_matched), c("specimenID","individualID"))
ampad_prot_matched[prot_gene_cols] <- lapply(ampad_prot_matched[prot_gene_cols], as.numeric)
prot_mat <- as.matrix(ampad_prot_matched[, prot_gene_cols])
X_ampad_prot_indiv <- rowsum(prot_mat, group = ampad_prot_matched$individualID, na.rm = TRUE) /
  rowsum((!is.na(prot_mat)) * 1, group = ampad_prot_matched$individualID)   # NOTE: parens required around !is.na() - operator precedence bug otherwise

met_spec_to_indiv <- ampad_met_groups |>
  filter(specimenID %in% ampad_met_fixed$specimenID, individualID %in% common_ampad_individuals) |>
  select(specimenID, individualID) |> distinct()
ampad_met_matched <- ampad_met_fixed |> filter(specimenID %in% met_spec_to_indiv$specimenID)
ampad_met_matched$individualID <- met_spec_to_indiv$individualID[match(ampad_met_matched$specimenID, met_spec_to_indiv$specimenID)]
met_feature_cols <- setdiff(colnames(ampad_met_matched), c("specimenID","individualID"))
ampad_met_matched[met_feature_cols] <- lapply(ampad_met_matched[met_feature_cols], as.numeric)
met_mat2 <- as.matrix(ampad_met_matched[, met_feature_cols])
X_ampad_met_indiv <- rowsum(met_mat2, group = ampad_met_matched$individualID, na.rm = TRUE) /
  rowsum((!is.na(met_mat2)) * 1, group = ampad_met_matched$individualID)

common_final <- intersect(rownames(X_ampad_prot_indiv), rownames(X_ampad_met_indiv))
X_ampad_prot_final <- X_ampad_prot_indiv[common_final, ]
X_ampad_met_final  <- X_ampad_met_indiv[common_final, ]

Y_ampad_df <- data.frame(individualID = common_final) |>
  left_join(ampad_indiv_meta |> select(individualID, ADoutcome, sex), by = "individualID")
Y_ampad <- factor(Y_ampad_df$ADoutcome, levels = c("Control","AD"))

X_ampad_prot_final <- X_ampad_prot_final[, !colnames(X_ampad_prot_final) %in% clinical_col_names]
X_ampad_met_final  <- X_ampad_met_final[,  !colnames(X_ampad_met_final)  %in% clinical_col_names]

X_ampad_prot_clean <- impute_median(as.matrix(X_ampad_prot_final))
X_ampad_met_clean  <- impute_median(as.matrix(X_ampad_met_final))
X_ampad_prot_dedup <- dedupe_columns(X_ampad_prot_clean)
colnames(X_ampad_prot_dedup) <- sub("^\\|", "", colnames(X_ampad_prot_dedup))   # strip leading "|" from "SYMBOL|UniProt" names

sex_num_ampad <- ifelse(Y_ampad_df$sex == "male", 1, ifelse(Y_ampad_df$sex == "female", 0, NA))
prot_sex_r_a <- sex_cor(X_ampad_prot_dedup, sex_num_ampad)
met_sex_r_a  <- sex_cor(X_ampad_met_clean,  sex_num_ampad)
prot_keep_a <- names(prot_sex_r_a)[abs(prot_sex_r_a) < 0.6 | is.na(prot_sex_r_a)]
met_keep_a  <- names(met_sex_r_a)[abs(met_sex_r_a)   < 0.6 | is.na(met_sex_r_a)]

X_ampad_prot_top <- top_var(X_ampad_prot_dedup[, prot_keep_a, drop = FALSE], 15)
X_ampad_met_top  <- top_var(X_ampad_met_clean[, met_keep_a, drop = FALSE], 20)

cat("AMP-AD: N =", length(common_final), "| prot:", ncol(X_ampad_prot_top), "met:", ncol(X_ampad_met_top), "\n")


## ============================================================================
## 3. ROSMAP — Proteomics + Metabolomics
##    (Unsupervised top-variance selection - NOT the original DIABLO-selected
##     features, which inflated apparent cross-omics correlation via
##     supervised feature selection. See analysis notes.)
## ============================================================================

load(paste0(base_rosmap, "diablo3_model.RData"))   # provides X3 (list: proteomics, metabolomics, variants), Y3 (diagnosis factor)

rosmap_clinical <- read.csv(paste0(base_rosmap, "ROSMAP Proteomics/ROSMAP_clinical_harmonized.csv"), stringsAsFactors = FALSE)

# --- clean gene symbols: X3$proteomics colnames are "SYMBOL.UNIPROTID" ---
colnames(X3$proteomics) <- make.unique(sub("\\..*", "", colnames(X3$proteomics)))

Y3_df <- data.frame(individualID = rownames(X3$proteomics)) |>
  left_join(rosmap_clinical |> select(individualID, sex) |> distinct(individualID, .keep_all = TRUE), by = "individualID")
sex_num_rosmap <- ifelse(Y3_df$sex == "male", 1, ifelse(Y3_df$sex == "female", 0, NA))

prot_sex_r_rosmap <- sex_cor(X3$proteomics, sex_num_rosmap)
met_sex_r_rosmap  <- sex_cor(X3$metabolomics, sex_num_rosmap)
prot_keep_rosmap <- names(prot_sex_r_rosmap)[abs(prot_sex_r_rosmap) < 0.6 | is.na(prot_sex_r_rosmap)]
met_keep_rosmap  <- names(met_sex_r_rosmap)[abs(met_sex_r_rosmap)   < 0.6 | is.na(met_sex_r_rosmap)]

prot_new_top <- top_var(X3$proteomics[, prot_keep_rosmap, drop = FALSE], 15)
met_new_top  <- top_var(X3$metabolomics[, met_keep_rosmap, drop = FALSE], 20)

# --- map metabolite IDs to chemical names via the ROSMAP Metabolon data dictionary ---
met_dict <- read.csv(paste0(base_rosmap, "ROSMAP Metabolomics/Metabolon/QC/ROSMAP Metabolon HD4 Data Dictionary.csv"), stringsAsFactors = FALSE)
met_id_map <- met_dict |> mutate(col_name = paste0("X", CHEM_ID)) |> select(col_name, CHEMICAL_NAME)
col_lookup <- setNames(met_id_map$CHEMICAL_NAME, met_id_map$col_name)
met_names_new <- ifelse(colnames(met_new_top) %in% names(col_lookup) & !is.na(col_lookup[colnames(met_new_top)]),
                         col_lookup[colnames(met_new_top)], colnames(met_new_top))
colnames(met_new_top) <- unname(met_names_new)   # unname() required - circlize breaks on named-vector colnames

# --- for circos/plotting use plain matrices with fully stripped attributes ---
prot_mat_clean <- as.matrix(prot_new_top)
met_mat_clean  <- as.matrix(met_new_top)

cat("ROSMAP: N =", length(Y3), "| prot:", ncol(prot_mat_clean), "met:", ncol(met_mat_clean), "\n")


## ============================================================================
## 4. CROSS-OMICS SIGNIFICANCE TESTING
##
##    Different cohorts get different tests depending on data structure:
##    - MSBB has paired proteomics + RNA-seq -> "cis-pair" test (protein vs.
##      its own gene's transcript level; one test per gene, not a full grid)
##    - AMP-AD / ROSMAP pair proteomics with metabolomics (no shared gene
##      identity) -> full N x M Spearman correlation grid with FDR correction
##
##    Tried and found NOT to help for MSBB/ROSMAP's null results: Bonferroni
##    vs FDR, Pearson vs Spearman, and max-statistic permutation testing (2000
##    perms). These null results reflect genuinely weak cross-omics coupling
##    at this sample size/platform, not a correction-method artifact.
## ============================================================================

# --- AMP-AD & ROSMAP: full grid Spearman + FDR ---
spear_ampad  <- r_to_p_spearman_grid(X_ampad_prot_top, X_ampad_met_top)
spear_rosmap <- r_to_p_spearman_grid(prot_mat_clean, met_mat_clean)

sig_ampad_spear  <- sum(p.adjust(as.vector(spear_ampad$p),  method = "BH") < 0.05)
sig_rosmap_spear <- sum(p.adjust(as.vector(spear_rosmap$p), method = "BH") < 0.05)

# --- MSBB: cis-pair test (requires org.Hs.eg.db for ENSG -> gene symbol mapping) ---
# BiocManager::install("org.Hs.eg.db")   # run once if not installed
library(org.Hs.eg.db)
library(AnnotationDbi)

rna_ensg_ids <- colnames(X_msbb_rna_matched)
ensg_clean <- sub("\\..*", "", rna_ensg_ids)
symbol_map <- AnnotationDbi::select(org.Hs.eg.db, keys = ensg_clean, keytype = "ENSEMBL", columns = "SYMBOL")
symbol_map_unique <- symbol_map[!duplicated(symbol_map$ENSEMBL) & !is.na(symbol_map$SYMBOL), ]
ensg_to_symbol <- setNames(symbol_map_unique$SYMBOL, symbol_map_unique$ENSEMBL)

msbb_protein_symbols <- colnames(X_msbb_prot_matched)
rna_symbols_for_matched_genes <- ensg_to_symbol[ensg_clean]
cis_pairs <- data.frame(ensg = names(rna_symbols_for_matched_genes), rna_col = rna_ensg_ids,
                         symbol = rna_symbols_for_matched_genes, stringsAsFactors = FALSE)
cis_pairs <- cis_pairs[!is.na(cis_pairs$symbol) & cis_pairs$symbol %in% msbb_protein_symbols, ]

common_prot_rna <- intersect(rownames(X_msbb_prot_matched), rownames(X_msbb_rna_matched))
prot_for_cis <- X_msbb_prot_matched[common_prot_rna, ]
rna_for_cis  <- X_msbb_rna_matched[common_prot_rna, ]

cis_r <- cis_p <- numeric(nrow(cis_pairs))
for (i in seq_len(nrow(cis_pairs))) {
  test <- suppressWarnings(cor.test(prot_for_cis[, cis_pairs$symbol[i]], rna_for_cis[, cis_pairs$rna_col[i]],
                                     use = "pairwise.complete.obs"))
  cis_r[i] <- test$estimate; cis_p[i] <- test$p.value
}
cis_pairs$r <- cis_r; cis_pairs$p <- cis_p
cis_pairs <- cis_pairs[!is.na(cis_pairs$p), ]
cis_pairs$p_bonf <- p.adjust(cis_pairs$p, method = "bonferroni")
cis_pairs$p_fdr  <- p.adjust(cis_pairs$p, method = "BH")

cat("\nSignificance summary (FDR < 0.05):\n")
cat("AMP-AD (Spearman grid):", sig_ampad_spear, "/", length(spear_ampad$p), "\n")
cat("MSBB (cis-pair):", sum(cis_pairs$p_fdr < 0.05), "/", nrow(cis_pairs), "\n")
cat("ROSMAP (Spearman grid):", sig_rosmap_spear, "/", length(spear_rosmap$p), "\n")


## ============================================================================
## 5A. FIGURE 1 — HEATMAPS (ComplexHeatmap), one panel per cohort
## ============================================================================

build_cohort_heatmap <- function(feature_mats, layer_labels, y, group_colors,
                                  title, legend = TRUE) {
  # feature_mats: named list of individual-x-feature matrices (already top-var
  #   selected), one per omics layer, all sharing the same row order (samples)
  # layer_labels: character vector, one label per element of feature_mats
  # y: diagnosis factor, same row order as feature_mats
  ht_mats <- lapply(feature_mats, function(m) t(scale(m)))
  combined_mat <- do.call(rbind, ht_mats)
  row_split <- rep(layer_labels, times = sapply(ht_mats, nrow))

  samp_order <- order(y)
  combined_mat <- combined_mat[, samp_order]
  y_ord <- y[samp_order]

  col_anno <- HeatmapAnnotation(Group = y_ord, col = list(Group = group_colors),
                                 show_legend = legend, annotation_name_side = "left")

  Heatmap(combined_mat, name = "z-score",
          col = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")),
          top_annotation = col_anno, row_split = row_split,
          cluster_columns = FALSE, show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 7),
          column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
          row_title_gp = gpar(fontsize = 9))
}

ht_rosmap <- build_cohort_heatmap(list(Proteomics = prot_mat_clean, Metabolomics = met_mat_clean),
                                   c("Proteomics","Metabolomics"), Y3, group_col_3_alt,
                                   "ROSMAP\n(Proteomics + Metabolomics)")
ht_msbb   <- build_cohort_heatmap(list(Proteomics = X_prot_top, GeneExpression = X_rna_top, Epigenomics = X_atac_top),
                                   c("Proteomics","GeneExpression","Epigenomics"), Y_all3, group_col_3,
                                   "MSBB\n(Proteomics + Gene Expression + Epigenomics)", legend = FALSE)
ht_ampad  <- build_cohort_heatmap(list(Proteomics = X_ampad_prot_top, Metabolomics = X_ampad_met_top),
                                   c("Proteomics","Metabolomics"), Y_ampad, group_col_2,
                                   "AMP-AD\n(Proteomics + Metabolomics)", legend = FALSE)

for (nm in c("rosmap","msbb","ampad")) {
  png(paste0(out_dir, "panel_", nm, "_heatmap.png"), width = 1400, height = 2400, res = 200)
  draw(get(paste0("ht_", nm)), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
}
# NOTE: assemble the three PNGs into one figure outside R (PowerPoint/Figma) -
# stitching separately-rendered ComplexHeatmap panels with magick::image_append
# works but does not perfectly align row heights across panels with differing
# feature counts. I found it easiest to use Canva's website to assemble 
# this figure


## ============================================================================
## 5B. FIGURE 2 — CIRCOS PLOTS (circlize), one panel per cohort
## ============================================================================

build_cohort_circos <- function(feature_mats, layer_labels, y, group_colors,
                                 layer_colors, cutoff, out_file, title_line1, title_line2) {
  # feature_mats: named list, e.g. list(Proteomics = X, Metabolomics = Y)
  all_features <- unlist(lapply(feature_mats, colnames), use.names = FALSE)
  feature_layer <- rep(names(feature_mats), times = sapply(feature_mats, ncol))
  layer_idx <- setNames(feature_layer, all_features)
  full_mat <- do.call(cbind, feature_mats)

  # dominant diagnosis group per feature (highest group mean)
  dom_group <- function(mat, y) {
    means <- t(apply(mat, 2, function(x) tapply(x, y, mean, na.rm = TRUE)))
    as.character(unname(apply(means, 1, function(r) names(which.max(r)))))
  }
  feature_dom <- unlist(lapply(feature_mats, dom_group, y = y))
  names(feature_dom) <- all_features

  # cross-layer-only correlation + link table
  cor_all <- cor(full_mat, use = "pairwise.complete.obs")
  pairs_all <- which(abs(cor_all) >= 0, arr.ind = TRUE)
  pairs_all <- pairs_all[pairs_all[,1] < pairs_all[,2], , drop = FALSE]
  cross_layer <- layer_idx[rownames(cor_all)[pairs_all[,1]]] != layer_idx[colnames(cor_all)[pairs_all[,2]]]
  pairs_all <- pairs_all[cross_layer, , drop = FALSE]
  cross_r <- cor_all[pairs_all]

  keep <- abs(cross_r) >= cutoff
  link_df <- data.frame(
    from = as.character(rownames(cor_all)[pairs_all[keep,1]]),
    to   = as.character(colnames(cor_all)[pairs_all[keep,2]]),
    r    = as.numeric(cross_r[keep]),
    stringsAsFactors = FALSE
  )
  link_df$col <- ifelse(link_df$r > 0, "#B2182B60", "#2166AC60")

  circos.clear()
  while (dev.cur() > 1) dev.off()

  png(out_file, width = 3000, height = 3450, res = 300)
  par(oma = c(0, 0, 8, 0)); par(mar = c(1, 1, 1, 1))
  circos.par(gap.degree = 1, start.degree = 90, cell.padding = c(0,0,0,0), track.margin = c(0.033, 0.033))
  circos.initialize(sectors = factor(all_features, levels = all_features), xlim = c(0, 1))

  circos.track(ylim = c(-3, 3), track.height = 0.15, panel.fun = function(x, yv) {
    sec <- get.cell.meta.data("sector.index")
    vals <- as.numeric(scale(full_mat[, sec]))
    dom <- feature_dom[sec]
    col <- if (dom == "AD") "#D73027" else if (dom == "MCI") "#E6820A" else "#4575B4"
    circos.lines(seq(0, 1, length.out = length(vals)), vals, col = col, lwd = 0.8)
  })
  circos.track(ylim = c(0,1), track.height = 0.07, panel.fun = function(x, yv) {
    sec <- get.cell.meta.data("sector.index")
    circos.text(0.5, 0.5, sec, facing = "clockwise", niceFacing = TRUE, cex = 0.45)
  }, bg.border = NA)
  circos.track(ylim = c(0,1), track.height = 0.08, panel.fun = function(x, yv) {
    sec <- get.cell.meta.data("sector.index")
    circos.rect(0, 0, 1, 1, col = group_colors[feature_dom[sec]], border = NA)
  })
  circos.track(ylim = c(0,1), track.height = 0.04, panel.fun = function(x, yv) {
    sec <- get.cell.meta.data("sector.index")
    circos.rect(0, 0, 1, 1, col = layer_colors[layer_idx[sec]], border = NA)
  })
  for (i in seq_len(nrow(link_df))) {
    circos.link(link_df$from[i], 0.5, link_df$to[i], 0.5, col = link_df$col[i], lwd = 1.2)
  }

  mtext(title_line1, side = 3, line = 4, cex = 1.1, font = 2, outer = TRUE)
  mtext(title_line2, side = 3, line = 2, cex = 0.9, outer = TRUE)
  legend("bottomleft", legend = names(group_colors), fill = group_colors, title = "Dominant Group", bty = "n", cex = 0.7)
  legend("bottomright", legend = c("Positive","Negative"), col = c("#B2182B","#2166AC"), lwd = 2, title = "Correlation", bty = "n", cex = 0.7)
  legend("topleft", legend = names(layer_colors), fill = layer_colors, title = "Omic Layer", bty = "n", cex = 0.7)
  dev.off()

  list(link_df = link_df, cross_r = cross_r, all_features = all_features)
}

layer_col_2 <- c("Proteomics" = "grey20", "Metabolomics" = "grey60")
layer_col_3 <- c("Proteomics" = "grey20", "GeneExpression" = "grey60", "Epigenomics" = "grey85")

circos_rosmap <- build_cohort_circos(
  list(Proteomics = prot_mat_clean, Metabolomics = met_mat_clean),
  c("Proteomics","Metabolomics"), Y3, group_col_3_alt, layer_col_2,
  cutoff = 0.15, out_file = paste0(out_dir, "Figure_ROSMAP_circos.png"),
  title_line1 = "Proteomics | Metabolomics (ROSMAP)",
  title_line2 = paste0("n = ", length(Y3), " | 15 proteins, 20 metabolites"))

circos_msbb <- build_cohort_circos(
  list(Proteomics = X_prot_top, GeneExpression = X_rna_top, Epigenomics = X_atac_top),
  c("Proteomics","GeneExpression","Epigenomics"), Y_all3, group_col_3, layer_col_3,
  cutoff = 0.35, out_file = paste0(out_dir, "Figure_MSBB_circos.png"),
  title_line1 = "Proteomics | Gene Expression | Epigenomics (MSBB)",
  title_line2 = "15 proteins, 20 genes, 10 peaks")

circos_ampad <- build_cohort_circos(
  list(Proteomics = X_ampad_prot_top, Metabolomics = X_ampad_met_top),
  c("Proteomics","Metabolomics"), Y_ampad, group_col_2, layer_col_2,
  cutoff = 0.12, out_file = paste0(out_dir, "Figure_AMPAD_circos.png"),
  title_line1 = "Proteomics | Metabolomics (AMP-AD)",
  title_line2 = paste0("n = ", length(Y_ampad), " | 15 proteins, 20 metabolites"))
# NOTE: cutoffs (0.15 / 0.35 / 0.12) were tuned per-cohort to yield a
# comparably legible number of circos links given each cohort's sample size -
# see Section 4 for why these are NOT equivalent to statistical significance
# thresholds (visualized explicitly in the Panel D/E summary below).


## ============================================================================
## 5C. FIGURE 1 PANEL D — TOP-100 PROTEIN OVERLAP ACROSS COHORTS (Venn)
##
##    IMPORTANT: this uses a SEPARATE, larger top-N (100, not 15) selected
##    directly from each cohort's FULL protein panel (not the top-15 subset
##    used in the heatmap/circos panels above). It also requires correcting
##    each cohort's protein ID format to a common gene-symbol convention
##    before comparing - without this fix, overlap comes out at ~0 for all
##    three cohorts, which is a naming-format artifact, not a real biological
##    finding (see notes inline below).
## ============================================================================

# --- ROSMAP: ID format is "SYMBOL.UNIPROTID" -> already fixed in Section 3
#     (colnames(X3$proteomics) cleaned via sub("\\..*", "", ...))
rosmap_all_symbols <- colnames(X3$proteomics)

# --- MSBB: already plain gene symbols, no fix needed ---
msbb_all_proteins <- colnames(X_prot_full_dedup)

# --- AMP-AD: mixed format - "SYMBOL|UNIPROTID" OR bare UniProt accession with
#     no symbol at all. Bare-accession entries cannot be compared by symbol
#     and must be dropped (this is NOT missing data, they are just
#     un-annotated in this dataset). ---
ampad_has_symbol <- grepl("^[A-Za-z0-9]+\\|", colnames(X_ampad_prot_dedup))
X_ampad_prot_symbolonly <- X_ampad_prot_dedup[, ampad_has_symbol]
colnames(X_ampad_prot_symbolonly) <- make.unique(sub("\\|.*", "", colnames(X_ampad_prot_symbolonly)))

n_top_overlap <- 100
rosmap_proteins_big <- colnames(top_var(X3$proteomics, n_top_overlap))
msbb_proteins_big   <- colnames(top_var(X_prot_full_dedup[, prot_keep, drop = FALSE], n_top_overlap))
ampad_proteins_big  <- colnames(top_var(X_ampad_prot_symbolonly[common_final, ], n_top_overlap))

protein_sets_big <- list(ROSMAP = rosmap_proteins_big, MSBB = msbb_proteins_big, `AMP-AD` = ampad_proteins_big)

library(VennDiagram)   # ggvenn/ggVennDiagram had version-compatibility errors in testing; VennDiagram was reliable
venn.plot <- venn.diagram(
  x = protein_sets_big, filename = NULL,
  fill = c("#66A182", "#6A5ACD", "#D4A017"), alpha = 0.4, cex = 1.3, fontface = "bold",
  cat.cex = 1.4, cat.fontface = "bold",
  main = "Overlap of Top 100 Proteins by Variance Across Cohorts", main.cex = 1.6
)
png(paste0(out_dir, "Figure_Venn_top100_proteins.png"), width = 3200, height = 3200, res = 300)
grid::grid.draw(venn.plot)
dev.off()


## ============================================================================
## 5D. FIGURE 2 PANEL D/E — SIGNIFICANCE SUMMARY TABLE + |r| VIOLIN PLOT
##
##    KNOWN BUG (fixed below): spear_ampad$r and spear_rosmap$r are 15x20
##    MATRICES (from r_to_p_spearman_grid()), not flat vectors, while
##    cis_pairs$r IS a flat vector (8480-long). Passing a matrix straight into
##    data.frame(abs_r = abs(matrix)) silently expands it into columns
##    "abs_r.1".."abs_r.20" instead of one "abs_r" column. bind_rows() then
##    only matches on the literal column name "abs_r" - which only the
##    MSBB frame has - so AMP-AD/ROSMAP rows end up entirely NA in the
##    combined violin data, producing two empty violins with no error thrown.
##    FIX: always flatten matrices with as.vector() before building the frame.
## ============================================================================

summary_table <- data.frame(
  Cohort = c("ROSMAP", "MSBB", "AMP-AD"),
  N = c(length(Y3), length(common_prot_rna), n_ampad_final <- length(common_final)),
  Omic_Layers = c("Proteomics + Metabolomics", "Proteomics + RNA-seq (cis)", "Proteomics + Metabolomics"),
  Pairs_Tested = c(length(spear_rosmap$r), nrow(cis_pairs), length(spear_ampad$r)),
  Test_Used = c("Spearman + FDR", "Pearson (cis-pair) + FDR", "Spearman + FDR"),
  FDR_Threshold = rep("q < 0.05", 3),
  N_Significant = c(sig_rosmap_spear, sum(cis_pairs$p_fdr < 0.05), sig_ampad_spear),
  stringsAsFactors = FALSE
)
summary_table$Pct_Significant <- paste0(round(100 * summary_table$N_Significant / summary_table$Pairs_Tested, 2), "%")

tbl_grob <- gridExtra::tableGrob(summary_table, rows = NULL,
                                  theme = gridExtra::ttheme_minimal(
                                    core = list(fg_params = list(fontsize = 9)),
                                    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))))
png(paste0(out_dir, "Figure_PanelD_table.png"), width = 2400, height = 900, res = 300)
grid::grid.draw(tbl_grob)
dev.off()

cohort_colors <- c("ROSMAP" = "#66A182", "MSBB" = "#6A5ACD", "AMP-AD" = "#D4A017")   # matches Venn colors above

r_dist_df <- dplyr::bind_rows(
  data.frame(Cohort = "ROSMAP", abs_r = abs(as.vector(spear_rosmap$r))),   # as.vector() REQUIRED - see note above
  data.frame(Cohort = "MSBB",   abs_r = abs(cis_pairs$r)),                  # already a flat vector
  data.frame(Cohort = "AMP-AD", abs_r = abs(as.vector(spear_ampad$r)))      # as.vector() REQUIRED - see note above
)
r_dist_df$Cohort <- factor(r_dist_df$Cohort, levels = c("ROSMAP", "MSBB", "AMP-AD"))
# sanity check after fixing: table(r_dist_df$Cohort) should read 300 / 8480 / 300

p_violin <- ggplot(r_dist_df, aes(x = Cohort, y = abs_r, fill = Cohort)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = cohort_colors, guide = "none") +
  labs(title = "Distribution of |Correlation| Across All Tested Pairs",
       subtitle = "Pre-FDR-filter; shows full tested-pair distribution per cohort",
       x = NULL, y = "|r| (absolute correlation)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(paste0(out_dir, "Figure_PanelE_violin.png"), plot = p_violin, width = 6, height = 5, dpi = 300, bg = "white")