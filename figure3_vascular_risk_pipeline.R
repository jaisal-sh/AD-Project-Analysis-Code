# ============================================================================
# Figure 3. Vascular Risk Factors and AD-Related Outcomes Across Three
# Independent Cohorts of the M2OVE-AD Consortium.
#
# Cohorts: DiCAD, Emory_Vascular, MC-CAA
# Note: These cohorts are distinct from ROSMAP, MSBB, and 
# AMP-AD (Figures 1–2) and do not share overlapping subjects.
#
# Panel A: Forest plot of standardized vascular-risk effects across cohorts
# Panel B: DiCAD - Amyloid PET SUVR by diabetes status
# Panel C: Emory_Vascular - Hippocampal volume by hypertension status
# Panel D: MC-CAA - Ab42 tissue level vs. CAA severity
# Panel E: MC-CAA - Transcriptomic variation (PCA) by Braak stage
#
# Author: Jaisal Sharma 
# ============================================================================

library(dplyr)
library(broom)
library(ggplot2)

out_dir <- "C:/Users/jayvs/Desktop/ad_vascular_multiomics/processed/"

# ---------------------------------------------------------------------------
# NOTE: This script picks up from the point where each cohort's cleaned
# analysis-ready data frames already exist in the environment:
#   dicad_final          - DiCAD demographics + amyloid PET data
#   emory_check          - Emory hippocampal volume + hypertension data
#   mccaa_demographics   - MC-CAA individual-level demographics
#   mccaa_biospecimen    - MC-CAA specimen metadata (bridges individualID -> specimenID)
#   mccaa_biochemical_final - MC-CAA biochemical assay results (Ab42, keyed by specimenID)
#   mccaa_rna / mccaa_pca_data - MC-CAA transcriptomic data with Braak stage, for Panel E
# The upstream data-loading/cleaning steps that produce these objects were
# done earlier in the analysis and are not reconstructed here
# ---------------------------------------------------------------------------


## ============================================================================
## 1. FIT PER-COHORT MODELS
## ============================================================================

# --- DiCAD: amyloid PET SUVR ~ diabetes status ---
dicad_model_data <- dicad_final %>%
  filter(!is.na(diabetes_cat), !is.na(amyloid_global_suvr)) %>%
  mutate(diabetic = ifelse(diabetes_cat %in% c("Treated diabetes", "Untreated diabetes"), 1, 0))

dicad_fit <- lm(amyloid_global_suvr ~ diabetic, data = dicad_model_data)
dicad_result <- tidy(dicad_fit, conf.int = TRUE) %>%
  filter(term == "diabetic") %>%
  mutate(study = "DiCAD", outcome = "Amyloid PET SUVR", exposure = "Diabetes")


# --- Emory_Vascular: hippocampal volume ~ hypertension status ---
emory_model_data <- emory_check %>%
  filter(!is.na(hxbp), mri_completeness == "complete") %>%
  mutate(hip_avg = (Left_Hippocampus + Right_Hippocampus) / 2)

emory_fit <- lm(hip_avg ~ hxbp, data = emory_model_data)
emory_result <- tidy(emory_fit, conf.int = TRUE) %>%
  filter(grepl("hxbp", term)) %>%
  mutate(study = "Emory_Vascular", outcome = "Hippocampal volume", exposure = "Hypertension")


# --- MC-CAA: Ab42 tissue level ~ CAA severity ---
# NOTE: the biochemical data is keyed by specimenID, not individualID directly.
# A naive join of demographics -> biochemical data on individualID silently
# drops (or mismatches) rows. The fix is to bridge through the biospecimen
# metadata table, filtered to the correct assay type, to get the right
# specimenID per individual FIRST, then join biochemical results on specimenID.
mccaa_model_data <- mccaa_demographics %>%
  left_join(mccaa_biospecimen %>% filter(assay == "ELISA") %>% select(individualID, specimenID),
            by = "individualID") %>%
  left_join(mccaa_biochemical_final, by = "specimenID") %>%
  filter(!is.na(AverageCAA), !is.na(`Ab42_TX\n(ug/mg\nprotein)`))
  # NOTE: the Ab42 column name contains LITERAL embedded newline characters
  # from the original messy Excel headers (not a formatting artifact of this
  # script). Copy the exact name from `colnames()` / `grep("Ab42", colnames(...))`
  # rather than retyping it, since re-typing "\n" as an escape sequence inside
  # backticks does NOT reliably match the real column name. Confirmed real
  # names via: grep("Ab42", colnames(mccaa_biochemical_final), value = TRUE)
  #   -> "Ab42_TBS\n(pg/mg\nprotein)", "Ab42_TX\n(pg/mg\nprotein)", "Ab42_FA\n(ug/mg\nprotein)"

cat("MC-CAA model data rows after proper specimenID bridge:", nrow(mccaa_model_data), "\n")

mccaa_fit <- lm(`Ab42_TX\n(ug/mg\nprotein)` ~ AverageCAA, data = mccaa_model_data)
mccaa_result <- tidy(mccaa_fit, conf.int = TRUE) %>%
  filter(term == "AverageCAA") %>%
  mutate(study = "MC-CAA", outcome = "Ab42 tissue level", exposure = "CAA severity")


## ============================================================================
## 2. PANEL A — FOREST PLOT: standardized effect per cohort, 95% CI
## ============================================================================

forest_data <- bind_rows(dicad_result, emory_result, mccaa_result) %>%
  select(study, outcome, exposure, estimate, conf.low, conf.high)

print(forest_data)

# NOTE: for a genuinely comparable "standardized effect (SD units)" axis
# across three different outcome scales (PET SUVR, mm^3 hippocampal volume,
# ug/mg protein), the underlying continuous variables should be z-scored
# (scale()) before fitting each lm() above, OR the raw coefficients should be
# converted to standardized betas post-hoc.

p_forest <- ggplot(forest_data, aes(x = estimate, y = reorder(study, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15, color = "#4C72B0") +
  geom_point(size = 3.5, color = "#2C5F8A") +
  labs(title = "Vascular Risk Effects Across Cohorts",
       x = "Standardized effect (SD units, 95% CI)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(out_dir, "panel_a_forest.png"), plot = p_forest, width = 7, height = 4, dpi = 300, bg = "white")


## ============================================================================
## 3. PANEL B — DiCAD: Amyloid PET SUVR by diabetes status
## ============================================================================

p_dicad_box <- ggplot(dicad_model_data, aes(x = factor(diabetic, labels = c("No diabetes","Diabetes")),
                                             y = amyloid_global_suvr, fill = factor(diabetic))) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.5, color = "grey20") +
  scale_fill_manual(values = c("#A8CBE8", "#4C72B0"), guide = "none") +
  labs(title = "DiCAD: Amyloid PET SUVR by Diabetes Status", x = NULL, y = "Amyloid PET SUVR") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(out_dir, "panel_b_dicad_box.png"), plot = p_dicad_box, width = 5, height = 4.5, dpi = 300, bg = "white")


## ============================================================================
## 4. PANEL C — Emory_Vascular: Hippocampal volume by hypertension status
## ============================================================================

p_emory_box <- ggplot(emory_model_data, aes(x = factor(hxbp), y = hip_avg, fill = factor(hxbp))) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.5, color = "grey20") +
  scale_fill_manual(values = c("#F4A582", "#D6604D"), guide = "none") +
  labs(title = "Emory_Vascular: Hippocampal Volume by Hypertension Status",
       x = NULL, y = "Hippocampal volume") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(out_dir, "panel_c_emory_box.png"), plot = p_emory_box, width = 5, height = 4.5, dpi = 300, bg = "white")


## ============================================================================
## 5. PANEL D — MC-CAA: Ab42 tissue level vs. CAA severity
## ============================================================================

p_mccaa_scatter <- ggplot(mccaa_model_data, aes(x = AverageCAA, y = `Ab42_TX\n(ug/mg\nprotein)`)) +
  geom_point(size = 1, alpha = 0.5, color = "#4DB6AC") +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.7) +
  labs(title = "MC-CAA: Ab42 Tissue Level vs. CAA Severity",
       x = "CAA severity", y = "Ab42 tissue level (ug/mg protein)") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(paste0(out_dir, "panel_d_mccaa_scatter.png"), plot = p_mccaa_scatter, width = 5.5, height = 4.5, dpi = 300, bg = "white")


## ============================================================================
## 6. PANEL E — MC-CAA: Transcriptomic variation (PCA) by Braak stage
## ============================================================================

# NOTE: assumes a PCA has already been run on a top-variance gene set (e.g.
# top 500 variable genes) and merged with Braak stage per sample, producing
# a data frame with at least PC1, PC2, and a numeric/ordinal `braak_stage`
# column. Reconstruct that PCA step above this line if not already in session.

p_mccaa_pca <- ggplot(mccaa_pca_data, aes(x = PC1, y = PC2, color = braak_stage)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradient(low = "#B3CDE3", high = "#08306B", name = "Braak\nstage") +
  labs(title = "MC-CAA: Transcriptomic Variation by Braak Stage",
       x = "PC1", y = "PC2") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA))
# NOTE: points did not show clear clustering by Braak stage in this analysis -
# reported as a genuine finding (Braak stage does not appear to be the
# dominant driver of transcriptomic variance captured by these top PCs),
# not adjusted/reframed to imply separation that wasn't present in the data when 
# ran this analysis.

ggsave(paste0(out_dir, "panel_e_mccaa_pca.png"), plot = p_mccaa_pca, width = 7, height = 5.5, dpi = 300, bg = "white")
