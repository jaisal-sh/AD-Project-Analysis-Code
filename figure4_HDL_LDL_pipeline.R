# ==============================================================================
# Figure 4. Lipid biology across clinical, tissue, and cross-tissue levels
# in vascular-risk and AD-pathology cohorts.
#
# Panel A: HDL/LDL by diabetes status (DiCAD)
# Panel B: Brain lipidomics volcano plot by Braak stage (ROSMAP_Lipidomics_Emory)
# Panel C: Plasma lipidomics volcano plot by Braak stage (ROSMAP_Lipidomics_Emory)
# Panel D: Brain-plasma lipid concordance scatter
# Panel E: LPE 18:1 boxplots, brain vs plasma, by Braak stage
#
# Data sources (AD Knowledge Portal / Synapse):
#   - DiCAD: dicad_laboratory_july2020.csv, dicad_clinicalmedications_july2020.csv,
#            dicad_metabolomics_july2020.csv (diabetes_cat), individual metadata
#   - ROSMAP_Lipidomics_Emory: ROSMAP_Brain_Lipidomic data_batch correction_combined
#            Pos and Neg.csv, ROSMAP_Plasma_lipidomic data_batch correction_ Neg.csv,
#            ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv, ROSMAP_clinical.csv
#
# No correction for multiple comparisons applied to the exploratory lipid screens
# (panels B-C); see figure caption / manuscript note for rationale.
#
# Author: Jaisal Sharma
# ==============================================================================

library(tidyverse)   # dplyr, tidyr, readr, ggplot2
library(broom)        # tidy() model output extraction
library(ggrepel)       # non-overlapping labels, panel D

# ------------------------------------------------------------------------------
# 0. Paths - adjust to your local directory structure
# ------------------------------------------------------------------------------
data_dir  <- "downloads/dicad"
lipid_dir <- "downloads/lipidomics_brain_vs_plasma"
out_dir   <- "processed"

# ==============================================================================
# PANEL A - DiCAD: HDL and LDL by diabetes status, adjusted for statin use
# ==============================================================================

laboratory   <- read_csv(file.path(data_dir, "dicad_laboratory_july2020.csv"))
dicad_meds   <- read_csv(file.path(data_dir, "dicad_clinicalmedications_july2020.csv"))

# flag statin use across all 32 medication slots per subject
med_desc_cols <- grep("^med_desc_", colnames(dicad_meds), value = TRUE)

dicad_meds_flagged <- dicad_meds %>%
  mutate(on_statin = apply(across(all_of(med_desc_cols)), 1, function(row) {
    any(grepl("STATIN", row, ignore.case = TRUE), na.rm = TRUE)
  })) %>%
  select(dicad_id, on_statin) %>%
  distinct(dicad_id, .keep_all = TRUE)

# NOTE: dicad_final / dicad_demographics assumed already built upstream
# (see dicad-cleaning pipeline) - reload session or rebuild diabetes_cat here.
lipid_panel_data <- laboratory %>%
  select(Dicad_ID, HDL, LDL) %>%
  rename(dicad_id = Dicad_ID) %>%
  mutate(dicad_id = as.character(dicad_id)) %>%
  left_join(dicad_final %>% distinct(dicad_id, diabetes_cat) %>%
              mutate(dicad_id = as.character(dicad_id)), by = "dicad_id") %>%
  filter(!is.na(diabetes_cat), !is.na(HDL), !is.na(LDL)) %>%
  mutate(diabetic = ifelse(diabetes_cat %in% c("Treated diabetes", "Untreated diabetes"), 1, 0)) %>%
  left_join(dicad_meds_flagged %>% mutate(dicad_id = as.character(dicad_id)), by = "dicad_id") %>%
  mutate(on_statin = ifelse(is.na(on_statin), FALSE, on_statin))

# unadjusted t-test
t.test(LDL ~ diabetic, data = lipid_panel_data)

# adjusted for statin use
summary(lm(LDL ~ diabetic + on_statin, data = lipid_panel_data))

# figure panel: HDL and LDL side by side, faceted
lipid_panel_long <- lipid_panel_data %>%
  select(dicad_id, diabetic, HDL, LDL) %>%
  pivot_longer(cols = c(HDL, LDL), names_to = "lipid_type", values_to = "value")

p_lipid_panel <- ggplot(lipid_panel_long,
                         aes(x = factor(diabetic, labels = c("No diabetes", "Diabetes")),
                             y = value, fill = lipid_type)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~lipid_type, scales = "free_y") +
  scale_fill_manual(values = c("HDL" = "#1baf7a", "LDL" = "#e34948")) +
  labs(x = NULL, y = "mg/dL", title = "DiCAD: HDL and LDL by diabetes status") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "panel_hdl_ldl_dicad.png"), p_lipid_panel,
       width = 7, height = 5, dpi = 300)

# ==============================================================================
# PANEL B - Brain lipidomics volcano plot (Braak stage)
# ==============================================================================

rosmap_clinical      <- read_csv(file.path(lipid_dir, "ROSMAP_clinical.csv"))
rosmap_lip_biospecimen <- read_csv(file.path(lipid_dir, "ROSMAP_Lipidomics_Emory_biospecimen_metadata.csv"))

rosmap_brain_lipid <- read_csv(
  file.path(lipid_dir, "ROSMAP_Brain_Lipidomic data_batch correction_combined Pos and Neg.csv")
) %>% rename(specimenID = `Biospecimen ID`)

rosmap_brain_lipid_final <- rosmap_brain_lipid %>%
  left_join(rosmap_lip_biospecimen %>% select(individualID, specimenID), by = "specimenID") %>%
  relocate(individualID, specimenID) %>%
  filter(!is.na(individualID))

brain_lipid_cols <- setdiff(colnames(rosmap_brain_lipid_final),
                             c("individualID", "specimenID", "sequence order", "batch"))

brain_lipid_model_data <- rosmap_brain_lipid_final %>%
  mutate(individualID = as.character(individualID)) %>%
  left_join(rosmap_clinical %>% mutate(individualID = as.character(individualID)) %>%
              select(individualID, braaksc), by = "individualID") %>%
  filter(!is.na(braaksc)) %>%
  mutate(braak_group = ifelse(braaksc >= 4, "High Braak", "Low Braak"))

# per-lipid t-test; log2FC = difference of means, since values are already
# log2-scale after batch correction (NOT a log-ratio of raw means)
brain_volcano_results <- lapply(brain_lipid_cols, function(lipid) {
  vals_high <- brain_lipid_model_data[[lipid]][brain_lipid_model_data$braak_group == "High Braak"]
  vals_low  <- brain_lipid_model_data[[lipid]][brain_lipid_model_data$braak_group == "Low Braak"]
  test <- tryCatch(t.test(vals_high, vals_low), error = function(e) NULL)
  if (is.null(test)) return(NULL)
  log2fc <- mean(vals_high, na.rm = TRUE) - mean(vals_low, na.rm = TRUE)
  data.frame(lipid = lipid, log2fc = log2fc, p_value = test$p.value)
}) %>% bind_rows() %>%
  mutate(neg_log10_p = -log10(p_value), significant = p_value < 0.05)

p_brain_volcano <- ggplot(brain_volcano_results, aes(x = log2fc, y = neg_log10_p, color = significant)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "#0C447C", "FALSE" = "grey70")) +
  labs(x = "log2 fold change (High Braak vs Low Braak)", y = "-log10(p-value)",
       title = "Brain lipidomics: Braak stage association") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "panel_brain_volcano.png"), p_brain_volcano,
       width = 6.5, height = 5, dpi = 300)

# ==============================================================================
# PANEL C - Plasma lipidomics volcano plot (Braak stage)
# ==============================================================================

rosmap_plasma_lipid <- read_csv(
  file.path(lipid_dir, "ROSMAP_Plasma_lipidomic data_batch correction_ Neg.csv")
) %>% rename(specimenID = `Biospecimen ID`)

rosmap_plasma_lipid_final <- rosmap_plasma_lipid %>%
  left_join(rosmap_lip_biospecimen %>% select(individualID, specimenID), by = "specimenID") %>%
  relocate(individualID, specimenID) %>%
  filter(!is.na(individualID))

plasma_lipid_cols <- setdiff(colnames(rosmap_plasma_lipid_final),
                              c("individualID", "specimenID", "injection order", "batch"))

plasma_lipid_model_data <- rosmap_plasma_lipid_final %>%
  mutate(individualID = as.character(individualID)) %>%
  left_join(rosmap_clinical %>% mutate(individualID = as.character(individualID)) %>%
              select(individualID, braaksc), by = "individualID") %>%
  filter(!is.na(braaksc)) %>%
  mutate(braak_group = ifelse(braaksc >= 4, "High Braak", "Low Braak"))

plasma_volcano_results <- lapply(plasma_lipid_cols, function(lipid) {
  vals_high <- plasma_lipid_model_data[[lipid]][plasma_lipid_model_data$braak_group == "High Braak"]
  vals_low  <- plasma_lipid_model_data[[lipid]][plasma_lipid_model_data$braak_group == "Low Braak"]
  test <- tryCatch(t.test(vals_high, vals_low), error = function(e) NULL)
  if (is.null(test)) return(NULL)
  log2fc <- mean(vals_high, na.rm = TRUE) - mean(vals_low, na.rm = TRUE)
  data.frame(lipid = lipid, log2fc = log2fc, p_value = test$p.value)
}) %>% bind_rows() %>%
  mutate(neg_log10_p = -log10(p_value), significant = p_value < 0.05)

p_plasma_volcano <- ggplot(plasma_volcano_results, aes(x = log2fc, y = neg_log10_p, color = significant)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "#e34948", "FALSE" = "grey70")) +
  labs(x = "log2 fold change (High Braak vs Low Braak)", y = "-log10(p-value)",
       title = "Plasma lipidomics: Braak stage association") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "panel_plasma_volcano.png"), p_plasma_volcano,
       width = 6.5, height = 5, dpi = 300)

# ==============================================================================
# PANEL D - Brain-plasma lipid concordance
# ==============================================================================

# lipid names differ in formatting between brain and plasma files;
# strip spaces/parens/underscores for approximate name matching
brain_volcano_results$lipid_clean  <- gsub("[ ()]", "", brain_volcano_results$lipid)
plasma_volcano_results$lipid_clean <- gsub("[ ()_]", "", plasma_volcano_results$lipid)

concordance_data <- brain_volcano_results %>%
  select(lipid_clean, brain_log2fc = log2fc, brain_p = p_value) %>%
  inner_join(plasma_volcano_results %>%
               select(lipid_clean, plasma_log2fc = log2fc, plasma_p = p_value),
             by = "lipid_clean")

p_concordance <- ggplot(concordance_data, aes(x = brain_log2fc, y = plasma_log2fc, label = lipid_clean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 3, color = "#0C447C") +
  ggrepel::geom_text_repel(size = 3.5) +
  labs(x = "Brain log2FC (High vs Low Braak)", y = "Plasma log2FC (High vs Low Braak)",
       title = "Brain-plasma lipid concordance") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "panel_concordance.png"), p_concordance,
       width = 6, height = 5.5, dpi = 300)

# ==============================================================================
# PANEL E - LPE 18:1: brain and plasma by Braak stage (strongest cross-tissue hit)
# ==============================================================================

brain_lpe_data  <- brain_lipid_model_data %>%
  select(individualID, braak_group, value = `LPE 18:1`) %>%
  mutate(tissue = "Brain")

plasma_lpe_data <- plasma_lipid_model_data %>%
  select(individualID, braak_group, value = `LPE(18:1)`) %>%
  mutate(tissue = "Plasma")

lpe_combined <- bind_rows(brain_lpe_data, plasma_lpe_data)

p_lpe_panel <- ggplot(lpe_combined, aes(x = braak_group, y = value, fill = braak_group)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~tissue, scales = "free_y") +
  scale_fill_manual(values = c("Low Braak" = "#8FC1E3", "High Braak" = "#0C447C")) +
  labs(x = NULL, y = "LPE 18:1 (log2 intensity)",
       title = "LPE 18:1: brain and plasma by Braak stage") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "panel_lpe_boxplots.png"), p_lpe_panel,
       width = 7, height = 5, dpi = 300)