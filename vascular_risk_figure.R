# vascular_risk_figure.R
# my goal: vascular risk factors vs AD-related outcomes in DiCAD, Emory_Vascular and MC-CAA,
# plus HDL/LDL by diabetes in DiCAD. one row per person everywhere, effects in SD units of the outcome
# per unit of exposure (diabetes yes/no, hypertension yes/no, 1 point of CAA severity)
# MC-CAA: the main estimate keeps every participant; one extreme Ab42 value is shown as a sensitivity
# e4 = any "4" in the APOE genotype string. age and sex come from the sheets named below.

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})
set.seed(1)

# ---- CONFIG ----
proj <- "C:/Users/jayvs/Desktop/ad_vascular_multiomics"
setwd(proj)
wb        <- "processed/M2OVE.xlsx"
dicad_dir <- "downloads/dicad"
mccaa_dir <- "downloads/mc-caa"
out_fig <- "figures"
out_res <- "results"
dir.create(out_fig, showWarnings = FALSE)
dir.create(out_res, showWarnings = FALSE)
cat("working in:", getwd(), "\n")

# ---- helpers ----
z <- function(x) as.numeric(scale(x))
as_num  <- function(x) suppressWarnings(as.numeric(as.character(x)))                  # "." and junk become NA
age_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))   # "90+" becomes 90
e4_flag <- function(g) { g <- as.character(g); ifelse(is.na(g) | g %in% c("", "."), NA_integer_, as.integer(grepl("4", g))) }

# standardised effect: outcome z-scored, exposure in its own units
fit_effect <- function(d, y, x, covars = character(0), study, exposure, outcome, model,
                       subset = "all participants") {
  d <- as.data.frame(d)
  # drop covariates that are mostly missing or constant so one empty column cannot wipe out the sample
  ok <- if (length(covars) == 0) logical(0) else
    sapply(covars, function(v) mean(!is.na(d[[v]])) >= 0.5 && length(unique(na.omit(d[[v]]))) > 1)
  dropped <- covars[!ok]
  covars  <- covars[ok]
  if (length(dropped) > 0) cat("  ", study, "|", model, "|", subset, "- dropped covariates:", paste(dropped, collapse = ", "), "\n")
  dd <- d[complete.cases(d[, c(y, x, covars)]), ]
  dd$yz <- z(dd[[y]])
  fit <- lm(as.formula(paste("yz ~", paste(c(x, covars), collapse = " + "))), data = dd)
  s <- summary(fit)$coefficients[x, ]
  ci <- confint(fit)[x, ]
  tibble(study = study, exposure = exposure, outcome = outcome, model = model, subset = subset,
         covariates = ifelse(length(covars) == 0, "none", paste(covars, collapse = " + ")),
         n = nrow(dd), estimate = unname(s["Estimate"]),
         conf.low = unname(ci[1]), conf.high = unname(ci[2]), p = unname(s["Pr(>|t|)"]))
}
# raw-unit effect (mg/dL for lipids)
fit_raw <- function(d, y, x, covars = character(0), model) {
  d <- as.data.frame(d)
  dd <- d[complete.cases(d[, c(y, x, covars)]), ]
  fit <- lm(as.formula(paste(y, "~", paste(c(x, covars), collapse = " + "))), data = dd)
  s <- summary(fit)$coefficients[x, ]
  ci <- confint(fit)[x, ]
  tibble(outcome = y, model = model, n = nrow(dd), estimate = unname(s["Estimate"]),
         conf.low = unname(ci[1]), conf.high = unname(ci[2]), p = unname(s["Pr(>|t|)"]))
}

# ============================================================================
# DiCAD: diabetes status, amyloid PET, age, sex, APOE
# ============================================================================
pet <- read_csv(file.path(dicad_dir, "DiCAD_amyloidPET_JULY2020.csv"), show_col_types = FALSE) %>%
  transmute(dicad_id = as.character(dicad_id), global_suvr = as_num(global_suvr))
cat("\nDiCAD PET: scans", nrow(pet), "| people", n_distinct(pet$dicad_id), "\n")
pet_p <- pet %>% group_by(dicad_id) %>%
  summarise(suvr = if (all(is.na(global_suvr))) NA_real_ else mean(global_suvr, na.rm = TRUE), .groups = "drop")

dm_raw <- read_csv(file.path(dicad_dir, "DiCAD_metabolomics_JULY2020.csv"), skip = 1,
                   locale = locale(encoding = "Windows-1252"), show_col_types = FALSE,
                   col_select = all_of(c("dicad_id", "diabetes_cat"))) %>%
  filter(!is.na(dicad_id), !is.na(diabetes_cat)) %>% mutate(dicad_id = as.character(dicad_id))
cat("DiCAD diabetes_cat: rows", nrow(dm_raw), "| people", n_distinct(dm_raw$dicad_id), "\n")
print(table(dm_raw$diabetes_cat))
dm <- dm_raw %>% group_by(dicad_id) %>%
  summarise(n_levels = n_distinct(diabetes_cat),
            diabetic = as.integer(any(diabetes_cat %in% c("Treated diabetes", "Untreated diabetes"))),
            .groups = "drop")
cat("people whose diabetes_cat differs across rows (coded as diabetic if ever diabetic):", sum(dm$n_levels > 1), "\n")

ca <- read_csv(file.path(dicad_dir, "DiCAD_clinicalAssessment_JULY2020.csv"), show_col_types = FALSE,
               col_select = all_of(c("dicad_id", "chartage"))) %>%
  transmute(dicad_id = as.character(dicad_id), chartage = as_num(chartage))
age_p <- ca %>% group_by(dicad_id) %>%
  summarise(age = if (all(is.na(chartage))) NA_real_ else min(chartage, na.rm = TRUE), .groups = "drop")

dd <- read_excel(wb, sheet = "dicad_demographics") %>%
  transmute(dicad_id = as.character(dicad_id), sex = factor(sex), e4 = e4_flag(apoeGenotype)) %>%
  distinct(dicad_id, .keep_all = TRUE)

dicad_amy <- dm %>% inner_join(pet_p, by = "dicad_id") %>%
  left_join(age_p, by = "dicad_id") %>% left_join(dd, by = "dicad_id")
stopifnot(!any(duplicated(dicad_amy$dicad_id)))
cat("DiCAD amyloid analysis: people", nrow(dicad_amy), "| with SUVR", sum(!is.na(dicad_amy$suvr)),
    "| missing age", sum(is.na(dicad_amy$age)), "| missing e4", sum(is.na(dicad_amy$e4)), "\n")

# HDL / LDL: one row per person, statin flag without matching nystatin
lab <- read_csv(file.path(dicad_dir, "DiCAD_laboratory_JULY2020.csv"), show_col_types = FALSE,
                col_select = all_of(c("Dicad_ID", "HDL", "LDL"))) %>%
  transmute(dicad_id = as.character(Dicad_ID), HDL = as_num(HDL), LDL = as_num(LDL)) %>%
  filter(!is.na(dicad_id))
cat("DiCAD lipids: rows", nrow(lab), "| people", n_distinct(lab$dicad_id), "\n")
lab_p <- lab %>% group_by(dicad_id) %>%
  summarise(HDL = if (all(is.na(HDL))) NA_real_ else mean(HDL, na.rm = TRUE),
            LDL = if (all(is.na(LDL))) NA_real_ else mean(LDL, na.rm = TRUE), .groups = "drop")

meds <- read_csv(file.path(dicad_dir, "DiCAD_clinicalMedications_JULY2020.csv"), show_col_types = FALSE)
med_cols <- grep("^med_desc_", names(meds), value = TRUE)
is_statin <- function(row) any(grepl("STATIN", row, ignore.case = TRUE) & !grepl("NYSTATIN", row, ignore.case = TRUE), na.rm = TRUE)
meds_f <- meds %>% mutate(dicad_id = as.character(dicad_id),
                          on_statin = apply(across(all_of(med_cols)), 1, is_statin)) %>%
  group_by(dicad_id) %>% summarise(on_statin = any(on_statin), .groups = "drop")

dicad_lipid <- dm %>% inner_join(lab_p, by = "dicad_id") %>% left_join(meds_f, by = "dicad_id") %>%
  left_join(age_p, by = "dicad_id") %>% left_join(dd, by = "dicad_id")
cat("DiCAD lipid analysis: people", nrow(dicad_lipid), "| no medication record (coded as not on statin):",
    sum(is.na(dicad_lipid$on_statin)), "\n")
dicad_lipid <- dicad_lipid %>% mutate(on_statin = ifelse(is.na(on_statin), FALSE, on_statin))
stopifnot(!any(duplicated(dicad_lipid$dicad_id)))

lip_models <- bind_rows(
  fit_raw(dicad_lipid, "HDL", "diabetic", character(0), "Unadjusted"),
  fit_raw(dicad_lipid, "HDL", "diabetic", "on_statin", "Statin-adjusted"),
  fit_raw(dicad_lipid, "HDL", "diabetic", c("on_statin", "age", "sex"), "Statin, age, sex-adjusted"),
  fit_raw(dicad_lipid, "LDL", "diabetic", character(0), "Unadjusted"),
  fit_raw(dicad_lipid, "LDL", "diabetic", "on_statin", "Statin-adjusted"),
  fit_raw(dicad_lipid, "LDL", "diabetic", c("on_statin", "age", "sex"), "Statin, age, sex-adjusted"))
print(lip_models)

# ============================================================================
# Emory_Vascular: hypertension vs hippocampal volume (one baseline visit)
# ============================================================================
em_dem <- read_excel(wb, sheet = "emory_demographics") %>%
  transmute(individualID = as.character(individualID), sex = factor(sex),
            age = age_num(samplingAge), e4 = e4_flag(apoeGenotype)) %>% distinct(individualID, .keep_all = TRUE)
em_cl <- read_excel(wb, sheet = "emory_clinical") %>%
  transmute(individualID = as.character(individualID), htn = as_num(hxbp))   # the "." codes become NA
em_mri <- read_excel(wb, sheet = "emory_mri")
if ("mri_completeness" %in% names(em_mri)) em_mri <- em_mri %>% filter(mri_completeness == "complete")
em_mri <- em_mri %>% transmute(individualID = as.character(individualID),
                               hip = (as_num(Left_Hippocampus) + as_num(Right_Hippocampus)) / 2,
                               icv = as_num(EstimatedTotalIntraCranialVol))
emory <- em_cl %>% inner_join(em_mri, by = "individualID") %>% left_join(em_dem, by = "individualID")
stopifnot(!any(duplicated(emory$individualID)))
cat("\nEmory: people", nrow(emory), "| with HTN and hippocampus", sum(!is.na(emory$htn) & !is.na(emory$hip)),
    "| missing age", sum(is.na(emory$age)), "| missing e4", sum(is.na(emory$e4)), "\n")

# ============================================================================
# MC-CAA: CAA severity vs Ab42 (TX fraction, pg/mg protein)
# ============================================================================
mc_dem <- read_excel(wb, sheet = "mccaa_demographics") %>%
  transmute(individualID = as.character(individualID), caa = as_num(AverageCAA), sex = factor(sex),
            age = age_num(ageDeath), e4 = e4_flag(apoeGenotype)) %>% distinct(individualID, .keep_all = TRUE)
bio <- read_excel(wb, sheet = "mccaa_biochemical")
ab42_col <- grep("^Ab42_TX", names(bio), value = TRUE)
stopifnot(length(ab42_col) == 1)
cat("\nAb42 column used:", gsub("\n", " ", ab42_col), "\n")
bio <- bio %>% transmute(specimenID = as.character(specimenID), ab42 = as_num(.data[[ab42_col]]))
bs <- read_csv(file.path(mccaa_dir, "MC-CAA_biospecimen_metadata.csv"), show_col_types = FALSE)
print(table(bs$assay, useNA = "ifany"))
br <- bs %>% filter(assay == "ELISA") %>%
  transmute(individualID = as.character(individualID), specimenID = as.character(specimenID))
mc <- bio %>% inner_join(br, by = "specimenID")
stopifnot(nrow(mc) > 0)
cat("MC-CAA Ab42: specimens", nrow(mc), "| people", n_distinct(mc$individualID), "\n")
mc_p <- mc %>% group_by(individualID) %>%
  summarise(ab42 = if (all(is.na(ab42))) NA_real_ else mean(ab42, na.rm = TRUE), .groups = "drop") %>%
  inner_join(mc_dem, by = "individualID") %>% filter(!is.na(ab42), !is.na(caa))
stopifnot(!any(duplicated(mc_p$individualID)))
mc_p$ab42_log2 <- log2(mc_p$ab42 + 1)

out_row <- mc_p[which.max(mc_p$ab42), ]
cat("largest Ab42:", out_row$ab42, "| individual:", out_row$individualID, "| CAA:", out_row$caa,
    "(my notes say individual 759, 15365, CAA 0, so check this matches)\n")
mc_x <- mc_p %>% filter(individualID != out_row$individualID)
cat("MC-CAA analysis: people", nrow(mc_p), "| after excluding the largest value", nrow(mc_x), "\n")

# ============================================================================
# forest plot data: unadjusted and adjusted per cohort
# ============================================================================
forest_df <- bind_rows(
  fit_effect(dicad_amy, "suvr", "diabetic", character(0), "DiCAD", "Diabetes", "amyloid PET", "Unadjusted"),
  fit_effect(dicad_amy, "suvr", "diabetic", c("age", "sex", "e4"), "DiCAD", "Diabetes", "amyloid PET", "Adjusted"),
  fit_effect(emory, "hip", "htn", character(0), "Emory_Vascular", "Hypertension", "hippocampal volume", "Unadjusted"),
  fit_effect(emory, "hip", "htn", c("age", "sex", "e4", "icv"), "Emory_Vascular", "Hypertension", "hippocampal volume", "Adjusted"),
  fit_effect(mc_p, "ab42", "caa", character(0), "MC-CAA", "CAA severity", "Ab42", "Unadjusted"),
  fit_effect(mc_p, "ab42", "caa", c("age", "sex", "e4"), "MC-CAA", "CAA severity", "Ab42", "Adjusted"),
  fit_effect(mc_x, "ab42", "caa", character(0), "MC-CAA", "CAA severity", "Ab42", "Unadjusted", "excluding largest value"),
  fit_effect(mc_x, "ab42", "caa", c("age", "sex", "e4"), "MC-CAA", "CAA severity", "Ab42", "Adjusted", "excluding largest value")
) %>% mutate(row = paste0(exposure, " and ", outcome, "\n(", study,
                          ifelse(subset == "all participants", "", paste0(", ", subset)), ")"))
# outlier-robust check: log2 Ab42, all participants
robust <- bind_rows(
  fit_effect(mc_p, "ab42_log2", "caa", character(0), "MC-CAA", "CAA severity", "log2 Ab42", "Unadjusted"),
  fit_effect(mc_p, "ab42_log2", "caa", c("age", "sex", "e4"), "MC-CAA", "CAA severity", "log2 Ab42", "Adjusted"))
print(forest_df %>% select(study, subset, model, n, estimate, conf.low, conf.high, p), n = Inf)
print(robust %>% select(study, outcome, model, n, estimate, conf.low, conf.high, p))

sp <- cor.test(mc_p$caa, mc_p$ab42, method = "spearman", exact = FALSE)
cat("MC-CAA Spearman, all participants: rho =", round(unname(sp$estimate), 3), "p =", signif(sp$p.value, 3), "\n")
write_csv(tibble(test = "Spearman, CAA severity vs Ab42, all participants", n = nrow(mc_p),
                 rho = unname(sp$estimate), p = sp$p.value),
          file.path(out_res, "vascular_mccaa_spearman.csv"))

# ============================================================================
# panels
# ============================================================================
forest_df <- forest_df %>% mutate(row = factor(row, levels = rev(unique(row))),
                                  model = factor(model, levels = c("Unadjusted", "Adjusted")))
pA <- ggplot(forest_df, aes(estimate, row, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.25,
                position = position_dodge(width = 0.6)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("Unadjusted" = "#4C72B0", "Adjusted" = "#DD8452"), name = NULL) +
  labs(x = "standardised effect (SD of outcome per unit of exposure, 95% CI)", y = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

two_box <- function(d, grp, y, labs2, ylab, cols) {
  d <- d %>% mutate(g = factor(.data[[grp]], levels = c(0, 1), labels = labs2)) %>% filter(!is.na(g), !is.na(.data[[y]]))
  nn <- d %>% count(g) %>% mutate(l = paste0(g, "\n(n = ", n, ")"))
  ggplot(d, aes(g, .data[[y]], fill = g)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, size = 0.9, alpha = 0.5, colour = "grey20") +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_discrete(labels = setNames(nn$l, as.character(nn$g))) +
    labs(x = NULL, y = ylab) + theme_minimal(base_size = 11)
}
pB <- two_box(dicad_amy, "diabetic", "suvr", c("NGT or pre-diabetes", "Diabetes"), "amyloid PET SUVR (DiCAD)", c("#A8CBE8", "#4C72B0"))
pC <- two_box(emory, "htn", "hip", c("No hypertension", "Hypertension"), "hippocampal volume (Emory_Vascular)", c("#F4A582", "#D6604D"))

pD <- ggplot(mc_p, aes(caa, ab42)) +
  geom_point(size = 1, alpha = 0.5, colour = "#4DB6AC") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", se = TRUE, linewidth = 0.7) +
  geom_smooth(data = mc_x, method = "lm", formula = y ~ x, colour = "#D95F0E", se = FALSE,
              linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 0.15, y = 7000, hjust = 0, vjust = 1, size = 3,
           label = paste0("largest value off scale\n(", round(out_row$ab42), ", CAA ", out_row$caa, ")")) +
  coord_cartesian(ylim = c(0, 7000)) +
  labs(x = "CAA severity", y = "Ab42, TX fraction (pg/mg protein)",
       caption = "black line: all participants; orange dashed: excluding the largest value") +
  theme_minimal(base_size = 11)

lip_long <- dicad_lipid %>%
  mutate(grp = factor(diabetic, levels = c(0, 1), labels = c("NGT or pre-diabetes", "Diabetes"))) %>%
  select(dicad_id, grp, HDL, LDL) %>% pivot_longer(c(HDL, LDL), names_to = "lipid_type", values_to = "value") %>%
  filter(!is.na(value), !is.na(grp))
ann <- lip_models %>% filter(model == "Statin-adjusted") %>%
  transmute(lipid_type = outcome,
            label = paste0("statin-adjusted difference: ", round(estimate, 1), " mg/dL (p = ", signif(p, 2), ")"))
pE <- ggplot(lip_long, aes(grp, value, fill = lipid_type)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.7, alpha = 0.4, colour = "grey20") +
  facet_wrap(~lipid_type, scales = "free_y") +
  geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("HDL" = "#1baf7a", "LDL" = "#e34948"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = NULL, y = "mg/dL (DiCAD)") + theme_minimal(base_size = 11)

fig_vasc <- pA / (pB | pC | pD) / pE + plot_layout(heights = c(1.1, 1, 1)) + plot_annotation(tag_levels = "A")
print(fig_vasc)

# ---- save ----
ggsave(file.path(out_fig, "vascular_risk_main.png"), fig_vasc, width = 12, height = 14, dpi = 300)
ggsave(file.path(out_fig, "vascular_risk_main.pdf"), fig_vasc, width = 12, height = 14)
write_csv(forest_df %>% select(-row), file.path(out_res, "vascular_risk_effects.csv"))
write_csv(robust,     file.path(out_res, "vascular_mccaa_log2_sensitivity.csv"))
write_csv(lip_models, file.path(out_res, "vascular_dicad_hdl_ldl_models.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_res, "sessionInfo_vascular_risk.txt"))
cat("done, files written to", file.path(getwd(), c(out_fig, out_res)), "\n")