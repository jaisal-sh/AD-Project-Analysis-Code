# metabolon_lipid_braak.R
# my goal: do lipid features carry stronger Braak associations than other metabolites in ROSMAP?
# result: class level is null, sub-pathways move in opposite directions (exploratory, post hoc)
# tissue: brain (file name ROSMAP_Metabolon_HD4_Brain514, all 512 clean IDs are in the raw file)
# region: DLPFC per the AD Knowledge Portal 18.7 release notes (514 donors)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
set.seed(1)

# ---- CONFIG ----
setwd("C:/Users/jayvs/Desktop/ROSMAP Multi-Omics")
metab_file <- "processed/metabolomics/ROSMAP_metabolon_HD4_clean.csv"
meta_file  <- "C:/Users/jayvs/Desktop/ad_vascular_multiomics/downloads/lipidomics_brain_vs_plasma/ROSMAP_clinical.csv"
annot_file <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Metabolomics/Metabolon/QC/ROSMAP Metabolon HD4 Data Dictionary.csv"
id_col      <- "individualID"
braak_col   <- "braaksc"
age_col     <- "age_death"
sex_col     <- "msex"
lipid_label <- "Lipid"
non_metab   <- c("BATCH_209", "BATCH_305", "BATCH_400", "BATCH_402",
                 "SAMPLE_TYPE", "MATRIX_DESIGNATION")
out_fig <- "figures"
out_res <- "results"
dir.create(out_fig, showWarnings = FALSE)
dir.create(out_res, showWarnings = FALSE)
cat("working in:", getwd(), "\n")

# ---- load (one time only) ----
metab <- read_csv(metab_file, show_col_types = FALSE)
meta  <- read_csv(meta_file,  show_col_types = FALSE)
stopifnot(!any(duplicated(metab[[id_col]])))   # repeat samples would double count

dict_full <- read_csv(annot_file, show_col_types = FALSE)
annot <- dict_full %>%
  transmute(feature = as.character(CHEM_ID), super_pathway = SUPER_PATHWAY,
            sub_pathway = SUB_PATHWAY) %>%
  distinct(feature, .keep_all = TRUE)
lib_map <- dict_full %>%   # each platform has its own batch column (LIB_ID 209/305/400/402)
  transmute(feature = as.character(CHEM_ID), batch_col = paste0("BATCH_", LIB_ID))

# ---- join clinical, Braak groups ----
dat <- metab %>%
  inner_join(meta %>% select(all_of(c(id_col, braak_col, age_col, sex_col))), by = id_col) %>%
  mutate(
    braak = as.numeric(.data[[braak_col]]),
    age   = as.numeric(gsub("\\+", "", as.character(.data[[age_col]]))),  # handles "90+"
    sex   = factor(.data[[sex_col]])
  ) %>%
  filter(!is.na(braak), !is.na(age), !is.na(sex)) %>%
  mutate(braak_group = factor(ifelse(braak >= 4, "High", "Low"), levels = c("Low", "High")))

all_feat <- setdiff(names(metab), c(id_col, non_metab))
features <- all_feat[sapply(all_feat, function(f) is.numeric(dat[[f]]) && mean(!is.na(dat[[f]])) >= 0.8)]

# check 1: counts before i trust anything downstream
cat("subjects:", nrow(dat), "| High:", sum(dat$braak_group == "High"),
    "| Low:", sum(dat$braak_group == "Low"), "\n")
cat("features before filter:", length(all_feat), "| kept:", length(features),
    "| kept with annotation:", sum(features %in% annot$feature), "\n")
if (sum(features %in% annot$feature) == 0) stop("annotation join failed, nothing matched")
stopifnot(all(lib_map$batch_col[match(features, lib_map$feature)] %in% names(dat)))

# check 2: did the 80% filter hit lipids and non-lipids evenly
filter_tab <- tibble(feature = all_feat, kept = all_feat %in% features) %>%
  left_join(annot, by = "feature") %>%
  mutate(class = case_when(is.na(super_pathway) ~ "No annotation",
                           super_pathway == lipid_label ~ "Lipid", TRUE ~ "Non-lipid")) %>%
  group_by(class) %>%
  summarise(total = n(), kept = sum(kept), pct_kept = round(100 * kept / total))
print(filter_tab)

# ---- main model: z-scored metabolite ~ Braak group + age + sex ----
res <- map_dfr(features, function(f) {
  d <- dat %>% mutate(y = as.numeric(scale(.data[[f]])))
  s <- summary(lm(y ~ braak_group + age + sex, data = d))$coefficients["braak_groupHigh", ]
  tibble(feature = f, beta = s["Estimate"], p = s["Pr(>|t|)"])
}) %>%
  mutate(fdr = p.adjust(p, method = "BH"), abs_beta = abs(beta)) %>%
  left_join(annot, by = "feature") %>%
  filter(!is.na(super_pathway)) %>%
  mutate(class = ifelse(super_pathway == lipid_label, "Lipid", "Non-lipid"))

print(res %>% group_by(class) %>%
        summarise(n = n(), n_fdr05 = sum(fdr < 0.05), median_abs_beta = median(abs_beta)))

w  <- wilcox.test(abs_beta ~ class, data = res)
ft <- fisher.test(table(res$class, factor(res$fdr < 0.05, levels = c(TRUE, FALSE))))
cat("Wilcoxon p:", signif(w$p.value, 3), "| Fisher p:", signif(ft$p.value, 3), "\n")

# ---- sensitivity: add run batch for each feature's platform ----
res_b <- map_dfr(features, function(f) {
  b <- lib_map$batch_col[match(f, lib_map$feature)]
  d <- dat %>% mutate(y = as.numeric(scale(.data[[f]])), batch = factor(.data[[b]]))
  s <- summary(lm(y ~ braak_group + age + sex + batch, data = d))$coefficients["braak_groupHigh", ]
  tibble(feature = f, beta_b = s["Estimate"], p_b = s["Pr(>|t|)"])
}) %>% mutate(fdr_b = p.adjust(p_b, "BH"))

cmp <- res %>% inner_join(res_b, by = "feature")
cat("beta correlation, unadjusted vs batch-adjusted:", round(cor(cmp$beta, cmp$beta_b), 3), "\n")
print(cmp %>% group_by(class) %>% summarise(sig_before = sum(fdr < 0.05), sig_after = sum(fdr_b < 0.05)))
wb <- wilcox.test(abs(beta_b) ~ class, data = cmp)
cat("Wilcoxon p with batch:", signif(wb$p.value, 3), "\n")

# ---- lipid sub-pathways (exploratory, post hoc) ----
lip <- res %>% filter(class == "Lipid")
sub_tab <- lip %>% group_by(sub_pathway) %>%
  summarise(n = n(), n_sig = sum(fdr < 0.05), pct_down = round(100 * mean(beta < 0)),
            med_beta = round(median(beta), 2), .groups = "drop") %>%
  arrange(med_beta)

# each sub-pathway (n >= 5) vs the other lipids, corrected across sub-pathways
sub_test <- lip %>% group_by(sub_pathway) %>% filter(n() >= 5) %>%
  group_modify(~ {
    inn <- .x$beta
    out <- lip$beta[lip$sub_pathway != .y$sub_pathway]
    tibble(n = length(inn), med_beta = round(median(inn), 2), p = wilcox.test(inn, out)$p.value)
  }) %>%
  ungroup() %>% mutate(fdr_sub = p.adjust(p, "BH")) %>% arrange(p)
print(sub_test, n = Inf)

# ---- Braak coding sensitivity: signs and the class test are comparable, magnitudes are not ----
fit_one <- function(f, coding) {
  d <- dat %>% mutate(y = as.numeric(scale(.data[[f]])),
                      g = switch(coding, cont = braak, ge3 = as.numeric(braak >= 3),
                                 ge4 = as.numeric(braak >= 4), ge5 = as.numeric(braak >= 5)))
  summary(lm(y ~ g + age + sex, data = d))$coefficients["g", "Estimate"]
}
focus <- c("Phosphatidylethanolamine (PE)", "Phosphatidylcholine (PC)", "Hexosylceramides (HCER)",
           "Lysophospholipid", "Sphingomyelins")  # three from the sub-pathway test, two from my lipidomics classes, so post hoc

sens <- map_dfr(c("cont", "ge3", "ge4", "ge5"), function(cd) {
  tibble(feature = features, coding = cd, beta = map_dbl(features, fit_one, coding = cd))
}) %>% inner_join(res %>% select(feature, class, sub_pathway, beta_main = beta), by = "feature")

sens_class <- sens %>% group_by(coding) %>%
  summarise(r_vs_main = round(cor(beta, beta_main), 3),
            wilcox_p = signif(wilcox.test(abs(beta[class == "Lipid"]), abs(beta[class == "Non-lipid"]))$p.value, 3),
            .groups = "drop")
print(sens_class)
sens_focus <- sens %>% filter(sub_pathway %in% focus) %>%
  group_by(coding, sub_pathway) %>% summarise(med_beta = round(median(beta), 2), .groups = "drop") %>%
  pivot_wider(names_from = coding, values_from = med_beta)
print(sens_focus)

# ---- panels ----
cols <- c("Lipid" = "#2C7FB8", "Non-lipid" = "grey60")
cnt <- res %>% group_by(class) %>%
  summarise(lab = paste0(first(class), "\n", sum(fdr < 0.05), "/", n(), " FDR < 0.05"))

pA <- ggplot(res, aes(class, abs_beta, fill = class)) +
  geom_violin(alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(aes(colour = fdr < 0.05), width = 0.08, size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey75"), name = "FDR < 0.05") +
  scale_x_discrete(labels = setNames(cnt$lab, cnt$class)) +
  annotate("text", x = 1.5, y = max(res$abs_beta) * 1.05,
           label = paste0("Wilcoxon p = ", signif(w$p.value, 2)), size = 3.5) +
  labs(x = NULL, y = "|standardised effect| (High vs Low Braak)") +
  theme_minimal(base_size = 12)

by_path <- res %>% group_by(super_pathway) %>%
  summarise(n = n(), n_sig = sum(fdr < 0.05), med = median(abs_beta), .groups = "drop") %>%
  mutate(label = paste0(n_sig, "/", n), super_pathway = fct_reorder(super_pathway, med))
pB <- res %>%
  mutate(super_pathway = factor(super_pathway, levels = levels(by_path$super_pathway))) %>%
  ggplot(aes(super_pathway, abs_beta)) +
  geom_boxplot(aes(fill = super_pathway == lipid_label), outlier.size = 0.6) +
  geom_text(data = by_path, aes(y = max(res$abs_beta) * 1.05, label = label), size = 3) +
  scale_fill_manual(values = c("TRUE" = "#2C7FB8", "FALSE" = "grey80"), guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "|standardised effect|", caption = "labels = FDR-significant / features tested") +
  theme_minimal(base_size = 12)

pC <- sub_test %>%
  mutate(short = sub_pathway %>%
           str_replace("Fatty Acid Metabolism \\(Acyl Carnitine, (.*)\\)", "Acyl carnitine, \\1") %>%
           str_replace(" \\(n3 and n6\\)", "") %>%
           str_replace("Long Chain (.*) Fatty Acid", "LC \\1 FA") %>%
           str_replace("Long Chain", "LC"),
         lab = paste0(short, " (n=", n, ")"),
         sig = fdr_sub < 0.05,
         lab = fct_reorder(lab, med_beta),
         dir = ifelse(med_beta < 0, "Lower in high Braak", "Higher in high Braak")) %>%
  ggplot(aes(med_beta, lab, fill = dir, alpha = sig)) +
  geom_col() +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35), guide = "none") +
  scale_fill_manual(values = c("Lower in high Braak" = "#2C7FB8",
                               "Higher in high Braak" = "#D95F0E"), name = NULL) +
  labs(x = "median standardised effect (High vs Low Braak)", y = NULL,
       caption = "bars = median effect within group; solid = differs from other lipids\n(FDR < 0.05 across 21 sub-pathways, exploratory); groups with n >= 5 only\nLC = long chain, FA = fatty acid; HCER direction depends on Braak cutoff (see supplement)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# feature-level effects in the focus groups (post hoc)
sC_dat <- lip %>% filter(sub_pathway %in% focus) %>%
  group_by(sub_pathway) %>% mutate(lab = paste0(sub_pathway, " (n=", n(), ")")) %>% ungroup()
lab_levels <- sC_dat %>% distinct(sub_pathway, lab) %>% arrange(match(sub_pathway, rev(focus))) %>% pull(lab)
sC_dat$lab <- factor(sC_dat$lab, levels = lab_levels)
pD <- ggplot(sC_dat, aes(beta, lab)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  geom_jitter(aes(colour = fdr < 0.05), height = 0.12, size = 1.3, alpha = 0.7) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey70"), name = "FDR < 0.05") +
  labs(x = "standardised effect per feature (High vs Low Braak)", y = NULL,
       caption = "groups selected post hoc (see text)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# supplement: batch sensitivity and volcano
pE <- ggplot(cmp, aes(beta, beta_b, colour = class)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "effect (age + sex)", y = "effect (+ run batch)",
       caption = paste0("r = ", round(cor(cmp$beta, cmp$beta_b), 3))) +
  theme_minimal(base_size = 11)
pF <- ggplot(res, aes(beta, -log10(fdr), colour = class)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "standardised effect (High vs Low Braak)", y = "-log10(FDR)") +
  theme_minimal(base_size = 11)

fig_main <- (pC | pD) + plot_annotation(tag_levels = "A")
fig_supp <- ((pA | pB) / (pE | pF)) + plot_annotation(tag_levels = "A")
print(fig_main)
print(fig_supp)

# ---- save everything the commit will need ----
ggsave(file.path(out_fig, "metabolon_lipid_braak_main.png"), fig_main, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(out_fig, "metabolon_lipid_braak_main.pdf"), fig_main, width = 13, height = 6.5)
ggsave(file.path(out_fig, "metabolon_lipid_braak_supp.png"), fig_supp, width = 12, height = 10, dpi = 300)
ggsave(file.path(out_fig, "metabolon_lipid_braak_supp.pdf"), fig_supp, width = 12, height = 10)
write_csv(res,        file.path(out_res, "rosmap_braak_metabolite_effects.csv"))
write_csv(filter_tab, file.path(out_res, "rosmap_metabolite_filter_by_class.csv"))
write_csv(sub_tab,    file.path(out_res, "rosmap_lipid_subpathway_summary.csv"))
write_csv(sub_test,   file.path(out_res, "rosmap_lipid_subpathway_test.csv"))
write_csv(cmp %>% select(feature, class, super_pathway, beta, fdr, beta_b, fdr_b),
          file.path(out_res, "rosmap_braak_metabolite_effects_batch_adjusted.csv"))
write_csv(sens_class, file.path(out_res, "rosmap_braak_coding_sensitivity_class.csv"))
write_csv(sens_focus, file.path(out_res, "rosmap_braak_coding_sensitivity_focus.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_res, "sessionInfo_metabolon_lipid_braak.txt"))
cat("done, files written to", file.path(getwd(), c(out_fig, out_res)), "\n")