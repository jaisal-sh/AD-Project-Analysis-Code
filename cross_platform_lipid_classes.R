# cross_platform_lipid_classes.R
# my goal: do the two ROSMAP brain platforms agree on Braak effects, class by class?
# lipidomics (DLPFC) vs Metabolon HD4 (brain, region not confirmed); both z-scored effects,
# High vs Low Braak, adjusted for age and sex
# caveats: participants probably overlap, class mapping is inexact, species within a class
# are correlated, so intervals and the correlation are descriptive only

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
set.seed(1)
#Gridlines removal
theme_minimal <- function(...) {
  ggplot2::theme_minimal(...) +
    theme(panel.grid = element_blank(),
          axis.line  = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"))
}

# ---- CONFIG ----
proj <- "C:/Users/jayvs/Desktop/ad_vascular_multiomics"
setwd(proj)
metab_res_file <- "C:/Users/jayvs/Desktop/ROSMAP Multi-Omics/results/rosmap_braak_metabolite_effects.csv"
dict_file <- "C:/Users/jayvs/Desktop/ROSMAP dataset/ROSMAP Metabolomics/Metabolon/QC/ROSMAP Metabolon HD4 Data Dictionary.csv"
lipid_res_file <- "results/rosmap_lipidomics_braak_adjusted.csv"
out_fig <- "figures"
out_res <- "results"
dir.create(file.path(out_res, "inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, showWarnings = FALSE)
cat("working in:", getwd(), "\n")

# keep a copy of the metabolon feature-level results here so this repo can rerun the figure
file.copy(metab_res_file, file.path(out_res, "inputs/metabolon_braak_effects.csv"), overwrite = TRUE)

# ---- load both result tables ----
metab <- read_csv(metab_res_file, col_types = cols(feature = col_character(), .default = col_guess()))
stopifnot(all(c("feature", "beta", "fdr", "sub_pathway") %in% names(metab)))
names_key <- read_csv(dict_file, show_col_types = FALSE) %>%
  transmute(feature = as.character(CHEM_ID), chem_name = CHEMICAL_NAME)
metab <- metab %>% left_join(names_key, by = "feature")

lipd <- read_csv(lipid_res_file, show_col_types = FALSE) %>% filter(tissue == "Brain")
stopifnot(nrow(lipd) > 0, all(c("species", "beta_z", "class") %in% names(lipd)))
cat("metabolon features:", nrow(metab), "| lipidomics brain species:", nrow(lipd), "\n")

# metabolon lysophospholipids include LPC etc., so keep only the GPE (LPE) species to match lipidomics LPE
metab <- metab %>% mutate(group = case_when(
  sub_pathway == "Lysophospholipid" & grepl("GPE \\(", chem_name) ~ "Lysophospholipid (GPE only)",
  sub_pathway == "Plasmalogen"      & grepl("GPE \\(", chem_name) ~ "Plasmalogen (GPE only)",
  TRUE ~ sub_pathway))

# ---- class mapping (lipidomics class <-> metabolon group) ----
cls_map <- tribble(
  ~class_lipidomics, ~metabolon_group,                  ~label,
  "PE",              "Phosphatidylethanolamine (PE)",   "PE",
  "PC",              "Phosphatidylcholine (PC)",        "PC",
  "PS",              "Phosphatidylserine (PS)",         "PS",
  "PI",              "Phosphatidylinositol (PI)",       "PI",
  "SM",              "Sphingomyelins",                  "Sphingomyelins",
  "Cer",             "Ceramides",                       "Ceramides",
  "Hex1Cer",         "Hexosylceramides (HCER)",         "Hexosylceramides",
  "Hex2Cer",         "Lactosylceramides (LCER)",        "Lactosylceramides",
  "plasmenylPE",     "Plasmalogen (GPE only)",          "Plasmalogens (PE)",
  "LPE",             "Lysophospholipid (GPE only)",     "LPE",
  "DG",              "Diacylglycerol",                  "DG"
)

lipd_cls  <- lipd  %>% inner_join(cls_map, by = c("class" = "class_lipidomics"))
metab_cls <- metab %>% inner_join(cls_map, by = c("group" = "metabolon_group"))

# bootstrap interval over species for the class median (descriptive: species are correlated)
boot_ci <- function(x, B = 2000) {
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  m <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(m, c(0.025, 0.975)))
}
summ <- function(d, val, platform) {
  d %>% group_by(label) %>%
    summarise(n = n(), med = median(.data[[val]]),
              lo = boot_ci(.data[[val]])[1], hi = boot_ci(.data[[val]])[2],
              n_fdr = sum(fdr < 0.05), .groups = "drop") %>%
    mutate(platform = platform)
}
pl_a <- "Lipidomics (brain, DLPFC)"
pl_b <- "Metabolon (brain)"
both <- bind_rows(summ(lipd_cls %>% mutate(fdr = fdr), "beta_z", pl_a),
                  summ(metab_cls, "beta", pl_b))

# keep classes with at least 3 species on both platforms
keep_cls <- both %>% group_by(label) %>%
  filter(n_distinct(platform) == 2, min(n) >= 3) %>% pull(label) %>% unique()
both <- both %>% filter(label %in% keep_cls)

# check 1: what got matched and how many species per class per platform
print(both %>% select(label, platform, n, n_fdr, med) %>% arrange(label, platform), n = Inf)
cat("classes compared:", length(keep_cls), "\n")
stopifnot("plasmalogen GPE fix is active" = "Plasmalogens (PE)" %in% both$label)
stopifnot(length(keep_cls) >= 5)   # fewer than this and the comparison is not worth plotting

wide <- both %>% select(label, platform, med) %>%
  pivot_wider(names_from = platform, values_from = med) %>%
  rename(lipidomics = all_of(pl_a), metabolon = all_of(pl_b)) %>%
  mutate(same_sign = sign(lipidomics) == sign(metabolon))

# check 2: which classes agree in direction
print(wide %>% arrange(lipidomics), n = Inf)
rho <- suppressWarnings(cor.test(wide$lipidomics, wide$metabolon, method = "spearman"))
cat("same direction in", sum(wide$same_sign), "of", nrow(wide), "classes",
    "| Spearman rho =", round(unname(rho$estimate), 2), "(descriptive, n classes =", nrow(wide), ")\n")

# ---- panels ----
pcols <- setNames(c("#2C7FB8", "#D95F0E"), c(pl_a, pl_b))
ord <- wide %>% arrange(lipidomics) %>% pull(label)
long <- both %>% mutate(label = factor(label, levels = ord))

pA <- ggplot(long, aes(med, label, colour = platform)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_line(aes(group = label), colour = "grey75") +
  geom_linerange(aes(xmin = lo, xmax = hi), alpha = 0.45) +
  geom_point(aes(size = n)) +
  scale_colour_manual(values = pcols, name = NULL) +
  scale_size_continuous(range = c(2, 5), name = "species in class") +
  labs(x = "median standardised effect (High vs Low Braak)", y = NULL,
       caption = "bars = 95% bootstrap interval over species (descriptive)\nclasses with >= 3 species on both platforms") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical", plot.caption = element_text(hjust = 0))

pB <- ggplot(wide, aes(lipidomics, metabolon, label = label, colour = same_sign)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 2.2) +
  geom_text(hjust = -0.12, vjust = 0.4, size = 3, check_overlap = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "#D95F0E"),
                      labels = c("TRUE" = "same direction", "FALSE" = "opposite direction"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "lipidomics, median standardised effect", y = "Metabolon, median standardised effect",
       caption = paste0("Spearman rho = ", round(unname(rho$estimate), 2), " across ", nrow(wide),
                        " classes (descriptive)")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig_cross <- (pA | pB) + plot_annotation(tag_levels = "A")
print(fig_cross)

# ---- save for the commit ----
ggsave(file.path(out_fig, "cross_platform_lipid_classes.png"), fig_cross, width = 12, height = 6, dpi = 300)
ggsave(file.path(out_fig, "cross_platform_lipid_classes.pdf"), fig_cross, width = 12, height = 6)
write_csv(both, file.path(out_res, "cross_platform_lipid_class_medians.csv"))
write_csv(wide, file.path(out_res, "cross_platform_lipid_class_agreement.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_res, "sessionInfo_cross_platform_lipid_classes.txt"))
cat("done, files written to", file.path(getwd(), c(out_fig, out_res)), "\n")