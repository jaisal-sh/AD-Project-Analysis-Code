# lipidomics_braak_adjusted.R
# my goal: redo the brain and plasma lipidomics Braak comparison with the same model as the
# Metabolon figure (age + sex, BH FDR), one row per person, and summarise by lipid class
# brain = DLPFC (pos + neg combined), plasma = negative ionization only, so classes that
# need positive mode (like LPC) may be missing from plasma for technical reasons
# plasma had 232 people with two technical replicate rows, averaged per person below
# brain batch is unevenly spread across Braak groups, so i report the batch-adjusted model too

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})
set.seed(1)

# ---- CONFIG ----
proj <- "C:/Users/jayvs/Desktop/ad_vascular_multiomics"
setwd(proj)
wb <- "processed/M2OVE.xlsx"
out_fig <- "figures"
out_res <- "results"
dir.create(out_fig, showWarnings = FALSE)
dir.create(out_res, showWarnings = FALSE)
cat("working in:", getwd(), "\n")

# ---- clinical ----
clin <- read_excel(wb, sheet = "rosmap_lipid_clinical") %>%
  select(individualID, msex, age_death, braaksc) %>% distinct()
stopifnot(!any(duplicated(clin$individualID)))   # one row per person or the join double counts

# ---- load one tissue, average technical replicates per person, join clinical, drop standards ----
prep_tissue <- function(sheet) {
  x <- read_excel(wb, sheet = sheet)
  num_cols <- setdiff(names(x)[sapply(x, is.numeric)], c("sequence order", "injection order", "batch"))
  n_rows <- nrow(x)
  dup_ids <- unique(x$individualID[duplicated(x$individualID)])
  
  if (length(dup_ids) > 0) {
    # how well do the two rows per person agree, across every species
    two <- x %>% filter(individualID %in% dup_ids) %>%
      group_by(individualID) %>% filter(n() == 2) %>% mutate(rep = row_number()) %>% ungroup()
    a <- two %>% filter(rep == 1) %>% arrange(individualID)
    b <- two %>% filter(rep == 2) %>% arrange(individualID)
    rc <- suppressWarnings(map_dbl(num_cols, function(cn)
      tryCatch(cor(a[[cn]], b[[cn]], use = "complete.obs"), error = function(e) NA_real_)))
    cat("\n", sheet, "| repeated people:", length(dup_ids),
        "| median replicate correlation, all species:", round(median(rc, na.rm = TRUE), 3),
        "| 10th percentile:", round(quantile(rc, 0.1, na.rm = TRUE), 3), "\n")
    # both rows of a person must share a batch or averaging would hide a batch difference
    nb <- x %>% filter(individualID %in% dup_ids) %>% group_by(individualID) %>%
      summarise(nb = n_distinct(batch), .groups = "drop")
    stopifnot(all(nb$nb == 1))
    # average on the log2 scale, one row per person
    x <- x %>% group_by(individualID) %>%
      summarise(batch = first(batch),
                across(all_of(num_cols), ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
                .groups = "drop")
  }
  stopifnot(!any(duplicated(x$individualID)))   # one row per person from here on
  cat(sheet, "| rows before:", n_rows, "| people after:", nrow(x), "\n")
  
  dat <- x %>%
    inner_join(clin, by = "individualID") %>%
    mutate(braak = as.numeric(braaksc),
           age   = as.numeric(gsub("\\+", "", as.character(age_death))),   # handles "90+"
           sex   = factor(msex),
           batch = factor(batch)) %>%
    filter(!is.na(braak), !is.na(age), !is.na(sex)) %>%
    mutate(braak_group = factor(ifelse(braak >= 4, "High", "Low"), levels = c("Low", "High")))
  cand <- num_cols
  std  <- cand[grepl("-d[0-9]+", cand)]          # deuterated internal standards, not lipids
  cand <- setdiff(cand, std)
  keep <- cand[sapply(cand, function(f) mean(!is.na(dat[[f]])) >= 0.8)]
  cat("  subjects:", nrow(dat),
      "| High:", sum(dat$braak_group == "High"), "| Low:", sum(dat$braak_group == "Low"), "\n")
  cat("  species before filter:", length(cand), "| kept (>= 80% non-missing):", length(keep), "\n")
  cat("  dropped as standards:", ifelse(length(std) == 0, "none", paste(std, collapse = ", ")), "\n")
  cat("  batch levels:", nlevels(droplevels(dat$batch)), "\n")
  list(dat = dat, species = keep)
}

brain  <- prep_tissue("rosmap_lipid_brain")
plasma <- prep_tissue("rosmap_lipid_plasma")

# ---- model: species ~ Braak group + age + sex. values are already log2, so the raw-scale
# coefficient is the adjusted log2 fold change; the z-scored one is comparable to Metabolon
fit_species <- function(p, with_batch = FALSE) {
  use_batch <- with_batch && nlevels(droplevels(p$dat$batch)) >= 2   # lm fails on a one-level factor
  rhs <- if (use_batch) "braak_group + age + sex + batch" else "braak_group + age + sex"
  map_dfr(p$species, function(f) {
    d <- p$dat %>% mutate(y_raw = .data[[f]], y_z = as.numeric(scale(.data[[f]])), batch = droplevels(batch))
    sr <- summary(lm(as.formula(paste("y_raw ~", rhs)), data = d))$coefficients["braak_groupHigh", ]
    sz <- summary(lm(as.formula(paste("y_z ~", rhs)), data = d))$coefficients["braak_groupHigh", ]
    tt <- t.test(d$y_raw[d$braak_group == "High"], d$y_raw[d$braak_group == "Low"])
    tibble(species = f, log2fc_adj = sr["Estimate"], beta_z = sz["Estimate"],
           p = sz["Pr(>|t|)"], p_unadj_ttest = tt$p.value)
  }) %>% mutate(fdr = p.adjust(p, method = "BH"))
}

get_class <- function(s) str_extract(s, "^[A-Za-z0-9]+")   # "Cer d16:1/18:0" and "Cer(d18:1_16:0)" both give "Cer"

res_all <- bind_rows(
  fit_species(brain)  %>% mutate(tissue = "Brain"),
  fit_species(plasma) %>% mutate(tissue = "Plasma")
) %>% mutate(class = get_class(species))

# check: adjusted vs the old unadjusted t-test counts (old brain 33/92, plasma 13/117 on 662 rows)
print(res_all %>% group_by(tissue) %>%
        summarise(species = n(),
                  nominal_ttest_p05 = sum(p_unadj_ttest < 0.05),
                  fdr05_adjusted = sum(fdr < 0.05),
                  pct_down_among_fdr_sig = round(100 * mean(log2fc_adj[fdr < 0.05] < 0))))

# ---- batch sensitivity (the files are already batch corrected, so i expect little change) ----
res_b <- bind_rows(
  fit_species(brain,  with_batch = TRUE) %>% mutate(tissue = "Brain"),
  fit_species(plasma, with_batch = TRUE) %>% mutate(tissue = "Plasma")
) %>% select(tissue, species, beta_z_b = beta_z, fdr_b = fdr)
cmp <- res_all %>% inner_join(res_b, by = c("tissue", "species"))
print(cmp %>% group_by(tissue) %>%
        summarise(beta_corr = round(cor(beta_z, beta_z_b), 3),
                  sig_before = sum(fdr < 0.05), sig_after = sum(fdr_b < 0.05)))

# ---- direction summary and batch balance (the numbers i quote in the legend) ----
dir_sum <- cmp %>% group_by(tissue) %>%
  summarise(n = n(), n_neg = sum(beta_z < 0), pct_neg_main = round(100 * mean(beta_z < 0)),
            pct_neg_batch_adj = round(100 * mean(beta_z_b < 0)),
            median_beta_main = round(median(beta_z), 3),
            median_beta_batch_adj = round(median(beta_z_b), 3), .groups = "drop")
print(dir_sum)
batch_tab <- bind_rows(
  as.data.frame(table(batch = brain$dat$batch,  braak = brain$dat$braak_group))  %>% mutate(tissue = "Brain"),
  as.data.frame(table(batch = plasma$dat$batch, braak = plasma$dat$braak_group)) %>% mutate(tissue = "Plasma"))
print(batch_tab %>% pivot_wider(names_from = braak, values_from = Freq) %>% arrange(tissue, batch), n = Inf)

# ---- class summary ----
class_sum <- res_all %>% group_by(tissue, class) %>%
  summarise(n = n(), n_sig = sum(fdr < 0.05), pct_down = round(100 * mean(log2fc_adj < 0)),
            med_log2fc = round(median(log2fc_adj), 3), .groups = "drop") %>%
  filter(n >= 3)
print(class_sum %>% arrange(tissue, med_log2fc), n = Inf)

# ---- panels ----
tcols <- c("Brain" = "#2C7FB8", "Plasma" = "#D95F0E")
volc <- function(d, ttl) {
  ggplot(d, aes(log2fc_adj, -log10(p), colour = fdr < 0.05)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "grey70"), name = "FDR < 0.05") +
    labs(title = ttl, x = "adjusted log2 fold change (High vs Low Braak)", y = "-log10(p)") +
    theme_minimal(base_size = 11)
}
pA <- volc(res_all %>% filter(tissue == "Brain"),  "Brain lipids (DLPFC)")
pB <- volc(res_all %>% filter(tissue == "Plasma"), "Plasma lipids (negative ionization)")

pC <- class_sum %>%
  mutate(class = fct_reorder(class, med_log2fc, .fun = median)) %>%
  ggplot(aes(med_log2fc, class, colour = tissue, size = n)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(position = position_dodge(width = 0.6), alpha = 0.85) +
  geom_text(aes(label = paste0(n_sig, "/", n)), position = position_dodge(width = 0.6),
            hjust = -0.7, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = tcols, name = NULL) +
  scale_size_continuous(range = c(2, 6), name = "species in class") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "median adjusted log2 fold change (High vs Low Braak)", y = NULL,
       caption = "labels = FDR-significant / species tested; classes with >= 3 species; adjusted for age and sex") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

fig_lipid <- (pA | pB) / pC + plot_layout(heights = c(1, 1.3)) + plot_annotation(tag_levels = "A")
print(fig_lipid)

# ---- save for the commit ----
ggsave(file.path(out_fig, "lipidomics_braak_adjusted_main.png"), fig_lipid, width = 11, height = 10, dpi = 300)
ggsave(file.path(out_fig, "lipidomics_braak_adjusted_main.pdf"), fig_lipid, width = 11, height = 10)
write_csv(res_all,   file.path(out_res, "rosmap_lipidomics_braak_adjusted.csv"))
write_csv(class_sum, file.path(out_res, "rosmap_lipidomics_braak_class_summary.csv"))
write_csv(cmp %>% select(tissue, species, class, beta_z, fdr, beta_z_b, fdr_b),
          file.path(out_res, "rosmap_lipidomics_braak_batch_sensitivity.csv"))
write_csv(dir_sum,   file.path(out_res, "rosmap_lipidomics_braak_direction_summary.csv"))
write_csv(batch_tab, file.path(out_res, "rosmap_lipidomics_braak_batch_by_group.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_res, "sessionInfo_lipidomics_braak_adjusted.txt"))
cat("done, files written to", file.path(getwd(), c(out_fig, out_res)), "\n")