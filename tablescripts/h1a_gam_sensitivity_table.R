############################################################
## Publication-ready H1A RE-GAM sensitivity table
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

input_path <- here("results", "H1A Final", "sensitivity", "H1A_RE_GAM_sensitivity.csv")
contrast_path <- here("results", "H1A Final", "sensitivity", "H1A_RE_GAM_percentile_contrasts.csv")
output_path <- here("tables", "H1A_RE_GAM_sensitivity_table.docx")

if (!file.exists(input_path)) {
  stop("Missing H1A RE-GAM sensitivity CSV. Run the H1A sensitivity analysis first.")
}

if (!file.exists(contrast_path)) {
  stop("Missing H1A percentile-contrast CSV. Re-run the updated H1A sensitivity analysis in H1A_Final.Rmd to generate H1A_RE_GAM_percentile_contrasts.csv.")
}

fmt_num <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits)
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

sig_label <- function(x) {
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

gam_tbl <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  filter(model %in% c("adjusted", "unadjusted")) %>%
  filter(term_type == "smooth") %>%
  mutate(
    Model = case_when(
      model == "adjusted" ~ "Adjusted",
      model == "unadjusted" ~ "Unadjusted",
      TRUE ~ model
    ),
    Exposure = case_when(
      X_label == "Ferritin (µg/L)" ~ "Ferritin",
      X_label == "Hemoglobin (g/dL)" ~ "Hemoglobin",
      X_label == "Hepcidin (ng/mL)" ~ "Hepcidin",
      X_label == "sTfR (mg/L)" ~ "sTfR",
      X_label == "RBP (µmol/L)" ~ "RBP",
      X_label == "Vitamin B12 (pmol/L)" ~ "Vitamin B12",
      X_label == "Folate (nmol/L)" ~ "Folate",
      TRUE ~ X_label
    ),
    EDF = fmt_num(edf, 2),
    `Ref df` = fmt_num(ref.df, 2),
    Statistic = fmt_num(statistic, 2),
    `Shape p` = fmt_p(pval),
    Sig = sig_label(pval),
    N = as.integer(as.numeric(N))
  ) %>%
  select(
    Model, Exposure, EDF, `Ref df`,
    Statistic, `Shape p`, Sig, N
  )

contrast_tbl <- readr::read_csv(contrast_path, show_col_types = FALSE) %>%
  filter(model %in% c("adjusted", "unadjusted")) %>%
  mutate(
    Model = case_when(
      model == "adjusted" ~ "Adjusted",
      model == "unadjusted" ~ "Unadjusted",
      TRUE ~ model
    ),
    Exposure = case_when(
      X_label == "Ferritin (µg/L)" ~ "Ferritin",
      X_label == "Hemoglobin (g/dL)" ~ "Hemoglobin",
      X_label == "Hepcidin (ng/mL)" ~ "Hepcidin",
      X_label == "sTfR (mg/L)" ~ "sTfR",
      X_label == "RBP (µmol/L)" ~ "RBP",
      X_label == "Vitamin B12 (pmol/L)" ~ "Vitamin B12",
      X_label == "Folate (nmol/L)" ~ "Folate",
      TRUE ~ X_label
    ),
    P10 = fmt_num(q10, 2),
    P90 = fmt_num(q90, 2),
    Diff = fmt_num(point_diff, 2),
    `Diff 95% CI` = paste0(fmt_num(lb_diff, 2), ", ", fmt_num(ub_diff, 2)),
    `Contrast p` = fmt_p(pval_diff),
    `Contrast sig` = sig_label(pval_diff),
    N = as.integer(as.numeric(N))
  ) %>%
  select(
    Model, Exposure, P10, P90,
    Diff, `Diff 95% CI`, `Contrast p`, `Contrast sig`, N
  )

combined_tbl <- gam_tbl %>%
  left_join(
    contrast_tbl,
    by = c("Model", "Exposure", "N")
  ) %>%
  mutate(
    Model = factor(Model, levels = c("Adjusted", "Unadjusted")),
    Exposure = factor(Exposure, levels = sort(unique(Exposure)))
  ) %>%
  arrange(Model, Exposure) %>%
  mutate(
    Model = recode(as.character(Model), Adjusted = "Adj", Unadjusted = "Unadj"),
    Exposure = as.character(Exposure)
  )

ft_combined <- flextable(combined_tbl) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  set_table_properties(layout = "fixed") %>%
  align(j = c("Model", "Exposure"), align = "left", part = "all") %>%
  align(
    j = c("EDF", "Ref df", "Statistic", "Shape p", "Sig",
          "P10", "P90", "Diff", "Diff 95% CI",
          "Contrast p", "Contrast sig", "N"),
    align = "center",
    part = "all"
  ) %>%
  fontsize(size = 6.2, part = "all") %>%
  padding(padding = 0.8, part = "all") %>%
  line_spacing(space = 0.8, part = "all") %>%
  width(j = "Model", width = 0.55) %>%
  width(j = "Exposure", width = 1.0) %>%
  width(j = "EDF", width = 0.42) %>%
  width(j = "Ref df", width = 0.48) %>%
  width(j = "Statistic", width = 0.5) %>%
  width(j = "Shape p", width = 0.48) %>%
  width(j = "Sig", width = 0.4) %>%
  width(j = "P10", width = 0.48) %>%
  width(j = "P90", width = 0.48) %>%
  width(j = "Diff", width = 0.48) %>%
  width(j = "Diff 95% CI", width = 0.8) %>%
  width(j = "Contrast p", width = 0.48) %>%
  width(j = "Contrast sig", width = 0.45) %>%
  width(j = "N", width = 0.35) %>%
  fit_to_width(max_width = 10.2) %>%
  add_footer_lines(
    values = c(
      "Each row combines the shape test and the P90 vs P10 contrast for the same biomarker. EDF near 1 suggests little nonlinearity."
    )
  )

landscape_section <- prop_section(
  page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.3, bottom = 0.3, left = 0.35, right = 0.35)
)

doc <- read_docx()
doc <- body_add_par(doc, "H1A RE-GAM sensitivity analysis", style = "Normal")
doc <- body_end_block_section(doc, block_section(landscape_section))
doc <- body_add_flextable(doc, ft_combined)

print(doc, target = output_path)

cat("H1A RE-GAM sensitivity table saved to:", output_path, "\n")
