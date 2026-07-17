############################################################
## Publication-ready H1 result tables
## H1A and H1B, adjusted and unadjusted
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

paths <- list(
  h1a_adj = here("results", "H1A Final", "H1Aadj_results.csv"),
  h1a_unadj = here("results", "H1A Final", "H1Aunadj_results.csv"),
  h1b_adj = here("results", "H1B", "h1badj_results.csv"),
  h1b_unadj = here("results", "H1B", "h1bunadj_results.csv")
)

output_path <- here("tables", "H1_results_tables.docx")

h1a_labels <- c(
  logFER_inf = "Ferritin",
  hemo_eeel = "Hemoglobin",
  hepcidin_eeel = "Hepcidin",
  logsTFR_inf = "Soluble transferrin receptor",
  logRBP_inf = "Retinol-binding protein",
  logB12 = "Vitamin B12",
  logfolate = "Folate",
  TS_t3_Z = "Telomere length z-score"
)

h1b_labels <- c(
  lowfer_inf = "Low ferritin",
  highsTFR = "High soluble transferrin receptor",
  irondefanemia_inf = "Iron deficiency anemia",
  vitadef_inf = "Vitamin A deficiency",
  folatedef = "Folate deficiency",
  B12def = "Vitamin B12 deficiency",
  B12marg = "Marginal vitamin B12",
  TS_t3_Z = "Telomere length z-score"
)

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

clean_result_table <- function(path, exposure_labels, table_type) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    select(X, Y, coef, lb, ub, pval, N) %>%
    mutate(
      Exposure = recode(X, !!!exposure_labels),
      Outcome = recode(Y, !!!exposure_labels),
      Estimate = fmt_num(coef, 2),
      `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
      `P value` = fmt_p(pval),
      N = as.integer(N)
    ) %>%
    select(Exposure, Outcome, Estimate, `95% CI`, `P value`, N) %>%
    arrange(Exposure) %>%
    mutate(Table = table_type)
}

make_ft <- function(df, title_text) {
  flextable(df %>% select(-Table)) %>%
    set_header_labels(
      Exposure = "Exposure",
      Outcome = "Outcome",
      Estimate = "Beta",
      `95% CI` = "95% CI",
      `P value` = "P value",
      N = "N"
    ) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    align(align = "left", j = c("Exposure", "Outcome"), part = "all") %>%
    align(align = "center", j = c("Estimate", "95% CI", "P value", "N"), part = "all") %>%
    fontsize(size = 9, part = "all") %>%
    padding(padding = 2, part = "all") %>%
    width(j = "Exposure", width = 2.2) %>%
    width(j = "Outcome", width = 1.8) %>%
    width(j = "Estimate", width = 0.8) %>%
    width(j = "95% CI", width = 1.4) %>%
    width(j = "P value", width = 0.8) %>%
    width(j = "N", width = 0.6) %>%
    fit_to_width(max_width = 8.8) %>%
    add_footer_lines(values = title_text)
}

h1a_adj_tbl <- clean_result_table(
  paths$h1a_adj,
  exposure_labels = h1a_labels,
  table_type = "H1A adjusted"
)

h1a_unadj_tbl <- clean_result_table(
  paths$h1a_unadj,
  exposure_labels = h1a_labels,
  table_type = "H1A unadjusted"
)

h1b_adj_tbl <- clean_result_table(
  paths$h1b_adj,
  exposure_labels = h1b_labels,
  table_type = "H1B adjusted"
)

h1b_unadj_tbl <- clean_result_table(
  paths$h1b_unadj,
  exposure_labels = h1b_labels,
  table_type = "H1B unadjusted"
)

ft_h1a_adj <- make_ft(
  h1a_adj_tbl,
  "Adjusted associations of micronutrient concentrations with telomere length z-score (H1A)"
)

ft_h1a_unadj <- make_ft(
  h1a_unadj_tbl,
  "Unadjusted associations of micronutrient concentrations with telomere length z-score (H1A)"
)

ft_h1b_adj <- make_ft(
  h1b_adj_tbl,
  "Adjusted associations of micronutrient deficiency indicators with telomere length z-score (H1B)"
)

ft_h1b_unadj <- make_ft(
  h1b_unadj_tbl,
  "Unadjusted associations of micronutrient deficiency indicators with telomere length z-score (H1B)"
)

doc <- read_docx()
doc <- body_add_par(doc, "H1A and H1B result tables", style = "heading 1")
doc <- body_add_par(doc, "Table 1. H1A adjusted results", style = "heading 2")
doc <- body_add_flextable(doc, ft_h1a_adj)
doc <- body_add_par(doc, "Table 2. H1A unadjusted results", style = "heading 2")
doc <- body_add_flextable(doc, ft_h1a_unadj)
doc <- body_add_par(doc, "Table 3. H1B adjusted results", style = "heading 2")
doc <- body_add_flextable(doc, ft_h1b_adj)
doc <- body_add_par(doc, "Table 4. H1B unadjusted results", style = "heading 2")
doc <- body_add_flextable(doc, ft_h1b_unadj)

print(doc, target = output_path)

cat("H1 results tables saved to:", output_path, "\n")
