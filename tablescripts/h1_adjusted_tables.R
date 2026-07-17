############################################################
## Publication-ready adjusted-only H1 tables
## Separate H1A and H1B outputs
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

paths <- list(
  h1a_adj = here("results", "H1A Final", "H1Aadj_results.csv"),
  h1b_adj = here("results", "H1B", "h1badj_results.csv")
)

output_h1a <- here("tables", "H1A_adjusted_results_table.docx")
output_h1b <- here("tables", "H1B_adjusted_results_table.docx")

missing_files <- paths[!vapply(paths, file.exists, logical(1))]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required adjusted H1 result files:",
      paste(unname(missing_files), collapse = ", ")
    )
  )
}

h1a_meta <- tribble(
  ~X, ~Exposure,
  "logFER_inf", "Ferritin*",
  "hemo_eeel", "Hemoglobin",
  "hepcidin_eeel", "Hepcidin",
  "logsTFR_inf", "Soluble transferrin receptor*",
  "logRBP_inf", "Retinol-binding protein*",
  "logB12", "Vitamin B12*",
  "logfolate", "Folate*"
)

h1b_meta <- tribble(
  ~X, ~Exposure,
  "lowfer_inf", "Low ferritin (<12 ug/L)",
  "highsTFR", "High soluble transferrin receptor (>8.3 mg/L)",
  "irondefanemia_inf", "Iron deficiency anemia (hemoglobin <11 g/dL and low ferritin or high sTfR)",
  "vitadef_inf", "Vitamin A deficiency (RBP <0.83 umol/L)",
  "folatedef", "Folate deficiency (<10 nmol/L)",
  "B12def", "Vitamin B12 deficiency (<150 pmol/L)",
  "B12marg", "Marginal vitamin B12 (150 to <221 pmol/L)"
)

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

clean_h1a_table <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    select(X, Y, coef, lb, ub, pval, N) %>%
    left_join(h1a_meta, by = "X") %>%
    mutate(
      Outcome = "Telomere length z-score",
      Beta = fmt_num(coef, 2),
      `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
      `P value` = fmt_p(pval),
      N = as.integer(N)
    ) %>%
    select(Exposure, Outcome, Beta, `95% CI`, `P value`, N)
}

clean_h1b_table <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    select(X, Y, coef, lb, ub, pval, N) %>%
    left_join(h1b_meta, by = "X") %>%
    mutate(
      Outcome = "Telomere length z-score",
      Beta = fmt_num(coef, 2),
      `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
      `P value` = fmt_p(pval),
      N = as.integer(N)
    ) %>%
    select(Exposure, Outcome, Beta, `95% CI`, `P value`, N)
}

make_ft <- function(df, widths) {
  ft <- flextable(df) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    align(j = c(1, 2), align = "left", part = "all") %>%
    align(j = c(3, 4, 5, 6), align = "center", part = "all") %>%
    fontsize(size = 8.8, part = "all") %>%
    padding(padding = 2, part = "all") %>%
    line_spacing(space = 0.95, part = "all")

  for (nm in names(widths)) {
    ft <- width(ft, j = nm, width = widths[[nm]])
  }

  ft %>% fit_to_width(max_width = 9.5)
}

h1a_tbl <- clean_h1a_table(paths$h1a_adj)
h1b_tbl <- clean_h1b_table(paths$h1b_adj)

ft_h1a <- make_ft(
  h1a_tbl,
  widths = c(
    Exposure = 2.4,
    Outcome = 2.0,
    Beta = 0.8,
    `95% CI` = 1.3,
    `P value` = 0.8,
    N = 0.6
  )
) %>%
  add_footer_lines(
    values = c(
      "Adjusted associations of micronutrient concentrations with telomere length z-score at 28 months (H1A).",
      "* Log-transformed exposure. Beta estimates reflect a 1-unit increase in the modeled log-scale biomarker."
    )
  )

ft_h1b <- make_ft(
  h1b_tbl,
  widths = c(
    Exposure = 2.8,
    Outcome = 2.0,
    Beta = 0.8,
    `95% CI` = 1.3,
    `P value` = 0.8,
    N = 0.6
  )
) %>%
  add_footer_lines(
    values = c(
      "Adjusted associations of micronutrient deficiency indicators with telomere length z-score at 28 months (H1B)."
    )
  )

doc_h1a <- read_docx() %>%
  body_add_par("Table. H1A adjusted analysis", style = "heading 1") %>%
  body_add_flextable(ft_h1a)

doc_h1b <- read_docx() %>%
  body_add_par("Table. H1B adjusted analysis", style = "heading 1") %>%
  body_add_flextable(ft_h1b)

print(doc_h1a, target = output_h1a)
print(doc_h1b, target = output_h1b)

cat("H1A adjusted table saved to:", output_h1a, "\n")
cat("H1B adjusted table saved to:", output_h1b, "\n")
