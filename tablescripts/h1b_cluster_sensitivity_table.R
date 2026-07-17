############################################################
## Simple H1B cluster random-effect sensitivity table
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

input_path <- here("results", "H1B", "H1B_sensitivity_cluster_random_effect.csv")
output_path <- here("tables", "H1B_sensitivity_cluster_random_effect_table.docx")

if (!file.exists(input_path)) {
  stop("Missing H1B cluster random-effect sensitivity CSV.")
}

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

tbl_h1b_cluster <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  transmute(
    Exposure = X_label,
    Outcome = Y_label,
    Beta = fmt_num(coef, 2),
    `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
    `P value` = fmt_p(pval),
    N = as.integer(N)
  )

ft_h1b_cluster <- flextable(tbl_h1b_cluster) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  align(j = c("Exposure", "Outcome"), align = "left", part = "all") %>%
  align(j = c("Beta", "95% CI", "P value", "N"), align = "center", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding = 2, part = "all") %>%
  width(j = "Exposure", width = 2.3) %>%
  width(j = "Outcome", width = 2.0) %>%
  width(j = "Beta", width = 0.8) %>%
  width(j = "95% CI", width = 1.4) %>%
  width(j = "P value", width = 0.8) %>%
  width(j = "N", width = 0.6) %>%
  fit_to_width(max_width = 8.5)

doc <- read_docx() %>%
  body_add_par("Table. H1B cluster random-effect sensitivity analysis", style = "heading 1") %>%
  body_add_flextable(ft_h1b_cluster)

print(doc, target = output_path)

cat("H1B cluster random-effect table saved to:", output_path, "\n")
