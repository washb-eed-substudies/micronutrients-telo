############################################################
## Publication-ready H2B co-deficiency tables
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

input_path <- here("results", "H2B", "H2B_pair_prevalence_fisher.csv")
output_path <- here("tables", "H2B_tables.docx")
legacy_output_path <- here("tables", "H2B_co_deficiency_table.docx")

if (!file.exists(input_path)) {
  stop("Missing H2B results CSV. Run H2B_Final.Rmd first.")
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

sig_label <- function(raw_sig, fdr_sig) {
  case_when(
    !is.na(fdr_sig) & fdr_sig != "" ~ paste0(fdr_sig, " FDR"),
    !is.na(raw_sig) & raw_sig != "" ~ paste0(raw_sig, " raw"),
    TRUE ~ "NS"
  )
}

pair_title <- function(x) {
  x %>%
    str_replace_all(" deficiency anemia", " deficiency anemia") %>%
    str_replace_all("Vitamin A deficiency", "Vitamin A deficiency") %>%
    str_replace_all("High sTfR", "High sTfR") %>%
    str_replace_all("Low ferritin", "Low ferritin")
}

h2b_raw <- readr::read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    pair_label = pair_title(pair_label),
    sig_display = sig_label(sig_raw, sig_fdr)
  ) %>%
  arrange(p_value, desc(obs_exp_ratio))

h2b_main_tbl <- h2b_raw %>%
  transmute(
    `Co-deficiency pair` = pair_label,
    `Pairwise complete N` = as.integer(n_complete),
    `Co-deficient n (%)` = paste0(as.integer(n_both), " (", fmt_pct(prev_both, 1), ")"),
    `Expected co-deficiency (%)` = fmt_pct(expected_both_indep, 1),
    `Observed/expected ratio` = fmt_num(obs_exp_ratio, 2),
    `Odds ratio (95% CI)` = paste0(
      fmt_num(odds_ratio, 2), " (", fmt_num(or_lb, 2), ", ", fmt_num(or_ub, 2), ")"
    ),
    `P value` = fmt_p(p_value),
    `FDR-adjusted P` = fmt_p(p_fdr),
    Significance = sig_display
  )

h2b_marginal_tbl <- h2b_raw %>%
  transmute(
    `Co-deficiency pair` = pair_label,
    `Pairwise complete N` = as.integer(n_complete),
    `Deficiency 1 n (%)` = paste0(as.integer(n_def1), " (", fmt_pct(prev_def1, 1), ")"),
    `Deficiency 2 n (%)` = paste0(as.integer(n_def2), " (", fmt_pct(prev_def2, 1), ")"),
    `Co-deficient n (%)` = paste0(as.integer(n_both), " (", fmt_pct(prev_both, 1), ")")
  )

none_fdr_sig <- all(is.na(h2b_raw$sig_fdr) | h2b_raw$sig_fdr == "")

ft_main <- flextable(h2b_main_tbl) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  autofit() %>%
  align(j = "Co-deficiency pair", align = "left", part = "all") %>%
  align(j = 2:ncol(h2b_main_tbl), align = "center", part = "all") %>%
  fontsize(size = 8.5, part = "all") %>%
  padding(padding = 2, part = "all") %>%
  width(j = "Co-deficiency pair", width = 2.8) %>%
  width(j = "Pairwise complete N", width = 0.8) %>%
  width(j = "Co-deficient n (%)", width = 1.0) %>%
  width(j = "Expected co-deficiency (%)", width = 1.1) %>%
  width(j = "Observed/expected ratio", width = 1.0) %>%
  width(j = "Odds ratio (95% CI)", width = 1.5) %>%
  width(j = "P value", width = 0.8) %>%
  width(j = "FDR-adjusted P", width = 0.9) %>%
  width(j = "Significance", width = 0.8) %>%
  fit_to_width(max_width = 9.5) %>%
  add_footer_lines(
    values = c(
      "Primary H2B pairwise co-deficiency results based on Fisher's exact test.",
      "All analyses were restricted to children with non-missing telomere length z-score.",
      "Pairwise complete N = children within the telomere-complete analytic sample with non-missing values for both deficiencies in a given pair.",
      "Observed/expected ratio >1 indicates co-deficiency more frequent than expected under independence.",
      if (none_fdr_sig) {
        "No co-deficiency pair remained statistically significant after FDR correction."
      } else {
        "Significance labels prioritize FDR-adjusted significance when present; otherwise they show raw significance."
      }
    )
  )

ft_marginal <- flextable(h2b_marginal_tbl) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  autofit() %>%
  align(j = "Co-deficiency pair", align = "left", part = "all") %>%
  align(j = 2:ncol(h2b_marginal_tbl), align = "center", part = "all") %>%
  fontsize(size = 8.5, part = "all") %>%
  padding(padding = 2, part = "all") %>%
  width(j = "Co-deficiency pair", width = 3.0) %>%
  fit_to_width(max_width = 9.5) %>%
  add_footer_lines(
    values = c(
      "Marginal deficiency counts can overlap because co-deficient children contribute to both component deficiency totals."
    )
  )

doc <- read_docx() %>%
  body_add_par("Table 1. H2B pairwise co-deficiency analysis", style = "heading 1") %>%
  body_end_section_landscape() %>%
  body_add_flextable(ft_main) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Table 2. H2B component deficiency counts within each pair", style = "heading 1") %>%
  body_add_flextable(ft_marginal)

print(doc, target = output_path)
print(doc, target = legacy_output_path)

cat("H2B tables saved to:", output_path, "\n")
