############################################################
## Publication-ready H2A result tables
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

paths <- list(
  adj = here("results", "H2Aadj_results.csv"),
  unadj = here("results", "H2Aunadj_results.csv"),
  primary_cat = here("results", "H2A_primary_categorical_model.csv"),
  missingness = here("results", "H2A_missingness_sensitivity.csv")
)

output_path <- here("tables", "H2A_results_tables.docx")
data_path <- "/Users/sumuccu/Desktop/bangladesh-merged-eed-mn-data.csv"

missing_files <- paths[!vapply(paths, file.exists, logical(1))]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required H2A result files:",
      paste(unname(missing_files), collapse = ", ")
    )
  )
}

if (!file.exists(data_path)) {
  stop("Missing underlying dataset needed for H2A domain-size table.")
}

h2a_labels <- c(
  mm_domain_primary = "Micronutrient domain burden score",
  mm_domain_primary_relaxed = "Missing-data sensitivity score",
  mm_noiron = "Non-iron deficiency burden",
  mm_bdef = "B-vitamin deficiency burden",
  mm_aanemia = "Vitamin A and iron deficiency anemia burden",
  TS_t3_Z = "Telomere length z-score"
)

domain_labels <- c(
  irondef_inf = "Iron deficiency domain",
  vit_A_def = "Vitamin A deficiency domain",
  B12def = "Vitamin B12 deficiency domain",
  folatedef = "Folate deficiency domain",
  irondefanemia_inf = "Iron deficiency anemia domain"
)

composite_components <- list(
  "Micronutrient deficiency burden score" = c("irondef_inf", "vit_A_def", "B12def", "folatedef"),
  "Non-iron deficiency burden" = c("vit_A_def", "B12def", "folatedef"),
  "B-vitamin deficiency burden" = c("B12def", "folatedef"),
  "Vitamin A and iron deficiency anemia burden" = c("vit_A_def", "irondefanemia_inf")
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

sig_label <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

fmt_ci <- function(lb, ub, digits = 1, suffix = "") {
  paste0(
    formatC(lb, format = "f", digits = digits),
    ", ",
    formatC(ub, format = "f", digits = digits),
    suffix
  )
}

make_ft <- function(df) {
  ft <- flextable(df) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    autofit() %>%
    fontsize(size = 9, part = "all") %>%
    padding(padding = 2, part = "all") %>%
    fit_to_width(max_width = 8.8)

  left_cols <- seq_len(min(2, ncol(df)))
  center_cols <- setdiff(seq_len(ncol(df)), left_cols)

  ft <- align(ft, j = left_cols, align = "left", part = "all")

  if (length(center_cols) > 0) {
    ft <- align(ft, j = center_cols, align = "center", part = "all")
  }

  ft
}

h2a_raw <- read.csv(data_path)

# Restrict all H2A descriptive and result tables to the telomere-available analytic sample.
h2a_raw <- h2a_raw %>%
  filter(!is.na(TS_t3_Z))

make_def_count <- function(dat, vars, min_nonmissing = length(vars)) {
  observed_n <- rowSums(!is.na(dat[, vars]))
  score <- rowSums(dat[, vars], na.rm = TRUE)
  score[observed_n < min_nonmissing] <- NA_real_
  score
}

h2a_raw$mm_domain_primary <- make_def_count(
  h2a_raw,
  vars = c("irondef_inf", "vit_A_def", "B12def", "folatedef")
)
h2a_raw$mm_domain_primary_relaxed <- make_def_count(
  h2a_raw,
  vars = c("irondef_inf", "vit_A_def", "B12def", "folatedef"),
  min_nonmissing = 3
)
h2a_raw$mm_noiron <- make_def_count(
  h2a_raw,
  vars = c("vit_A_def", "B12def", "folatedef")
)
h2a_raw$mm_bdef <- make_def_count(
  h2a_raw,
  vars = c("B12def", "folatedef")
)
h2a_raw$mm_aanemia <- make_def_count(
  h2a_raw,
  vars = c("vit_A_def", "irondefanemia_inf")
)
h2a_raw$mm_domain_primary_cat <- case_when(
  is.na(h2a_raw$mm_domain_primary) ~ NA_character_,
  h2a_raw$mm_domain_primary == 0 ~ "0 deficiencies",
  h2a_raw$mm_domain_primary == 1 ~ "1 deficiency",
  h2a_raw$mm_domain_primary >= 2 ~ "2+ deficiencies"
)

binom_ci_pct <- function(num, denom) {
  ci <- stats::binom.test(num, denom)$conf.int * 100
  c(lb = ci[1], ub = ci[2])
}

domain_size_tbl <- tibble(
  variable = c("irondef_inf", "vit_A_def", "B12def", "folatedef", "irondefanemia_inf")
) %>%
  mutate(
    Domain = domain_labels[variable],
    `Non-missing domain N` = vapply(variable, function(v) sum(!is.na(h2a_raw[[v]])), numeric(1)),
    `Deficient n (%)` = vapply(
      variable,
      function(v) paste0(
        sum(h2a_raw[[v]] == 1, na.rm = TRUE),
        " (",
        formatC(100 * mean(h2a_raw[[v]] == 1, na.rm = TRUE), format = "f", digits = 1),
        "%)"
      ),
      character(1)
    ),
    `Deficient % 95% CI` = vapply(
      variable,
      function(v) {
        denom <- sum(!is.na(h2a_raw[[v]]))
        num <- sum(h2a_raw[[v]] == 1, na.rm = TRUE)
        ci <- binom_ci_pct(num, denom)
        fmt_ci(ci["lb"], ci["ub"], digits = 1, suffix = "%")
      },
      character(1)
    ),
    `Domain N with telomere outcome` = vapply(
      variable,
      function(v) sum(!is.na(h2a_raw[[v]]) & !is.na(h2a_raw$TS_t3_Z)),
      numeric(1)
    ),
    `Deficient n with outcome (%)` = vapply(
      variable,
      function(v) {
        denom <- sum(!is.na(h2a_raw[[v]]) & !is.na(h2a_raw$TS_t3_Z))
        num <- sum(h2a_raw[[v]] == 1 & !is.na(h2a_raw$TS_t3_Z), na.rm = TRUE)
        paste0(
          num,
          " (",
          formatC(100 * num / denom, format = "f", digits = 1),
          "%)"
        )
      },
      character(1)
    ),
    `Deficient % with outcome 95% CI` = vapply(
      variable,
      function(v) {
        denom <- sum(!is.na(h2a_raw[[v]]) & !is.na(h2a_raw$TS_t3_Z))
        num <- sum(h2a_raw[[v]] == 1 & !is.na(h2a_raw$TS_t3_Z), na.rm = TRUE)
        ci <- binom_ci_pct(num, denom)
        fmt_ci(ci["lb"], ci["ub"], digits = 1, suffix = "%")
      },
      character(1)
    )
  ) %>%
  select(
    Domain,
    `Non-missing domain N`,
    `Deficient n (%)`,
    `Deficient % 95% CI`,
    `Domain N with telomere outcome`,
    `Deficient n with outcome (%)`,
    `Deficient % with outcome 95% CI`
  )

component_domain_tbl <- purrr::imap_dfr(
  composite_components,
  function(vars, composite_name) {
    composite_complete_n <- sum(stats::complete.cases(h2a_raw[, vars]))

    tibble(variable = vars) %>%
      mutate(
        Composite = composite_name,
        Domain = domain_labels[variable],
        `Domain non-missing N` = vapply(variable, function(v) sum(!is.na(h2a_raw[[v]])), numeric(1)),
        `Deficient n (%)` = vapply(
          variable,
          function(v) paste0(
            sum(h2a_raw[[v]] == 1, na.rm = TRUE),
            " (",
            formatC(100 * mean(h2a_raw[[v]] == 1, na.rm = TRUE), format = "f", digits = 1),
            "%)"
          ),
          character(1)
        ),
        `Complete composite N` = composite_complete_n
      ) %>%
      select(Composite, Domain, `Domain non-missing N`, `Deficient n (%)`, `Complete composite N`)
  }
)

composite_n_tbl <- readr::read_csv(paths$adj, show_col_types = FALSE) %>%
  transmute(
    Composite = h2a_labels[X],
    `Model N` = as.integer(N)
  )

summarize_burden_levels <- function(dat, score_var, max_domains, label, outcome_var = "TS_t3_Z") {
  score_sym <- rlang::sym(score_var)

  dat %>%
    filter(!is.na(!!score_sym)) %>%
    count(level = !!score_sym, name = "n") %>%
    complete(level = 0:max_domains, fill = list(n = 0)) %>%
    mutate(
      composite = label,
      total_n = sum(n),
      percent = 100 * n / total_n
    ) %>%
    rowwise() %>%
    mutate(
      level_label = paste0(level, ifelse(level == 1, " deficiency", " deficiencies")),
      prop_ci = list(stats::binom.test(n, total_n)$conf.int),
      percent_lb = 100 * prop_ci[[1]][1],
      percent_ub = 100 * prop_ci[[1]][2],
      n_with_outcome = sum(!is.na(dat[[score_var]]) & dat[[score_var]] == level & !is.na(dat[[outcome_var]]))
    ) %>%
    ungroup()
}

burden_level_tbl <- bind_rows(
  summarize_burden_levels(h2a_raw, "mm_domain_primary", 4, "Micronutrient domain burden score"),
  summarize_burden_levels(h2a_raw, "mm_domain_primary_relaxed", 4, "Missing-data sensitivity score"),
  summarize_burden_levels(h2a_raw, "mm_noiron", 3, "Non-iron deficiency burden"),
  summarize_burden_levels(h2a_raw, "mm_bdef", 2, "B-vitamin deficiency burden"),
  summarize_burden_levels(h2a_raw, "mm_aanemia", 2, "Vitamin A and iron deficiency anemia burden")
) %>%
  transmute(
    Composite = composite,
    `Deficiency level` = level_label,
    `Level n / total N (%)` = paste0(
      n, "/", total_n, " (", formatC(percent, format = "f", digits = 1), "%)"
    ),
    `Percent 95% CI` = fmt_ci(percent_lb, percent_ub, digits = 1, suffix = "%"),
    `N with telomere outcome` = as.integer(n_with_outcome)
  )

primary_tbl <- readr::read_csv(paths$primary_cat, show_col_types = FALSE) %>%
  transmute(
    Comparison = category,
    Outcome = h2a_labels[Y],
    Estimate = fmt_num(coef, 2),
    `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
    `P value` = fmt_p(pval),
    Significance = sig_label(pval),
    N = as.integer(N)
  )

primary_desc_tbl <- h2a_raw %>%
  filter(!is.na(mm_domain_primary_cat), !is.na(TS_t3_Z)) %>%
  group_by(mm_domain_primary_cat) %>%
  summarise(
    n = n(),
    mean_telomere_z = mean(TS_t3_Z, na.rm = TRUE),
    sd_telomere_z = sd(TS_t3_Z, na.rm = TRUE),
    se_telomere_z = sd_telomere_z / sqrt(n),
    mean_lb = mean_telomere_z - 1.96 * se_telomere_z,
    mean_ub = mean_telomere_z + 1.96 * se_telomere_z,
    .groups = "drop"
  ) %>%
  transmute(
    Category = mm_domain_primary_cat,
    `Category N` = as.integer(n),
    `Mean telomere z-score` = fmt_num(mean_telomere_z, 2),
    `Mean 95% CI` = paste0(fmt_num(mean_lb, 2), ", ", fmt_num(mean_ub, 2))
  )

alt_adj_tbl <- readr::read_csv(paths$adj, show_col_types = FALSE) %>%
  transmute(
    Exposure = h2a_labels[X],
    Outcome = h2a_labels[Y],
    Estimate = fmt_num(coef, 2),
    `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
    `P value` = fmt_p(pval),
    Significance = sig_label(pval),
    N = as.integer(N)
  ) %>%
  filter(Exposure != "Missing-data sensitivity score") %>%
  arrange(match(
    Exposure,
    c(
      "Micronutrient domain burden score",
      "Non-iron deficiency burden",
      "B-vitamin deficiency burden",
      "Vitamin A and iron deficiency anemia burden"
    )
  ))

missingness_tbl <- readr::read_csv(paths$missingness, show_col_types = FALSE) %>%
  transmute(
    Analysis = analysis,
    Outcome = h2a_labels[Y],
    Estimate = fmt_num(coef, 2),
    `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
    `P value` = fmt_p(pval),
    Significance = sig_label(pval),
    N = as.integer(N)
  )

alt_unadj_tbl <- readr::read_csv(paths$unadj, show_col_types = FALSE) %>%
  transmute(
    Exposure = h2a_labels[X],
    Outcome = h2a_labels[Y],
    Estimate = fmt_num(coef, 2),
    `95% CI` = paste0(fmt_num(lb, 2), ", ", fmt_num(ub, 2)),
    `P value` = fmt_p(pval),
    Significance = sig_label(pval),
    N = as.integer(N)
  ) %>%
  filter(Exposure != "Missing-data sensitivity score") %>%
  arrange(match(
    Exposure,
    c(
      "Micronutrient domain burden score",
      "Non-iron deficiency burden",
      "B-vitamin deficiency burden",
      "Vitamin A and iron deficiency anemia burden"
    )
  ))

ft_primary <- make_ft(primary_tbl) %>%
  add_footer_lines(
    values = c(
      "Preliminary H2A analysis comparing children with 1 deficiency or 2+ deficiencies to the reference group with no deficient micronutrient domains."
    )
  )

ft_domain_sizes <- make_ft(domain_size_tbl) %>%
  add_footer_lines(
    values = c(
      "All H2A tables are restricted to children with non-missing telomere length z-score (TS_t3_Z).",
      "Domain-specific Ns show the number of children with non-missing values for each deficiency domain.",
      "Because the analytic sample is already restricted to telomere-complete children, domain-specific Ns represent exposure completeness within the telomere-available cohort."
    )
  )

ft_component_domains <- make_ft(component_domain_tbl) %>%
  add_footer_lines(
    values = c(
      "This table maps each H2A burden score to the deficiency domains used to construct it.",
      "Complete composite N is the number of telomere-complete children with all required component domains observed for that score."
    )
  )

ft_composite_n <- make_ft(composite_n_tbl) %>%
  add_footer_lines(
    values = c(
      "All model Ns are drawn from the telomere-complete analytic cohort.",
      "Composite-model N varies because each burden score uses a different set of component domains and missing-data rules."
    )
  )

ft_burden_levels <- make_ft(burden_level_tbl) %>%
  add_footer_lines(
    values = c(
      "Deficiency-level distributions are based on the telomere-complete analytic cohort.",
      "Deficiency-level counts show how many children had 0, 1, 2, 3, or 4 deficient domains within each H2A composite score.",
      "Percent confidence intervals are exact binomial 95% confidence intervals."
    )
  )

ft_primary_desc <- make_ft(primary_desc_tbl) %>%
  add_footer_lines(
    values = c(
      "Descriptive telomere means are shown with Wald 95% confidence intervals within each primary-burden category."
    )
  )

ft_alt_adj <- make_ft(alt_adj_tbl) %>%
  add_footer_lines(
    values = c(
      "Main H2A analysis showing adjusted estimates for the micronutrient domain burden scores."
    )
  )

ft_missing <- make_ft(missingness_tbl) %>%
  add_footer_lines(
    values = c(
      "Sensitivity analysis comparing the main micronutrient domain burden score with the missing-data sensitivity score."
    )
  )

ft_alt_unadj <- make_ft(alt_unadj_tbl) %>%
  add_footer_lines(
    values = c(
      "Unadjusted H2A estimates shown as a supplementary table."
    )
  )

doc <- read_docx() %>%
  body_add_par("H2A result tables", style = "heading 1") %>%
  body_add_par("Table 1. H2A deficiency domain sample sizes", style = "heading 2") %>%
  body_add_flextable(ft_domain_sizes) %>%
  body_add_par("Table 2. H2A component domains within each burden score", style = "heading 2") %>%
  body_add_flextable(ft_component_domains) %>%
  body_add_par("Table 3. H2A deficiency-level distribution within each burden score", style = "heading 2") %>%
  body_add_flextable(ft_burden_levels) %>%
  body_add_par("Table 4. H2A composite-model sample sizes", style = "heading 2") %>%
  body_add_flextable(ft_composite_n) %>%
  body_add_par("Table 5. Primary categorical burden descriptive summary", style = "heading 2") %>%
  body_add_flextable(ft_primary_desc) %>%
  body_add_par("Table 6. Preliminary categorical burden analysis", style = "heading 2") %>%
  body_add_flextable(ft_primary) %>%
  body_add_par("Table 7. Main analysis: adjusted micronutrient domain burden scores", style = "heading 2") %>%
  body_add_flextable(ft_alt_adj) %>%
  body_add_par("Table 8. Sensitivity analysis: missing-data", style = "heading 2") %>%
  body_add_flextable(ft_missing) %>%
  body_add_par("Supplementary Table. Unadjusted burden score analyses", style = "heading 2") %>%
  body_add_flextable(ft_alt_unadj)

print(doc, target = output_path)

cat("H2A results tables saved to:", output_path, "\n")
