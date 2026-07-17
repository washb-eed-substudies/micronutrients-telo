############################################################
## Table 1: Thesis-ready enrollment table
## Child and anthropometric characteristics only
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

data_path <- "/Users/sumuccu/Desktop/bangladesh-merged-eed-mn-data.csv"
output_path <- here("tables", "Enrollment_table_t3.docx")

data <- readr::read_csv(data_path, show_col_types = FALSE)

data_enroll <- data %>%
  mutate(
    age_months_t3 = ageday_bt3 / 30.4,
    birthord = suppressWarnings(as.numeric(birthord)),
    female = case_when(
      sex == 0 ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    TS_t3_Z = as.numeric(TS_t3_Z)
  )

# Final cohort for the thesis enrollment table:
# children with non-missing telomere outcome at 28 months.
data_cohort <- data_enroll %>%
  filter(!is.na(TS_t3_Z))

cat("Original N:", nrow(data_enroll), "\n")
cat("Analysis cohort N (non-missing telomere):", nrow(data_cohort), "\n")

continuous_summary <- function(data, var, label) {
  x <- data[[var]]
  nonmissing_n <- sum(!is.na(x))

  tibble(
    Characteristic = label,
    `N non-missing` = nonmissing_n,
    `Median (IQR) or n (%)` = if (nonmissing_n == 0) {
      "NA"
    } else {
      sprintf(
        "%.2f (%.2f, %.2f)",
        median(x, na.rm = TRUE),
        quantile(x, 0.25, na.rm = TRUE),
        quantile(x, 0.75, na.rm = TRUE)
      )
    }
  )
}

binary_summary <- function(data, var, label) {
  x <- data[[var]]
  nonmissing_n <- sum(!is.na(x))
  yes_n <- sum(x == 1, na.rm = TRUE)
  pct <- if (nonmissing_n == 0) NA_real_ else 100 * yes_n / nonmissing_n

  tibble(
    Characteristic = label,
    `N non-missing` = nonmissing_n,
    `Median (IQR) or n (%)` = if (nonmissing_n == 0) {
      "NA"
    } else {
      sprintf("%d (%.1f%%)", yes_n, pct)
    }
  )
}

categorical_summary <- function(data, var, label) {
  x <- data[[var]]
  nonmissing_n <- sum(!is.na(x))
  counts <- sort(table(x), decreasing = TRUE)

  tibble(
    Characteristic = label,
    `N non-missing` = nonmissing_n,
    `Median (IQR) or n (%)` = if (nonmissing_n == 0) {
      "NA"
    } else {
      paste(
        paste0(names(counts), ": ", as.integer(counts), " (", round(100 * counts / nonmissing_n, 1), "%)"),
        collapse = "; "
      )
    }
  )
}

section_header <- function(title) {
  tibble(
    Characteristic = title,
    `N non-missing` = NA_integer_,
    `Median (IQR) or n (%)` = ""
  )
}

child_section <- bind_rows(
  section_header("Child characteristics"),
  continuous_summary(data_cohort, "age_months_t3", "Age at 28 months, months"),
  binary_summary(data_cohort, "female", "Female sex"),
  continuous_summary(data_cohort, "birthord", "Birth order"),
  categorical_summary(data_cohort, "tr", "Study arm")
)

anthro_section <- bind_rows(
  section_header("Anthropometry at 28 months"),
  continuous_summary(data_cohort, "laz_t3", "Length-for-age z-score"),
  continuous_summary(data_cohort, "waz_t3", "Weight-for-age z-score"),
  continuous_summary(data_cohort, "whz_t3", "Weight-for-length z-score"),
  continuous_summary(data_cohort, "hcz_t3", "Head circumference-for-age z-score")
)

maternal_household_section <- bind_rows(
  section_header("Maternal and household characteristics"),
  continuous_summary(data_cohort, "momage", "Maternal age, years"),
  continuous_summary(data_cohort, "momheight", "Maternal height, cm"),
  continuous_summary(data_cohort, "momedu", "Maternal education, years"),
  categorical_summary(data_cohort, "hfiacat", "Household food insecurity"),
  continuous_summary(data_cohort, "HHwealth_scaled", "Household wealth score"),
  continuous_summary(data_cohort, "Ncomp", "Household size"),
  continuous_summary(data_cohort, "Nlt18", "Household members younger than 18 years")
)

recent_illness_section <- bind_rows(
  section_header("Recent child illness at 28 months"),
  binary_summary(data_cohort, "diar7d_t3", "Diarrhea in past 7 days"),
  binary_summary(data_cohort, "ari7d_t3", "Acute respiratory illness in past 7 days"),
  binary_summary(data_cohort, "fever7d_t3", "Fever in past 7 days")
)

telomere_section <- bind_rows(
  section_header("Telomere outcome"),
  continuous_summary(data_cohort, "TS_t3_Z", "Telomere length z-score at 28 months")
)

tbl_enrollment <- bind_rows(
  child_section,
  anthro_section,
  maternal_household_section,
  recent_illness_section,
  telomere_section
)

section_rows <- which(is.na(tbl_enrollment$`N non-missing`))

ft_enrollment <- flextable(tbl_enrollment) %>%
  set_header_labels(
    Characteristic = "Characteristic",
    `N non-missing` = "N non-missing",
    `Median (IQR) or n (%)` = "Median (IQR) or n (%)"
  ) %>%
  theme_vanilla() %>%
  bold(part = "header") %>%
  bold(i = section_rows, bold = TRUE, part = "body") %>%
  bg(i = section_rows, bg = "#EAEAEA", part = "body") %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = c(2, 3), align = "center", part = "all") %>%
  valign(valign = "center", part = "all") %>%
  width(j = 1, width = 3.0) %>%
  width(j = 2, width = 1.0) %>%
  width(j = 3, width = 3.3) %>%
  fontsize(size = 9, part = "all") %>%
  padding(padding = 3, part = "all") %>%
  line_spacing(space = 0.95, part = "all") %>%
  fit_to_width(max_width = 7.2)

doc <- read_docx()
doc <- body_add_par(
  doc,
  "Table 1. Enrollment and baseline characteristics of children included in the telomere analysis at 28 months",
  style = "heading 1"
)
doc <- body_add_par(
  doc,
  sprintf(
    "Values are presented as median (interquartile range) or n (%%). The analysis cohort included children with non-missing telomere length z-score at 28 months (N = %d).",
    nrow(data_cohort)
  ),
  style = "Normal"
)
doc <- body_add_flextable(doc, ft_enrollment)

print(doc, target = output_path)

cat("Enrollment table saved to:", output_path, "\n")
