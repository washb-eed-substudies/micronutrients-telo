############################################################
## Table 3: Descriptive statistics for micronutrients
## and telomere length at 28 months
############################################################

library(tidyverse)
library(flextable)
library(officer)
library(here)

data_path <- "/Users/sumuccu/Desktop/bangladesh-merged-eed-mn-data.csv"
output_path <- here("tables", "Table3_descriptive_statistics_t3.docx")

data <- readr::read_csv(data_path, show_col_types = FALSE)

data_table3 <- data %>%
  mutate(
    TS_t3 = as.numeric(TS_t3),
    TS_t3_Z = as.numeric(TS_t3_Z)
  ) %>%
  filter(!is.na(TS_t3_Z))

cat("Original N:", nrow(data), "\n")
cat("Table 3 cohort N (non-missing telomere z-score):", nrow(data_table3), "\n")

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

section_header <- function(title) {
  tibble(
    Characteristic = title,
    `N non-missing` = NA_integer_,
    `Median (IQR) or n (%)` = ""
  )
}

micronutrient_continuous_section <- bind_rows(
  section_header("Micronutrient concentrations at 28 months"),
  continuous_summary(data_table3, "FER_eeel", "Ferritin, ug/L"),
  continuous_summary(data_table3, "hemo_eeel", "Hemoglobin, g/dL"),
  continuous_summary(data_table3, "hepcidin_eeel", "Hepcidin, ng/mL"),
  continuous_summary(data_table3, "sTFR_eeel", "Soluble transferrin receptor, mg/L"),
  continuous_summary(data_table3, "RBP_eeel", "Retinol-binding protein, umol/L"),
  continuous_summary(data_table3, "B12_eeel", "Vitamin B12, pmol/L"),
  continuous_summary(data_table3, "folate_eeel", "Folate, nmol/L")
)

micronutrient_binary_section <- bind_rows(
  section_header("Micronutrient deficiency indicators at 28 months"),
  binary_summary(data_table3, "lowfer_inf", "Low ferritin (<12 ug/L)"),
  binary_summary(data_table3, "highsTFR", "High soluble transferrin receptor (>8.3 mg/L)"),
  binary_summary(data_table3, "irondef_inf", "Iron deficiency (ferritin <12 ug/L or sTfR >8.3 mg/L)"),
  binary_summary(data_table3, "irondefanemia_inf", "Iron deficiency anemia (hemoglobin <11 g/dL and low ferritin or high sTfR)"),
  binary_summary(data_table3, "vit_A_def", "Vitamin A deficiency (RBP <0.83 umol/L)"),
  binary_summary(data_table3, "B12def", "Vitamin B12 deficiency (<150 pmol/L)"),
  binary_summary(data_table3, "B12marg", "Marginal vitamin B12 (150 to <221 pmol/L)"),
  binary_summary(data_table3, "folatedef", "Folate deficiency (<10 nmol/L)")
)

telomere_section <- bind_rows(
  section_header("Telomere length at 28 months"),
  continuous_summary(data_table3, "TS_t3", "Telomere length"),
  continuous_summary(data_table3, "TS_t3_Z", "Telomere length z-score")
)

tbl_table3 <- bind_rows(
  micronutrient_continuous_section,
  micronutrient_binary_section,
  telomere_section
)

section_rows <- which(is.na(tbl_table3$`N non-missing`))
deficiency_rows <- which(
  grepl("\\(", tbl_table3$Characteristic) &
    !is.na(tbl_table3$`N non-missing`)
)

ft_table3 <- flextable(tbl_table3) %>%
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
  width(j = 1, width = 4.6) %>%
  width(j = 2, width = 1.0) %>%
  width(j = 3, width = 2.1) %>%
  fontsize(size = 8.8, part = "all") %>%
  padding(padding = 2, part = "all") %>%
  line_spacing(space = 0.95, part = "all") %>%
  fit_to_width(max_width = 7.6)

for (i in deficiency_rows) {
  label <- tbl_table3$Characteristic[i]
  base_label <- sub("\\s*\\(.*$", "", label)
  threshold_label <- sub("^.*?(\\(.*\\))$", "\\1", label)

  ft_table3 <- compose(
    ft_table3,
    i = i,
    j = "Characteristic",
    value = as_paragraph(
      as_chunk(base_label),
      " ",
      as_chunk(threshold_label, props = fp_text_default(italic = TRUE))
    ),
    part = "body"
  )
}

doc <- read_docx()
doc <- body_add_par(
  doc,
  "Table 3. Descriptive statistics for micronutrient biomarkers, micronutrient deficiency indicators, and telomere length at 28 months",
  style = "heading 1"
)
doc <- body_add_par(
  doc,
  sprintf(
    "Values are presented as median (interquartile range) or n (%%). The descriptive cohort included children with non-missing telomere length z-score at 28 months (N = %d).",
    nrow(data_table3)
  ),
  style = "Normal"
)
doc <- body_add_flextable(doc, ft_table3)

print(doc, target = output_path)

cat("Table 3 saved to:", output_path, "\n")
