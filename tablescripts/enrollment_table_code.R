############################################################
## Table 1: Enrollment Table (WASH Benefits Bangladesh)
## Micronutrients & Telomere at 28 months (t3)
############################################################

## Install needed packages once (uncomment if needed)
# install.packages(c("tidyverse", "gtsummary", "flextable", "officer"))

library(tidyverse)
library(gtsummary)
library(flextable)
library(officer)

############################################################
## 1. Read in your data
############################################################

# Update this path if needed ??? currently assumes file is on Desktop
data <- readr::read_csv(
  "/Users/sumuccu/Desktop/bangladesh-merged-eed-mn-data.csv"
)

############################################################
## 2. Create derived variables (age in months, female)
############################################################

data_enroll <- data %>%
  mutate(
    # Age in months at timepoint 3 (approx: 30.4 days/month)
    age_months_t3 = ageday_bt3 / 30.4,
    
    # Create a 'female' indicator because sex is coded 1 = male, 0 = female
    female = case_when(
      sex == 0 ~ 1,   # 0 = female
      sex == 1 ~ 0,   # 1 = male
      TRUE ~ NA_real_
    ),
    
    # Make sure telomere is numeric
    TS_t3_Z = as.numeric(TS_t3_Z)
  )

############################################################
## 3. Define analysis cohort:
##    Children with telomere + ???1 micronutrient at t3
############################################################

# Continuous micronutrient variables
mn_vars <- c(
  "logFER_inf",
  "hemo_eeel",
  "hepcidin_eeel",
  "logsTFR_inf",
  "logRBP_inf",
  "logB12",
  "logfolate"
)

data_cohort <- data_enroll %>%
  filter(
    !is.na(TS_t3_Z),                                   # must have telomere
    rowSums(!is.na(select(., all_of(mn_vars)))) >= 1   # at least 1 micronutrient
  )

# Ns at each step
cat("Original N:", nrow(data_enroll), "\n")
cat("N with non-missing telomere (TS_t3_Z):",
    sum(!is.na(data_enroll$TS_t3_Z)), "\n")
cat("Analysis cohort N (TS_t3_Z + ???1 micronutrient):",
    nrow(data_cohort), "\n")

############################################################
## 4. Helper: function to make a gtsummary section
############################################################

make_section <- function(data, vars, labels) {
  data %>%
    select(all_of(vars)) %>%
    tbl_summary(
      type = list(
        all_continuous()  ~ "continuous",
        all_dichotomous() ~ "dichotomous"
      ),
      statistic = list(
        # Single row for continuous: Median (IQR)
        all_continuous()  ~ "{median} ({p25}, {p75})",
        # Single row for binary: n (%)
        all_dichotomous() ~ "{n} ({p}%)"
      ),
      digits   = all_continuous() ~ 2,
      label    = labels,
      missing  = "no"            # no separate 'Missing' rows
    )
}

############################################################
## 5. Prepare variables and labels for each section (Option A)
############################################################

## Child characteristics
child_vars <- c("age_months_t3", "female")
child_labels <- list(
  age_months_t3 ~ "Age at 28 months (months)",
  female        ~ "Female"
)

## Anthropometry (28 months)
anthro_vars <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")
anthro_labels <- list(
  laz_t3 ~ "Length-for-age Z score",
  waz_t3 ~ "Weight-for-age Z score",
  whz_t3 ~ "Weight-for-length Z score",
  hcz_t3 ~ "Head circumference-for-age Z score"
)

## Micronutrients (continuous)
micro_cont_vars <- c(
  "logFER_inf",
  "hemo_eeel",
  "hepcidin_eeel",
  "logsTFR_inf",
  "logRBP_inf",
  "logB12",
  "logfolate"
)

micro_cont_labels <- list(
  logFER_inf    ~ "Ferritin (log-transformed)",
  hemo_eeel     ~ "Hemoglobin (g/dL)",
  hepcidin_eeel ~ "Hepcidin",
  logsTFR_inf   ~ "Soluble transferrin receptor (log-transformed)",
  logRBP_inf    ~ "Retinol-binding protein (log-transformed)",
  logB12        ~ "Vitamin B12 (log-transformed)",
  logfolate     ~ "Folate (log-transformed)"
)

## Micronutrient deficiencies (binary)
micro_def_vars <- c(
  "lowfer_inf",
  "irondefanemia_inf",
  "vitadef_inf",
  "folatedef",
  "folatehigh",
  "B12def"
)

micro_def_labels <- list(
  lowfer_inf        ~ "Low ferritin (<12 µg/L)",
  irondefanemia_inf ~ "Iron-deficiency anemia (Hb <11 g/dL)",
  vitadef_inf       ~ "Vitamin A deficiency (RBP <0.83 µmol/L)",
  folatedef         ~ "Folate deficiency (<10 nmol/L)",
  folatehigh        ~ "High folate",
  B12def            ~ "Vitamin B12 deficiency"
)

## Telomere length
telo_vars <- "TS_t3_Z"
telo_labels <- list(
  TS_t3_Z ~ "Telomere length Z-score at 28 months"
)

############################################################
## 6. Recode binary variables as factor (No/Yes)
############################################################

bin_vars <- c("female", micro_def_vars)

data_cohort <- data_cohort %>%
  mutate(
    across(
      all_of(bin_vars),
      ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))
    )
  )

############################################################
## 7. Build each section table and stack with headers
############################################################

tbl_child      <- make_section(data_cohort, child_vars,      child_labels)
tbl_anthro     <- make_section(data_cohort, anthro_vars,     anthro_labels)
tbl_micro_cont <- make_section(data_cohort, micro_cont_vars, micro_cont_labels)
tbl_micro_def  <- make_section(data_cohort, micro_def_vars,  micro_def_labels)
tbl_telo       <- make_section(data_cohort, telo_vars,       telo_labels)

tbl_enrollment <- tbl_stack(
  list(
    tbl_child,
    tbl_anthro,
    tbl_micro_cont,
    tbl_micro_def,
    tbl_telo
  ),
  group_header = c(
    "Child",
    "Anthropometry (28 months)",
    "Micronutrients (continuous)",
    "Micronutrient deficiencies",
    "Telomere length"
  )
) %>%
  modify_header(
    label  ~ "**Characteristic**",
    stat_0 ~ "**n (%) or Median (IQR)**"
  ) %>%
  bold_labels()

############################################################
## 8. Convert to flextable and style
############################################################

ft_enrollment <- tbl_enrollment %>%
  as_flex_table() %>%
  theme_vanilla() %>%
  autofit() %>%
  align(j = 2, align = "center", part = "all")

############################################################
## 9. Save to a LANDSCAPE Word document on Desktop
############################################################

output_path <- "/Users/sumuccu/MicroTelo/tables/Enrollment_table_t3.docx"

doc <- read_docx()

# Add table title
doc <- body_add_par(
  doc,
  "Table 1. Summary of Study Participants??? Characteristics",
  style = "heading 1"
)

# Switch to landscape for the following section
doc <- body_end_section_landscape(doc)

# Add the table
doc <- body_add_flextable(doc, ft_enrollment)

print(doc, target = output_path)

cat("Enrollment table saved to:", output_path, "\n")
