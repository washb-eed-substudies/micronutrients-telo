rm(list=ls())

#source(here::here("0-config.R"))
library(tidyverse)
library(flextable)
library(officer)
library(data.table)

d <- readRDS("/Users/farheenjamshed/Downloads/bangladesh-cleaned-master-data (2).RDS")

d <- d %>% filter(pregnancy_telo==1)

d <- d %>% mutate(vit_A_def = ifelse(RBP_inf_preg < 0.7, 1, 0),
                  vit_A_low = ifelse(RBP_inf_preg < 1.05, 1, 0))

filtering <- function(row){
  any(!is.na(row))
}

# MATERNAL PREGNANCY BIOMARKERS UNCOMMENT AND FILL IN THIS CODE (UNCOMMENT WITH CTRL+SHIFT+C ON PC)
exp <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf", "vit_D_def", "vit_A_def", "vit_A_low", "iron_def", "ln_preg_cort", "logCRP", "logAGP", "mom_t0_ln_ifn", "sumscore_t0_mom_Z", "ln_preg_estri") 
out <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z") 
d1 <- d[apply(select(d, all_of(exp)), 1, filtering),] # only has rows where we have exposure data for the mom
d1 <- d1[apply(select(d1, all_of(out)), 1, filtering),] # only has rows where we have both some exposure data and some outcome data (all kids included in analyses)

m <- d1 %>% distinct(dataid, .keep_all = T)

nperc <- function(vector){
  n <- sum(vector==1, na.rm=T)
  perc <- round(n/sum(!is.na(vector))*100)
  paste(n, " (", perc, "%)", sep="")}

mediqr <- function(vector){
  quantiles <- round(quantile(vector, na.rm=T), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

child <- c('agemth_bt2', 'agemth_bt3', 'sex', 'laz_t1','waz_t1','whz_t1','hcz_t1',
           'laz_t2','waz_t2','whz_t2','hcz_t2',
           'laz_t3','waz_t3','whz_t3','hcz_t3','diar7d_t2','diar7d_t3')

mom <- c('momage', 'gest_age_weeks', 'momheight', 'momeduy', 'cesd_sum_t2', 'cesd_sum_ee_t3', 'pss_sum_mom_t3', 'life_viol_any_t3')

hfiacat_ind <- ifelse(m$hfiacat=="Food Secure", 0, 1)
household <- c('hfiacat_ind')
sum_hfiacat_ind <- sum(na.omit(hfiacat_ind))

n_med_col <- NULL
for (var in c(child)) {
  if (var %in% c('sex', 'diar7d_t2', 'diar7d_t3', 'life_viol_any_t3') | is.factor(d[[var]])) {
    if (var == 'sex') {
      n <- sum(d$sex=='female', na.rm=T)
      perc <- round(n/sum(!is.na(d$sex))*100)
      n_med_col <- c(n_med_col, paste(n, " (", perc, "%)", sep=""))
    }else {
      n_med_col <- c(n_med_col, nperc(d[[var]]))
    }
  }else {
    n_med_col <- c(n_med_col, mediqr(d[[var]]))
  }
}

for (var in c(mom)) {
  if (var %in% c('life_viol_any_t3') | is.factor(m[[var]])) {
    n_med_col <- c(n_med_col, nperc(m[[var]]))
  }else {
    n_med_col <- c(n_med_col, mediqr(m[[var]]))
  }
}

n_med_col <- c(n_med_col, nperc(hfiacat_ind))


tbl1 <- data.table("C1" = c("Child", rep("", length(child)-1),
                            "Mother", rep("",length(mom)-1),
                            "Household", rep("",length(household)-1)),
                   "C2" = c("", "", "",
                            "Anthropometry (3 months)","","","",
                            "Anthropometry (14 months)","","","",
                            "Anthropometry (28 months)","","","", 
                            "Diarrhea (14 months)", "Diarrhea (28 months)",
                            "","",
                            "Anthropometry at enrollment", "Education",
                            "Depression (14 months)", "Depression (28 months)",
                            "Perceived stress (28 months)", 
                            "Intimate partner violence",
                            "Household Food Insecurity"),
                   "C3" = c("Age at Year 1 (months)",
                            "Age at Year 2 (months)",
                            "Female",
                            "Length-for-age Z score", 
                            "Weight-for-age Z score", "Weight-for-length Z score", 
                            "Head circumference-for-age Z score",
                            "Length-for-age Z score", 
                            "Weight-for-age Z score", "Weight-for-length Z score", 
                            "Head circumference-for-age Z score",
                            "Length-for-age Z score", 
                            "Weight-for-age Z score", "Weight-for-length Z score", 
                            "Head circumference-for-age Z score",
                            "Caregiver-reported 7-day recall", 
                            "Caregiver-reported 7-day recall", 
                            "Age (years)", "Gestational age (weeks)", "Height (cm)", 
                            "Schooling completed (years)",
                            "CESD-20* score", "CESD-20* score", "Perceived Stress Scale score", "Any lifetime exposure", "Food-insecure households"),
                   "C4" = n_med_col)

tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
tbl1flex <- set_header_labels(tbl1flex,
                              values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or Median (IQR)"))
tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
tbl1flex <- autofit(tbl1flex, part = "all")
tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
tbl1flex <- fit_to_width(tbl1flex, max_width=8)
tbl1flex %>% add_footer_row(top=F, values = "*CESD-20 = Center for Epidemiologic Studies Depression Scale Revised.", colwidths = 4)


#sum(m$hfiacat_ind)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)
save_as_docx("Table 1" = tbl1flex, path="~/Desktop/pregnancy-telo/tables/pregnancy-telo-enrollment-052324.docx", 
             pr_section = sect_properties) 

#table(d$momedu)