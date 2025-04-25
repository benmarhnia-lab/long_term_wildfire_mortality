####################
# Title: 1_descriptive.R
# Date Created: 4/17/2025
# Author: Tim B. Frankland; compiled by Chen Chen
# Purpose: provide cohort descriptives for the long-term wildfire smoke on mortality KPSC project
####################

library(data.table)
library(here)
library(readr)
library(tidyverse)
library(purrr)
library(dplyr)
library(gtsummary)
library(gt)
library(paletteer)
library(RColorBrewer)
library(ggplot2)
library(splines)
library(corrr)
library(rcompanion)
library(polycor)
library(psych)

indir1 <- "" ## work directory

## Table one (cohort characteristics and health outcomes)
####################
##Bring in data
dt_cln<-readRDS(here(indir1, "data","dth_dt_cln.rds"))

#Get first record for each participant
dt_cln_1st<- dt_cln %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup()


dt_cln_last<- dt_cln %>%
  group_by(studyid) %>%
  slice(n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death, censoring)

#Get patients that pass away in the study period
dth_cln_pat_cnt<- dt_cln_last %>%
  filter(death==1) %>%
  distinct(studyid) %>%
  count()

#Get deaths 
dt_cln_dth<- dt_cln_last %>%
  dplyr::select(studyid,death) %>%
  group_by(studyid) %>%
  filter(death==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
dt_cln_rmn2<-dt_cln_1st %>%
  anti_join(dt_cln_dth) %>%
  dplyr::select(studyid) 


#Consolidate into one data frame
dt_cln_typ<-dt_cln_dth %>%
  bind_rows(dt_cln_rmn2,.id="death_status")


#Put in data frame for table 1
dt_cln_tbl<-dt_cln_1st %>%
  left_join(dt_cln_typ)  %>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lst_2_flwup=case_when(censoring==1 & death==0 ~ 1,
                               .default=0),
         lost_2_flwup=relevel(as.factor(case_when(lst_2_flwup==0 ~ "No",
                                                  lst_2_flwup==1 ~ "Yes")),ref="No"),
         r_e=case_when(r_e!="Other" ~ r_e,
                       r_e=="Other" & r_e_full=="Unknown" ~ "Other-Unknown",
                       r_e=="Other" & r_e_full!="Unknown" ~ "Other-Known"))


#Summarize first record data for each participant in table 1
dt_cln_tbl %>%
  dplyr::select(age_at_qlfy, gender, r_e, mrtl_status, rqrs_intrprtr, smkn_status, pov_p, pop_den,death_status, lost_2_flwup) %>%
  tbl_summary(label= list(age_at_qlfy="Patient Age at Cohort Entry",
                          gender="Gender",
                          r_e="Race/Ethnicity",
                          mrtl_status="Marital Status",
                          rqrs_intrprtr="Requires Interpreter",
                          smkn_status="Smoking Status",
                          pov_p="Census Tract Poverty",
                          pop_den="Population Density",
                          death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}"))
####################

## Table one (exposure metrics)
####################
#Bring in data
dt_cln<-readRDS(here(indir1, "data", "dth_dt_cln.rds"))


dt_cln %>%
  select(mean_daily_peak_week,mean_wf_pm,mean_non_wf_pm,non_zero_days,weeks_gt_5,smoke_waves) %>%
  tbl_summary(
    label=list(
      mean_daily_peak_week="Mean Daily Peak Week",
      mean_wf_pm="Mean WF PM",
      mean_non_wf_pm="Mean Non-WF PM",
      non_zero_days="Non-Zero Days",
      weeks_gt_5="Weeks GT 5",
      smoke_waves="Smoke Waves"
    ),
    type= list(c(mean_non_wf_pm, mean_wf_pm, mean_daily_peak_week, non_zero_days) ~ "continuous",c(weeks_gt_5, smoke_waves) ~ "categorical"),
    statistic = list(c(mean_non_wf_pm, mean_wf_pm, mean_daily_peak_week, non_zero_days) ~ c(
      "{median} ({p25}, {p75})"))
  )

dt_cln %>%
  ggplot(aes(x=mean_non_wf_pm)) +
  geom_histogram()

dt_cln %>%
  ggplot(aes(x=mean_wf_pm)) +
  geom_histogram()

dt_cln %>%
  ggplot(aes(x=mean_daily_peak_week)) +
  geom_histogram()

dt_cln %>%
  ggplot(aes(x=weeks_gt_5)) +
  geom_histogram()

dt_cln %>%
  ggplot(aes(x=non_zero_days)) +
  geom_histogram()

dt_cln %>%
  ggplot(aes(x=smoke_waves)) +
  geom_histogram()

dt_cln %>%
  summarise(enframe(quantile(mean_non_wf_pm, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "mean_non_wf_pm"))

dt_cln %>%
  summarise(enframe(quantile(mean_wf_pm, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "mean_wf_pm"))

dt_cln %>%
  summarise(enframe(quantile(mean_daily_peak_week, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "mean_daily_peak_week"))

dt_cln %>%
  summarise(enframe(quantile(weeks_gt_5, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "weeks_gt_5"))

dt_cln %>%
  summarise(enframe(quantile(non_zero_days, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "non_zero_days"))

dt_cln %>%
  summarise(enframe(quantile(smoke_waves, c(0.0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1.0)), "quantile", "smoke_waves"))
####################

## Table S one (cumulative death and loss-to-follow-up, mean follow-up time)
####################

####################
