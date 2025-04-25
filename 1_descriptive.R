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

# Data used in this code:
# dth_dt_cln.rds includes all health and exposure data, with one record for each quarter-person.

## Table one (summary statistics for cohort characteristics and health outcomes)
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

## Table one (summary statistics for exposure metrics)
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

## Table S one (summary statistics for cohort and subgroup 
## cumulative death and loss-to-follow-up, follow-up time)
## stratified dataset creation is included in 3_stratified_analysis.R
####################
## Entire cohort
#Bring in data (only for follow-up time, the other two in the chunk "Table one (summary statistics for cohort characteristics and health outcomes)")
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_cln.rds")),
               select=c(studyid,obs_num,censoring,death,yr))

dt_cln<-dt_cln %>%
  mutate(time_period=obs_num-12)

dt_cln_last<- dt_cln %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() 


dt_cln_tbl<- dt_cln_last  %>%
  mutate(lst_2_flwup=case_when(censoring==1 & death==0 ~ 1,
                               .default=0),
         lost_2_flwup=relevel(as.factor(case_when(lst_2_flwup==0 ~ "No",
                                                  lst_2_flwup==1 ~ "Yes")),ref="No")
         )

dt_cln_last %>%
  dplyr::select(time_period,yr) %>%
  tbl_summary(label= list(
    time_period="# of quarters in cohort before lost to follow-up",
    yr="Year lost to follow-up"),
    type= all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "{mean} ({sd})",
      "{median} ({p25}, {p75})",
      "{min}, {max}"))



## Stratified by sex
## female
#Bring in data
mort_dt_f<-subset(readRDS(here(indir1, "data", "sex_stratified_censoring","mort_dt_f.rds")),
                  select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_f<-mort_dt_f %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_f_1st<- mort_dt_f %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr,time_period))


mort_dt_f_last<- mort_dt_f %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr,time_period)

mort_dt_f_lastjn<-mort_dt_f_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_f_dth<- mort_dt_f_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_f_rmn2<-mort_dt_f_1st %>%
  anti_join(mort_dt_f_dth, by = join_by(studyid)) %>%
  dplyr::select(studyid) 
#Get censoring 
mort_dt_f_cnsr<- mort_dt_f_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_f_cnsrrmn2<-mort_dt_f_1st %>%
  anti_join(mort_dt_f_cnsr, by = join_by(studyid)) %>%
  dplyr::select(studyid) 
#Consolidate into one data frame
mort_dt_f_typ<-mort_dt_f_dth %>%
  bind_rows(mort_dt_f_rmn2,.id="death_status")

mort_dt_f_typ2<-mort_dt_f_cnsr %>%
  bind_rows(mort_dt_f_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_f_tbl<-mort_dt_f_1st %>%
  left_join(mort_dt_f_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_f_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_f_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_f_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Females")



## male
mort_dt_m<-subset(readRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_m.rds")),
                  select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_m<-mort_dt_m %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_m_1st<- mort_dt_m %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr,time_period))


mort_dt_m_last<- mort_dt_m %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr,time_period)

mort_dt_m_lastjn<-mort_dt_m_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_m_dth<- mort_dt_m_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_m_rmn2<-mort_dt_m_1st %>%
  anti_join(mort_dt_m_dth, by = join_by(studyid)) %>%
  dplyr::select(studyid) 
#Get censoring 
mort_dt_m_cnsr<- mort_dt_m_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_m_cnsrrmn2<-mort_dt_m_1st %>%
  anti_join(mort_dt_m_cnsr, by = join_by(studyid)) %>%
  dplyr::select(studyid) 
#Consolidate into one data frame
mort_dt_m_typ<-mort_dt_m_dth %>%
  bind_rows(mort_dt_m_rmn2,.id="death_status")

mort_dt_m_typ2<-mort_dt_m_cnsr %>%
  bind_rows(mort_dt_m_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_m_tbl<-mort_dt_m_1st %>%
  left_join(mort_dt_m_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_m_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_m_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_m_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Males")

## stratified by race and ethnicity (same as sex except for the categories)

## stratified by age
## Age 75+
#Bring in data
mort_dt_75p<-subset(readRDS(here(indir1, "data", "age_stratified_censoring", 
                                 "mort_dt_75p_at_qlfy.rds")), ## used age stratified data based on age at cohort entry (codes not included)
                    select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_75p<-mort_dt_75p %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_75p_1st<- mort_dt_75p %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr,time_period))


mort_dt_75p_last<- mort_dt_75p %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr, time_period)

mort_dt_75p_lastjn<-mort_dt_75p_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_75p_dth<- mort_dt_75p_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_75p_rmn2<-mort_dt_75p_1st %>%
  anti_join(mort_dt_75p_dth,by=join_by(studyid)) %>%
  dplyr::select(studyid) 

#Get censoring 
mort_dt_75p_cnsr<- mort_dt_75p_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_75p_cnsrrmn2<-mort_dt_75p_1st %>%
  anti_join(mort_dt_75p_cnsr,by=join_by(studyid)) %>%
  dplyr::select(studyid) 

#Consolidate into one data frame
mort_dt_75p_typ<-mort_dt_75p_dth %>%
  bind_rows(mort_dt_75p_rmn2,.id="death_status")

mort_dt_75p_typ2<-mort_dt_75p_cnsr %>%
  bind_rows(mort_dt_75p_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_75p_tbl<-mort_dt_75p_1st %>%
  left_join(mort_dt_75p_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_75p_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_75p_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_75p_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Age 75+ at Cohort Entry")


## Age <75
mort_dt_u75<-subset(readRDS(here(indir1, "data", "age_stratified_censoring",
                                 "mort_dt_u75_at_qlfy.rds")), ## used age stratified data based on age at cohort entry (codes not included)
                    select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_u75<-mort_dt_u75 %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_u75_1st<- mort_dt_u75 %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr, time_period)) %>%
  arrange(studyid)


mort_dt_u75_last<- mort_dt_u75 %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr, time_period)  %>%
  arrange(studyid)

mort_dt_u75_lastjn<-mort_dt_u75_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_u75_dth<- mort_dt_u75_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_u75_rmn2<-mort_dt_u75_1st %>%
  anti_join(mort_dt_u75_dth,by=join_by(studyid)) %>%
  dplyr::select(studyid) 

#Get censoring 
mort_dt_u75_cnsr<- mort_dt_u75_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_u75_cnsrrmn2<-mort_dt_u75_1st %>%
  anti_join(mort_dt_u75_cnsr,by=join_by(studyid)) %>%
  dplyr::select(studyid) 

#Consolidate into one data frame
mort_dt_u75_typ<-mort_dt_u75_dth %>%
  bind_rows(mort_dt_u75_rmn2,.id="death_status")

mort_dt_u75_typ2<-mort_dt_u75_cnsr %>%
  bind_rows(mort_dt_u75_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_u75_tbl<-mort_dt_u75_1st %>%
  left_join(mort_dt_u75_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_u75_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_u75_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_u75_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Age <75 at Cohort Entry")

## stratified by poverty
## Poverty >=15%
#Bring in data
mort_dt_pov<-subset(readRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_pov.rds")),
                    select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_pov<-mort_dt_pov %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_pov_1st<- mort_dt_pov %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr,time_period))


mort_dt_pov_last<- mort_dt_pov %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr,time_period)

mort_dt_pov_lastjn<-mort_dt_pov_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_pov_dth<- mort_dt_pov_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_pov_rmn2<-mort_dt_pov_1st %>%
  anti_join(mort_dt_pov_dth, by = join_by(studyid)) %>%
  dplyr::select(studyid) 

#Get censoring 
mort_dt_pov_cnsr<- mort_dt_pov_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_pov_cnsrrmn2<-mort_dt_pov_1st %>%
  anti_join(mort_dt_pov_cnsr, by = join_by(studyid)) %>%
  dplyr::select(studyid) 

#Consolidate into one data frame
mort_dt_pov_typ<-mort_dt_pov_dth %>%
  bind_rows(mort_dt_pov_rmn2,.id="death_status")

mort_dt_pov_typ2<-mort_dt_pov_cnsr %>%
  bind_rows(mort_dt_pov_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_pov_tbl<-mort_dt_pov_1st %>%
  left_join(mort_dt_pov_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_pov_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_pov_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_pov_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Poverty >=15%")

## Poverty <15%
mort_dt_ls_pov<-subset(readRDS(here(indir1, "data", "poverty_stratified_censoring", "mort_dt_ls_pov.rds")),
                       select=c(studyid,obs_num,censor,death_nw,year))

mort_dt_ls_pov<-mort_dt_ls_pov %>%
  mutate(time_period=obs_num-12,
         yr=as.numeric(as.character(year)))

#Get first record for each participant
mort_dt_ls_pov_1st<- mort_dt_ls_pov %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-c(year,yr,time_period))


mort_dt_ls_pov_last<- mort_dt_ls_pov %>%
  group_by(studyid) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  dplyr::select(studyid, death_nw, censor,yr,time_period)

mort_dt_ls_pov_lastjn<-mort_dt_ls_pov_last %>%
  dplyr::select(studyid,yr,time_period)

#Get deaths 
mort_dt_ls_pov_dth<- mort_dt_ls_pov_last %>%
  dplyr::select(studyid,death_nw) %>%
  group_by(studyid) %>%
  filter(death_nw==1) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#No death
mort_dt_ls_pov_rmn2<-mort_dt_ls_pov_1st %>%
  anti_join(mort_dt_ls_pov_dth, by = join_by(studyid)) %>%
  dplyr::select(studyid) 

#Get censoring 
mort_dt_ls_pov_cnsr<- mort_dt_ls_pov_last %>%
  dplyr::select(studyid,censor,death_nw) %>%
  group_by(studyid) %>%
  filter(censor==1 & is.na(death_nw)) %>%
  ungroup() %>%
  dplyr::select(studyid) 

#no censoring
mort_dt_ls_pov_cnsrrmn2<-mort_dt_ls_pov_1st %>%
  anti_join(mort_dt_ls_pov_cnsr,by = join_by(studyid)) %>%
  dplyr::select(studyid) 

#Consolidate into one data frame
mort_dt_ls_pov_typ<-mort_dt_ls_pov_dth %>%
  bind_rows(mort_dt_ls_pov_rmn2,.id="death_status")

mort_dt_ls_pov_typ2<-mort_dt_ls_pov_cnsr %>%
  bind_rows(mort_dt_ls_pov_cnsrrmn2,.id="censor_status")


#Put in data frame for table 1
mort_dt_ls_pov_tbl<-mort_dt_ls_pov_1st %>%
  left_join(mort_dt_ls_pov_typ,by=join_by(studyid)) %>%
  left_join(mort_dt_ls_pov_typ2,by=join_by(studyid)) %>%
  left_join(mort_dt_ls_pov_lastjn,by=join_by(studyid))%>%
  mutate(death_status=factor(death_status,labels=c("Death before study end","End of follow-up not from death")),
         lost_2_flwup=factor(censor_status,labels=c("Lost to follow-up excluding death","Death or in cohort to study end")))


mort_dt_ls_pov_tbl %>%
  dplyr::select(death_status,lost_2_flwup,time_period,yr) %>%
  tbl_summary(label= list(death_status="Mortality",
                          lost_2_flwup="Lost to Follow-up",
                          time_period="# of quarters in cohort before lost to follow-up",
                          yr="Year lost to follow-up"),
              type= all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "{mean} ({sd})",
                "{median} ({p25}, {p75})",
                "{min}, {max}")) %>%
  modify_header(label="Death and Lost to Follow-up Summary for Poverty <15%")
####################
