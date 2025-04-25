####################
# Title: 3_stratified_analysis.R
# Date Created: 4/17/2025
# Author: Tim B. Frankland; compiled by Chen Chen
# Purpose: conduct stratified analysis for the long-term wildfire smoke 
# on mortality KPSC project. To note, we only showed analytical codes for the metric average 
# wildfire PM2.5 (mean_wf_pm), with differences for other metrics marked by "To note".
####################

library(data.table)
library(splines)
library(meta)
library(here)
library(readr)
library(tidyverse)
library(purrr)
library(broom)
library(gt)
library(gtsummary)
library(epiDisplay)
library(psych)
library(speedglm)

indir1 <- "" ## work directory

####################
## Sex-tratified analysis
####################

## Create sex-stratified dataset: To note, this is done for all metrics
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death, mean_non_wf_pm,mean_wf_pm, mean_daily_peak_week, non_zero_days, weeks_gt_5, smoke_waves, gender, censoring, age_boq_s  ,mrtl_status,r_e, rqrs_intrprtr,smkn_status,pov_p,pop_den))


#IQR
IQR(dt_cln$mean_non_wf_pm)
IQR(dt_cln$mean_wf_pm)
IQR(dt_cln$mean_daily_peak_week)
IQR(dt_cln$non_zero_days)

dt_cln<-dt_cln %>%
  mutate(mean_non_wf_pm_iqr=IQR(mean_non_wf_pm), 
         mean_non_wf_pm=mean_non_wf_pm/mean_non_wf_pm_iqr,
         mean_wf_pm_iqr=IQR(mean_wf_pm),
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         mean_daily_peak_week_iqr=IQR(mean_daily_peak_week),
         mean_daily_peak_week=mean_daily_peak_week/mean_daily_peak_week_iqr,
         non_zero_days_iqr=IQR(non_zero_days),
         non_zero_days=non_zero_days/non_zero_days_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_non_wf_pm_iqr,mean_wf_pm_iqr,mean_daily_peak_week_iqr,non_zero_days_iqr))

#Get values of these variables at cohort entry
dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den,age_boq_s ) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year,
         age_boq_s_0=age_boq_s) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0,age_boq_s_0)

#Join cohort entry variables and then recode censoring and death
dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

rm(dt_cln0)

#Compare original death and censor variables to new
dt_cln %>%
  dplyr::select(studyid,obs_num,death,censoring) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death,censoring)

dt_cln %>%
  dplyr::select(studyid,obs_num,death_nw,censor) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death_nw,censor)

#Create study time period variable to place into censoring weights function
dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

#Split the dataset into two gender cohorts
mort_dt_f<-dt_cln %>%
  filter(gender=="F")

mort_dt_m<-dt_cln %>%
  filter(gender=="M")

rm(dt_cln)

mort_dt_f %>%
  saveRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_f.rds"))

mort_dt_m %>%
  saveRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_m.rds"))
####################

## create sex-stratified IPCW
####################
## below are variable names that should be manually entered
tf_cv <- c("mrtl_status","r_e","rqrs_intrprtr","smkn_status")
tv_cv<- c("pov_p","pop_den","year","age_boq_s") 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
time_variable <- "study_time_period"

## for female
#Bring in data
mort_dt_f<-readRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_f.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_f, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_f, family = binomial(), sparse=FALSE)

mort_dt_f$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_f)

mort_dt_f$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_f)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_f$w_c_i_mwfpm <- 1/mort_dt_f$ps.censoring1_mwfpm
mort_dt_f<- mort_dt_f %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_f$w_c_mwfpm)


## stablized weights
mort_dt_f$sw_c_i_mwfpm <- mort_dt_f$ps.censoring2_mwfpm/mort_dt_f$ps.censoring1_mwfpm
mort_dt_f<- mort_dt_f %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_f$sw_c_mwfpm)

mort_dt_f %>%
  saveRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_f.rds"))


## For male
#Bring in data
rm(mort_dt_m)

mort_dt_f<-readRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_m.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_m, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_m, family = binomial(), sparse=FALSE)

mort_dt_m$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_m)

mort_dt_m$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_m)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_m$w_c_i_mwfpm <- 1/mort_dt_m$ps.censoring1_mwfpm
mort_dt_m<- mort_dt_m %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                       .default=w_c_mwfpm))
summary(mort_dt_m$w_c_mwfpm)


## stablized weights
mort_dt_m$sw_c_i_mwfpm <- mort_dt_m$ps.censoring2_mwfpm/mort_dt_m$ps.censoring1_mwfpm
mort_dt_m<- mort_dt_m %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                       .default=sw_c_mwfpm))
summary(mort_dt_m$sw_c_mwfpm)

mort_dt_f %>%
  saveRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_m.rds"))
####################

## run pooled logistic regression
####################
## below are variable names that should be manually entered
exposure_<-c("mean_wf_pm")
tf_cv <- c("mrtl_status","r_e","rqrs_intrprtr","smkn_status")
tv_cv<- c("pov_p","pop_den","year","age_boq_s") 
time_variable <- "study_time_period"

## female
#Bring in data
mort_dt_f<-subset(readRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_f.rds")),
                  select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                           censor, age_boq_s ,mrtl_status,r_e, rqrs_intrprtr,
                           smkn_status,pov_p,pop_den,sw_c_mwfpm))


f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
f_adj <- glm(f2_adj, data = mort_dt_f, family = "binomial", weights=sw_c_mwfpm)
summary(f_adj)
f_adj %>%
  glance()
f_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## male
#Bring in data
mort_dt_m<-subset(readRDS(here(indir1, "data", "sex_stratified_censoring", "mort_dt_m.rds")),
                  select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                           censor, age_boq_s ,mrtl_status,r_e, rqrs_intrprtr,
                           smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
m_adj <- glm(f2_adj, data = mort_dt_m, family = "binomial", weights=sw_c_mwfpm)
summary(m_adj)
m_adj %>%
  glance()
m_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## Heterogeneity test
f_adj_wf<-f_adj %>%
  tidy() %>%
  filter(term==exposure_)

m_adj_wf<-m_adj %>%
  tidy() %>%
  filter(term==exposure_)


adj_wf<-f_adj_wf %>%
  bind_rows(m_adj_wf)

poo <- sum(adj_wf$estimate/adj_wf$std.error^2)/(sum(1/adj_wf$std.error^2))
Q <- sum((adj_wf$estimate-poo)^2/adj_wf$std.error^2)
df <- nrow(adj_wf) - 1
pval.Q <- pchisq(Q, df=df, lower.tail = FALSE) 

cat("the Cochran's Q test statistics for the null hypothesis of no heterogeneity in effect estimates across strata", 
    round(Q, digits = 2), 
    "with p-value (should be compared to type 1 error of 0.1 due to low statistical power of the heterogeneity test) of",
    pval.Q, "\n")
####################



####################
## Age-stratified analysis
####################

## Create age-stratified dataset: To note, this is done for all metrics
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death, mean_non_wf_pm,mean_wf_pm, 
                        mean_daily_peak_week, non_zero_days, weeks_gt_5, smoke_waves, 
                        gender, censoring, age_beg_of_qtr ,mrtl_status,r_e, rqrs_intrprtr,
                        smkn_status,pov_p,pop_den))


#IQR
IQR(dt_cln$mean_non_wf_pm)
IQR(dt_cln$mean_wf_pm)
IQR(dt_cln$mean_daily_peak_week)
IQR(dt_cln$non_zero_days)

dt_cln<-dt_cln %>%
  mutate(mean_non_wf_pm_iqr=IQR(mean_non_wf_pm),
         mean_non_wf_pm=mean_non_wf_pm/mean_non_wf_pm_iqr,
         mean_wf_pm_iqr=IQR(mean_wf_pm),
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         mean_daily_peak_week_iqr=IQR(mean_daily_peak_week),
         mean_daily_peak_week=mean_daily_peak_week/mean_daily_peak_week_iqr,
         non_zero_days_iqr=IQR(non_zero_days),
         non_zero_days=non_zero_days/non_zero_days_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_non_wf_pm_iqr,mean_wf_pm_iqr,mean_daily_peak_week_iqr,non_zero_days_iqr))

#Get values of these variables at cohort entry
dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0)

#Join cohort entry variables and then recode censoring and death
dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

rm(dt_cln0)

#Compare original death and censor variables to new
dt_cln %>%
  dplyr::select(studyid,obs_num,death,censoring) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death,censoring)

dt_cln %>%
  dplyr::select(studyid,obs_num,death_nw,censor) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death_nw,censor)

#Create study time period variable to place into censoring weights function
dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

#Split the dataset into two age cohorts
mort_dt_75p<-dt_cln %>%
  filter(age_beg_of_qtr>=75)

mort_dt_u75<-dt_cln %>%
  filter(age_beg_of_qtr<75)

rm(dt_cln)
mort_dt_75p %>%
  saveRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_75p.rds"))

mort_dt_u75 %>%
  saveRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_u75.rds"))
####################

## create age-stratified IPCW
####################
## below are variable names that should be manually entered
tf_cv <- c("gender","mrtl_status","r_e","rqrs_intrprtr","smkn_status")
tv_cv<- c("pov_p","pop_den","year") 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
time_variable <- "study_time_period"

## Censoring Weight Calculations Age 75+ Mean Wildfire PM2.5
#Bring in data
mort_dt_75p<-readRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_75p.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_75p, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_75p, family = binomial(), sparse=FALSE)

mort_dt_75p$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_75p)

mort_dt_75p$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_75p)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_75p$w_c_i_mwfpm <- 1/mort_dt_75p$ps.censoring1_mwfpm
mort_dt_75p<- mort_dt_75p %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_75p$w_c_mwfpm)


## stablized weights
mort_dt_75p$sw_c_i_mwfpm <- mort_dt_75p$ps.censoring2_mwfpm/mort_dt_75p$ps.censoring1_mwfpm
mort_dt_75p<- mort_dt_75p %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_75p$sw_c_mwfpm)

mort_dt_75p %>%
  saveRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_75p.rds"))


## Censoring Weight Calculations Age <75 Mean Wildfire PM2.5
#Bring in data
mort_dt_u75<-readRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_u75.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"
f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_u75, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_u75, family = binomial(), sparse=FALSE)

mort_dt_u75$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_u75)

mort_dt_u75$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_u75)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_u75$w_c_i_mwfpm <- 1/mort_dt_u75$ps.censoring1_mwfpm
mort_dt_u75<- mort_dt_u75 %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_u75$w_c_mwfpm)

## stablized weights
mort_dt_u75$sw_c_i_mwfpm <- mort_dt_u75$ps.censoring2_mwfpm/mort_dt_u75$ps.censoring1_mwfpm
mort_dt_u75<- mort_dt_u75 %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_u75$sw_c_mwfpm)

mort_dt_u75 %>%
  saveRDS(here(indir1, "data", "age_stratified_censoring", "mort_dt_u75.rds"))
####################

## run pooled logistic regression
####################
## below are variable names that should be manually entered
tf_cv <- c("gender", "mrtl_status","r_e","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv<- c("pov_p","pop_den","year")  ## variable names for time-varying covariates 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"

## Age 75+
#Bring in data
mort_dt_75p<-subset(readRDS(here(indir1, "data", "age_stratified_censoring","mort_dt_75p.rds")),
                    select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                             censor, age_beg_of_qtr ,mrtl_status,r_e, rqrs_intrprtr,
                             smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
ap_adj <- glm(f2_adj, data = mort_dt_75p, family = "binomial", weights=sw_c_mwfpm)
summary(ap_adj)
ap_adj %>%
  glance()
ap_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## Age <75
#Bring in data
mort_dt_u75<-subset(readRDS(here(indir1, "data", "age_stratified_censoring","mort_dt_u75.rds")),
                    select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                             censor, age_beg_of_qtr ,mrtl_status,r_e, rqrs_intrprtr,
                             smkn_status,pov_p,pop_den,sw_c_mwfpm))
f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
au_adj <- glm(f2_adj, data = mort_dt_u75, family = "binomial", weights=sw_c_mwfpm)
summary(au_adj)
au_adj %>%
  glance()
au_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## Heterogeneity Test 
ap_adj_wf<-ap_adj %>%
  tidy() %>%
  filter(term==exposure_)

au_adj_wf<-au_adj %>%
  tidy() %>%
  filter(term==exposure_)


adj_wf<-ap_adj_wf %>%
  bind_rows(au_adj_wf)

poo <- sum(adj_wf$estimate/adj_wf$std.error^2)/(sum(1/adj_wf$std.error^2))
Q <- sum((adj_wf$estimate-poo)^2/adj_wf$std.error^2)
df <- nrow(adj_wf) - 1
pval.Q <- pchisq(Q, df=df, lower.tail = FALSE) 

cat("the Cochran's Q test statistics for the null hypothesis of no heterogeneity in effect estimates across strata", 
    round(Q, digits = 2), 
    "with p-value (should be compared to type 1 error of 0.1 due to low statistical power of the heterogeneity test) of",
    pval.Q, "\n")
####################



####################
## Race-stratified analysis
####################

## Create Race-stratified dataset: To note, this is done for all metrics
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data","dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death, mean_non_wf_pm,mean_wf_pm, 
                        mean_daily_peak_week, non_zero_days, weeks_gt_5, smoke_waves, 
                        gender, censoring, age_boq_s  ,mrtl_status,r_e, r_e_full, 
                        rqrs_intrprtr,smkn_status,pov_p,pop_den))


#IQR
IQR(dt_cln$mean_non_wf_pm)
IQR(dt_cln$mean_wf_pm)
IQR(dt_cln$mean_daily_peak_week)
IQR(dt_cln$non_zero_days)

dt_cln<-dt_cln %>%
  mutate(mean_non_wf_pm_iqr=IQR(mean_non_wf_pm),
         mean_non_wf_pm=mean_non_wf_pm/mean_non_wf_pm_iqr,
         mean_wf_pm_iqr=IQR(mean_wf_pm),
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         mean_daily_peak_week_iqr=IQR(mean_daily_peak_week),
         mean_daily_peak_week=mean_daily_peak_week/mean_daily_peak_week_iqr,
         non_zero_days_iqr=IQR(non_zero_days),
         non_zero_days=non_zero_days/non_zero_days_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_non_wf_pm_iqr,mean_wf_pm_iqr,mean_daily_peak_week_iqr,non_zero_days_iqr))

#Get values of these variables at cohort entry
dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den,age_boq_s) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year,
         age_boq_s_0=age_boq_s) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0,age_boq_s_0)

#Join cohort entry variables and then recode censoring and death
dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

rm(dt_cln0)

#Compare original death and censor variables to new
dt_cln %>%
  dplyr::select(studyid,obs_num,death,censoring) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death,censoring)

dt_cln %>%
  dplyr::select(studyid,obs_num,death_nw,censor) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death_nw,censor)

#Create study time period variable to place into censoring weights function
dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

#Split the dataset into race cohorts
mort_dt_asian<-dt_cln %>%
  filter(r_e=="Asian")

mort_dt_black<-dt_cln %>%
  filter(r_e=="Black")

mort_dt_hispanic<-dt_cln %>%
  filter(r_e=="Hispanic")

mort_dt_other_unkwn<-dt_cln %>%
  filter(r_e=="Other" & r_e_full=="Unknown")

mort_dt_other_kwn<-dt_cln %>%
  filter(r_e=="Other" & r_e_full!="Unknown")

mort_dt_white<-dt_cln %>%
  filter(r_e=="White")

rm(dt_cln)
mort_dt_asian %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_asian.rds"))

mort_dt_black %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_black.rds"))

mort_dt_hispanic %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_hispanic.rds"))

mort_dt_other %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_other.rds"))

mort_dt_other_unkwn %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_other_unkwn.rds"))

mort_dt_other_kwn %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_other_kwn.rds"))

mort_dt_white %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_white.rds"))
####################

## create Race-stratified IPCW
####################
## below are variable names that should be manually entered
tf_cv <- c("gender", "mrtl_status","rqrs_intrprtr","smkn_status")
tv_cv<- c("pov_p","pop_den","year","age_boq_s") 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
exposure_ <- "mean_non_wf_pm"
time_variable <- "study_time_period"

## Censoring Weight Calculations Asian Mean wildfire PM2.5
#Bring in data
mort_dt_asian<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_asian.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_asian, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_asian, family = binomial(), sparse=FALSE)

mort_dt_asian$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_asian)

mort_dt_asian$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_asian)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_asian$w_c_i_mwfpm <- 1/mort_dt_asian$ps.censoring1_mwfpm
mort_dt_asian<- mort_dt_asian %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_asian$w_c_mwfpm)


## stablized weights
mort_dt_asian$sw_c_i_mwfpm <- mort_dt_asian$ps.censoring2_mwfpm/mort_dt_asian$ps.censoring1_mwfpm
mort_dt_asian<- mort_dt_asian %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_asian$sw_c_mwfpm)

mort_dt_asian %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_asian.rds"))


## Censoring Weight Calculations Black mean wildfire PM2.5
mort_dt_black<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_black.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_black, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_black, family = binomial(), sparse=FALSE)

mort_dt_black$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_black)

mort_dt_black$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_black)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_black$w_c_i_mwfpm <- 1/mort_dt_black$ps.censoring1_mwfpm
mort_dt_black<- mort_dt_black %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_black$w_c_mwfpm)


## stablized weights
mort_dt_black$sw_c_i_mwfpm <- mort_dt_black$ps.censoring2_mwfpm/mort_dt_black$ps.censoring1_mwfpm
mort_dt_black<- mort_dt_black %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_black$sw_c_mwfpm)

mort_dt_black %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring","mort_dt_black.rds"))



## Censoring Weight Calculations Hispanic Mean wildfire PM2.5
mort_dt_hispanic<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_hispanic.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_hispanic, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_hispanic, family = binomial(), sparse=FALSE)

mort_dt_hispanic$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_hispanic)

mort_dt_hispanic$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_hispanic)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_hispanic$w_c_i_mwfpm <- 1/mort_dt_hispanic$ps.censoring1_mwfpm
mort_dt_hispanic<- mort_dt_hispanic %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_hispanic$w_c_mwfpm)


## stablized weights
mort_dt_hispanic$sw_c_i_mwfpm <- mort_dt_hispanic$ps.censoring2_mwfpm/mort_dt_hispanic$ps.censoring1_mwfpm
mort_dt_hispanic<- mort_dt_hispanic %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_hispanic$sw_c_mwfpm)

mort_dt_hispanic %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_hispanic.rds"))

## Censoring Weight Calculations Other Mean wildfire PM2.5
mort_dt_other<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_other, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_other, family = binomial(), sparse=FALSE)

mort_dt_other$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_other)

mort_dt_other$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_other)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_other$w_c_i_mwfpm <- 1/mort_dt_other$ps.censoring1_mwfpm
mort_dt_other<- mort_dt_other %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_other$w_c_mwfpm)


## stablized weights
mort_dt_other$sw_c_i_mwfpm <- mort_dt_other$ps.censoring2_mwfpm/mort_dt_other$ps.censoring1_mwfpm
mort_dt_other<- mort_dt_other %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_other$sw_c_mwfpm)

mort_dt_other %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other.rds"))


## Censoring Weight Calculations White Mean wildfire PM2.5
mort_dt_white<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_white.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_white, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_white, family = binomial(), sparse=FALSE)

mort_dt_white$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_white)

mort_dt_white$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_white)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_white$w_c_i_mwfpm <- 1/mort_dt_white$ps.censoring1_mwfpm
mort_dt_white<- mort_dt_white %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_white$w_c_mwfpm)


## stablized weights
mort_dt_white$sw_c_i_mwfpm <- mort_dt_white$ps.censoring2_mwfpm/mort_dt_white$ps.censoring1_mwfpm
mort_dt_white<- mort_dt_white %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_white$sw_c_mwfpm)

mort_dt_white %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_white.rds"))


## Censoring Weight Calculations Other-Known Mean wildfire PM2.5
mort_dt_other_kwn<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_kwn.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_other_kwn, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_other_kwn, family = binomial(), sparse=FALSE)

mort_dt_other_kwn$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_other_kwn)

mort_dt_other_kwn$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_other_kwn)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_other_kwn$w_c_i_mwfpm <- 1/mort_dt_other_kwn$ps.censoring1_mwfpm
mort_dt_other_kwn<- mort_dt_other_kwn %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_other_kwn$w_c_mwfpm)


## stablized weights
mort_dt_other_kwn$sw_c_i_mwfpm <- mort_dt_other_kwn$ps.censoring2_mwfpm/mort_dt_other_kwn$ps.censoring1_mwfpm
mort_dt_other_kwn<- mort_dt_other_kwn %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_other_kwn$sw_c_mwfpm)

mort_dt_other_kwn %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_kwn.rds"))


## Censoring Weight Calculations Other-Unknown Mean wildfire PM2.5
mort_dt_other_unkwn<-readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_unkwn.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_other_unkwn, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_other_unkwn, family = binomial(), sparse=FALSE)

mort_dt_other_unkwn$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_other_unkwn)

mort_dt_other_unkwn$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_other_unkwn)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_other_unkwn$w_c_i_mwfpm <- 1/mort_dt_other_unkwn$ps.censoring1_mwfpm
mort_dt_other_unkwn<- mort_dt_other_unkwn %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_other_unkwn$w_c_mwfpm)


## stablized weights
mort_dt_other_unkwn$sw_c_i_mwfpm <- mort_dt_other_unkwn$ps.censoring2_mwfpm/mort_dt_other_unkwn$ps.censoring1_mwfpm
mort_dt_other_unkwn<- mort_dt_other_unkwn %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_other_unkwn$sw_c_mwfpm)

mort_dt_other_unkwn %>%
  saveRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_unkwn.rds"))


####################

## run pooled logistic regression
####################
## below are variable names that should be manually entered
tf_cv <- c("gender", "mrtl_status","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv <- c("pov_p","pop_den","year","age_boq_s") ## variable names for time-varying covariates 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"

## Asian 
#Bring in data
mort_dt_asian<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_asian.rds")),
                      select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                               gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                               smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
a_adj <- glm(f2_adj, data = mort_dt_asian, family = "binomial", weights = sw_c_mwfpm)
summary(a_adj)
a_adj %>%
  glance()
a_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## Black
#Bring in data
mort_dt_black<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_black.rds")),
                      select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                               gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                               smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
b_adj <- glm(f2_adj, data = mort_dt_black, family = "binomial", weights = sw_c_mwfpm)
summary(b_adj)
b_adj %>%
  glance()
b_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## Hispanic 
#Bring in data
mort_dt_hispanic<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_hispanic.rds")),
                      select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                               gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                               smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
h_adj <- glm(f2_adj, data = mort_dt_hispanic, family = "binomial", weights = sw_c_mwfpm)
summary(h_adj)
h_adj %>%
  glance()
h_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## Other 
mort_dt_other<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other.rds")),
                          select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                                   gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                                   smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
o_adj <- glm(f2_adj, data = mort_dt_other, family = "binomial", weights = sw_c_mwfpm)
summary(o_adj)
o_adj %>%
  glance()
o_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## White 
#Bring in data
mort_dt_white<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_white.rds")),
                         select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                                  gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                                  smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
w_adj <- glm(f2_adj, data = mort_dt_white, family = "binomial", weights = sw_c_mwfpm)
summary(w_adj)
w_adj %>%
  glance()
w_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## Other Known 
mort_dt_other_kwn<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_kwn.rds")),
                      select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                               gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                               smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
o_k_adj <- glm(f2_adj, data = mort_dt_other_kwn, family = "binomial", weights = sw_c_mwfpm)
summary(o_k_adj)
o_k_adj %>%
  glance()
o_k_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## Other Unknown
mort_dt_other_unkwn<-subset(readRDS(here(indir1, "data", "race_stratified_censoring", "mort_dt_other_unkwn.rds")),
                      select=c(studyid,obs_num,year, death_nw,  mean_wf_pm, censor, 
                               gender,age_boq_s ,mrtl_status,r_e, r_e_full,rqrs_intrprtr,
                               smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
o_u_adj <- glm(f2_adj, data = mort_dt_other_unkwn, family = "binomial", weights = sw_c_mwfpm)
summary(o_u_adj)
o_u_adj %>%
  glance()
o_u_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))


## Heterogeneity Test 
a_adj_wf<-a_adj %>%
  tidy() %>%
  filter(term==exposure_)

b_adj_wf<-b_adj %>%
  tidy() %>%
  filter(term==exposure_)

h_adj_wf<-h_adj %>%
  tidy() %>%
  filter(term==exposure_)

o_adj_wf<-o_adj %>%
  tidy() %>%
  filter(term==exposure_)

w_adj_wf<-w_adj %>%
  tidy() %>%
  filter(term==exposure_)

adj_wf<-a_adj_wf %>%
  bind_rows(b_adj_wf,h_adj_wf,o_adj_wf,w_adj_wf)

poo <- sum(adj_wf$estimate/adj_wf$std.error^2)/(sum(1/adj_wf$std.error^2))
Q <- sum((adj_wf$estimate-poo)^2/adj_wf$std.error^2)
df <- nrow(adj_wf) - 1
pval.Q <- pchisq(Q, df=df, lower.tail = FALSE) 

cat("the Cochran's Q test statistics for the null hypothesis of no heterogeneity in effect estimates across strata", 
    round(Q, digits = 2), 
    "with p-value (should be compared to type 1 error of 0.1 due to low statistical power of the heterogeneity test) of",
    pval.Q, "\n")
####################



####################
## Poverty-stratified analysis
####################

## Create Poverty-stratified dataset: To note, this is done for all metrics
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data","dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death, mean_non_wf_pm,mean_wf_pm, 
                        mean_daily_peak_week, non_zero_days, weeks_gt_5, smoke_waves, 
                        gender, censoring, age_boq_s  ,mrtl_status,r_e, rqrs_intrprtr,
                        smkn_status,pov_p,pop_den))


#IQR
IQR(dt_cln$mean_non_wf_pm)
IQR(dt_cln$mean_wf_pm)
IQR(dt_cln$mean_daily_peak_week)
IQR(dt_cln$non_zero_days)

dt_cln<-dt_cln %>%
  mutate(mean_non_wf_pm_iqr=IQR(mean_non_wf_pm),
         mean_non_wf_pm=mean_non_wf_pm/mean_non_wf_pm_iqr,
         mean_wf_pm_iqr=IQR(mean_wf_pm),
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         mean_daily_peak_week_iqr=IQR(mean_daily_peak_week),
         mean_daily_peak_week=mean_daily_peak_week/mean_daily_peak_week_iqr,
         non_zero_days_iqr=IQR(non_zero_days),
         non_zero_days=non_zero_days/non_zero_days_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_non_wf_pm_iqr,mean_wf_pm_iqr,mean_daily_peak_week_iqr,non_zero_days_iqr))

#Get values of these variables at cohort entry
dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den,age_boq_s) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year,
         age_boq_s_0=age_boq_s) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0,age_boq_s_0)

#Join cohort entry variables and then recode censoring and death
dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

rm(dt_cln0)

#Compare original death and censor variables to new
dt_cln %>%
  dplyr::select(studyid,obs_num,death,censoring) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death,censoring)

dt_cln %>%
  dplyr::select(studyid,obs_num,death_nw,censor) %>%
  group_by(studyid) %>%
  filter(obs_num==max(obs_num)) %>%
  ungroup() %>%
  count(death_nw,censor)

#Create study time period variable to place into censoring weights function
dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

#Create mean poverty over time 
mean_pov<- dt_cln %>%
  dplyr::select(studyid, death_nw, pov_p) %>%
  filter(death_nw==0) %>%
  group_by(studyid) %>%
  mutate(mean_pov=mean(pov_p)) %>%
  ungroup() %>%
  dplyr::select(studyid,mean_pov) %>%
  distinct(studyid,mean_pov)

mean_pov_dth_1st<- dt_cln %>%
  filter(obs_num==13 & (death_nw==1 | is.na(death_nw))) %>%
  dplyr::select(studyid,pov_p) %>%
  rename(mean_pov=pov_p)

total_mean_pov<- mean_pov %>%
  bind_rows(mean_pov_dth_1st)

dt_cln<- dt_cln %>%
  inner_join(total_mean_pov, by="studyid")

#Split the dataset into two poverty cohorts
mort_dt_pov<-dt_cln %>%
  filter(mean_pov>=15)

mort_dt_ls_pov<-dt_cln %>%
  filter(mean_pov<15)

rm(dt_cln)
mort_dt_pov %>%
  saveRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_pov.rds"))

mort_dt_ls_pov %>%
  saveRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_ls_pov.rds"))
####################

## create Poverty-stratified IPCW
####################
## below are variable names that should be manually entered
tf_cv <- c("gender","mrtl_status","r_e","rqrs_intrprtr","smkn_status")
tv_cv<- c("pop_den","year", "age_boq_s")
tv_cv_0<- paste(tv_cv, "0", sep="_") 
time_variable <- "study_time_period"

## Censoring Weight Calculations for 15% or more poverty Mean Wildfire PM2.5
mort_dt_pov<-readRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_pov.rds"))

#Set exposure
exposure_ <- "mean_wf_pm"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_pov, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_pov, family = binomial(), sparse=FALSE)

mort_dt_pov$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_pov)

mort_dt_pov$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_pov)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_pov$w_c_i_mwfpm <- 1/mort_dt_pov$ps.censoring1_mwfpm
mort_dt_pov<- mort_dt_pov %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_pov$w_c_mwfpm)


## stablized weights
mort_dt_pov$sw_c_i_mwfpm <- mort_dt_pov$ps.censoring2_mwfpm/mort_dt_pov$ps.censoring1_mwfpm
mort_dt_pov<- mort_dt_pov %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_pov$sw_c_mwfpm)

mort_dt_pov %>%
  saveRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_pov.rds"))

## Censoring Weight Calculations for less than 15% poverty Mean Wildfire PM2.5
#Bring in data
mort_dt_ls_pov<-readRDS(here(indir1, "data", "poverty_stratified_censoring", "mort_dt_ls_pov.rds"))

#Set exposure
exposure_ <- "mean_daily_peak_week"

f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = mort_dt_ls_pov, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=mort_dt_ls_pov, family = binomial(), sparse=FALSE)

mort_dt_ls_pov$ps.censoring1_mwfpm <- predict(m.censoring, type="response", mort_dt_ls_pov)

mort_dt_ls_pov$ps.censoring2_mwfpm <- predict(m.censoring.empty, type="response", mort_dt_ls_pov)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
mort_dt_ls_pov$w_c_i_mwfpm <- 1/mort_dt_ls_pov$ps.censoring1_mwfpm
mort_dt_ls_pov<- mort_dt_ls_pov %>%
  group_by(studyid) %>% 
  mutate(w_c_mwfpm=cumprod(w_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(w_c_mwfpm=case_when(censor==1 ~ 0,
                             .default=w_c_mwfpm))
summary(mort_dt_ls_pov$w_c_mwfpm)


## stablized weights
mort_dt_ls_pov$sw_c_i_mwfpm <- mort_dt_ls_pov$ps.censoring2_mwfpm/mort_dt_ls_pov$ps.censoring1_mwfpm
mort_dt_ls_pov<- mort_dt_ls_pov %>%
  group_by(studyid) %>% 
  mutate(sw_c_mwfpm=cumprod(sw_c_i_mwfpm)) %>%
  ungroup() %>%
  mutate(sw_c_mwfpm=case_when(censor==1 ~ 0,
                              .default=sw_c_mwfpm))
summary(mort_dt_ls_pov$sw_c_mwfpm)

mort_dt_ls_pov %>%
  saveRDS(here(indir1, "data", "poverty_stratified_censoring", "mort_dt_ls_pov.rds"))
####################

## run pooled logistic regression
####################
## below are variable names that should be manually entered
tf_cv <- c("gender","mrtl_status","r_e","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv <- c("pop_den","year", "age_boq_s") ## variable names for time-varying covariates 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"

## Adjusted mean wildfire PM modeling for 15% or more poverty
#Bring in data
mort_dt_pov<-subset(readRDS(here(indir1, "data", "poverty_stratified_censoring","mort_dt_pov.rds")),
                    select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                             censor, age_boq_s ,mrtl_status,r_e, rqrs_intrprtr,
                             smkn_status,pov_p,pop_den,sw_c_mwfpm))


f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
p_adj <- glm(f2_adj, data = mort_dt_pov, family = "binomial", weights=sw_c_mwfpm)
summary(p_adj)
p_adj %>%
  glance()
p_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))

## Adjusted mean wildfire PM modeling for less than 15% poverty 
#Bring in data
mort_dt_ls_pov<-subset(readRDS(here(indir1, "data", "poverty_stratified_censoring", "mort_dt_ls_pov.rds")),
                       select=c(studyid,obs_num,year, death_nw, mean_wf_pm, gender, 
                                censor, age_boq_s ,mrtl_status,r_e, rqrs_intrprtr,
                                smkn_status,pov_p,pop_den,sw_c_mwfpm))

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
lp_adj <- glm(f2_adj, data = mort_dt_ls_pov, family = "binomial", weights=sw_c_mwfpm)
summary(lp_adj)
lp_adj %>%
  glance()
lp_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))
####################
