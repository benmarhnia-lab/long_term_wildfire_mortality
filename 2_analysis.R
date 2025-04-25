####################
# Title: 2_main_analysis.R
# Date Created: 4/17/2025
# Author: Tim B. Frankland; compiled by Chen Chen
# Purpose: conduct main analysis and nonlinear analysis for the long-term wildfire smoke 
# on mortality KPSC project. To note, we only showed codes for the metric average 
# wildfire PM2.5 (mean_wf_pm), with differences for other metrics marked by "To note"
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
library(mgcv)
library(dlnm)

indir1 <- "" ## work directory

# Data used in this code:
# dth_dt_cln.rds includes all health and exposure data, with one record for each quarter-person.
# dth_wf_exp.csv includes birthday of each individual 

####################
## Main analysis
####################

## calculate inverse probability censoring weights for the entire population
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_cln.rds")),
               select=c(studyid,obs_num,year, death,  mean_wf_pm, censoring, 
                        gender,age ,mrtl_status,r_e, rqrs_intrprtr,smkn_status,
                        pov_p,pop_den))


#IQR
IQR(dt_cln$mean_wf_pm) 

dt_cln<-dt_cln %>%
  mutate(mean_wf_pm_iqr=IQR(mean_wf_pm), ## To note, for weeks_gt_5 and smoke_waves, we did not do this transformation and ran the model with original exposure unit
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_wf_pm_iqr))

dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0) 


dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

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

dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 


rm(dt_cln0)

## below are variable names that should be manually entered
tf_cv <- c("gender", "age","mrtl_status","r_e","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv <- c("pov_p","pop_den","year") ## variable names for time-varying covariates 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"


## Censoring Weight Calculations 
f.censoring <- reformulate(c(tf_cv, tv_cv_0, tv_cv,
                             paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring <- speedglm(f.censoring, data = dt_cln, family = binomial(), sparse=FALSE)

f.censoring.empty <- reformulate(c(tf_cv, tv_cv_0,
                                   paste0(exposure_, "* as.factor(", time_variable, ")")
), response = "censor==0")

m.censoring.empty <- speedglm(f.censoring.empty, data=dt_cln, family = binomial(), sparse=FALSE)

dt_cln$ps.censoring1 <- predict(m.censoring, type="response", dt_cln)

dt_cln$ps.censoring2 <- predict(m.censoring.empty, type="response", dt_cln)

rm(m.censoring)
rm(m.censoring.empty)

## unstablized weights
dt_cln$w_c_i <- 1/dt_cln$ps.censoring1
dt_cln<- dt_cln %>%
  group_by(studyid) %>% 
  mutate(w_c=cumprod(w_c_i)) %>%
  ungroup() %>%
  mutate(w_c=case_when(censor==1 ~ 0,
                       .default=w_c))
summary(dt_cln$w_c)


## stablized weights
dt_cln$sw_c_i <- dt_cln$ps.censoring2/dt_cln$ps.censoring1
dt_cln<- dt_cln %>%
  group_by(studyid) %>% 
  mutate(sw_c=cumprod(sw_c_i)) %>%
  ungroup() %>%
  mutate(sw_c=case_when(censor==1 ~ 0,
                        .default=sw_c))
summary(dt_cln$sw_c)

dt_cln %>%
  dplyr::select(studyid,obs_num,w_c, sw_c) %>%
  saveRDS(here(indir1, "data", "censor_weights_mwfpm.rds"))
####################

## Create Time-Varying Age Dataset
####################
dt<-read_csv(here(indir1, "data", "dth_wf_exp.csv"))

dt_bd<-dt %>%
  dplyr::select(studyid,birthdate) %>%
  group_by(studyid) %>%
  slice(1) %>%
  ungroup() 

dt_cln<-readRDS(here(indir1, "data", "dth_dt_cln.rds"))

dt_cln_tvage<-dt_cln %>%
  left_join(dt_bd,by=join_by(studyid)) %>%
  mutate(birthdate=as_date(mdy(birthdate)))%>%
  mutate(qtr_date=case_when(quarter==1 ~ make_date(yr,1,1),
                            quarter==2 ~ make_date(yr,4,1),
                            quarter==3 ~ make_date(yr,7,1),
                            quarter==4 ~ make_date(yr,10,1)),
         age_beg_of_qtr=year(as.period(interval(start=birthdate,end=qtr_date))),
         age_boq_s=ns(age_beg_of_qtr,df=2))

saveRDS(dt_cln_tvage,here(indir1, "data", "dth_dt_expos_tvage.rds"))
####################

## run pooled logistic model
####################
#Bring in data
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death,  mean_wf_pm, censoring, 
                        gender,age_boq_s ,mrtl_status,r_e,r_e_full, rqrs_intrprtr,
                        smkn_status,pov_p,pop_den))


# IQR
IQR(dt_cln$mean_wf_pm)

dt_cln<-dt_cln %>%
  mutate(mean_wf_pm_iqr=IQR(mean_wf_pm),
         mean_wf_pm=mean_wf_pm/mean_wf_pm_iqr,
         pop_den=pop_den/1000) %>%
  dplyr::select(-c(mean_wf_pm_iqr))


dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0) 

dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0),
         r_e=case_when(r_e!="Other" ~ r_e,
                       r_e=="Other" & r_e_full=="Unknown" ~ "Other-Unknown",
                       r_e=="Other" & r_e_full!="Unknown" ~ "Other-Known"))

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

dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

cw<-readRDS(here(indir1, "data", "censor_weights_mwfpm.rds"))

dt_cln<-dt_cln %>%
  inner_join(cw,by=join_by(studyid,obs_num))

rm(dt_cln0)

## below are variable names that should be manually entered
tf_cv <- c("gender", "age_boq_s","mrtl_status","r_e","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv <- c("pov_p","pop_den","year") ## variable names for time-varying covariates 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"

f2_adj <- reformulate(c(tf_cv,tv_cv,exposure_), response = "death_nw")
m2_adj <- glm(f2_adj, data = dt_cln, family = "binomial", weights = sw_c)
summary(m2_adj)
m2_adj %>%
  glance()
m2_adj %>% 
  tbl_regression(exponentiate = TRUE, include=mean_wf_pm,
                 estimate_fun = purrr::partial(style_ratio, digits = 4),
                 pvalue_fun = purrr::partial(style_sigfig, digits = 4))
####################

####################
## Nonlinear analysis
####################

#Bring in data
####################
dt_cln<-subset(readRDS(here(indir1, "data", "dth_dt_expos_tvage.rds")),
               select=c(studyid,obs_num,year, death,  mean_wf_pm, censoring, gender,age_boq_s ,mrtl_status,r_e, rqrs_intrprtr,smkn_status,pov_p,pop_den))


dt_cln<-dt_cln %>%
  mutate(pop_den=pop_den/1000)

dt_cln0<-dt_cln %>%
  dplyr::select(studyid,obs_num,year,pov_p,pop_den) %>%
  group_by(studyid) %>%
  filter(obs_num==13) %>%
  mutate(pov_p_0=pov_p,
         pop_den_0=pop_den,
         year_0=year) %>%
  ungroup() %>%
  dplyr::select(studyid,pov_p_0,pop_den_0,year_0) 


dt_cln<-dt_cln %>%
  left_join(dt_cln0,by="studyid") %>%
  mutate(censor=case_when(death==0 & censoring==1 ~ 1, 
                          .default = 0),
         death_nw=case_when(death==0 & censoring==1 ~ NA,
                            death==1 ~ 1,
                            .default=0))

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

dt_cln<-dt_cln %>%
  group_by(studyid) %>%
  mutate(study_time_period=row_number()) %>%
  ungroup() 

cw<-readRDS(here(indir1, "results", "censor_weights_mwfpm.rds"))

dt_cln<-dt_cln %>%
  inner_join(cw,by=join_by(studyid,obs_num))

dt_nw<-dt_cln %>%
  filter(!is.na(death_nw))

rm(dt_cln0)

## below are variable names that should be manually entered
tf_cv <- c("gender", "mrtl_status","r_e","rqrs_intrprtr","smkn_status") ## variable names for time-fixed covariates
tv_cv <- c("pov_p","pop_den","year","age_boq_s") ## variable names for time-varying covariates 
tv_cv_0<- paste(tv_cv, "0", sep="_") 
exposure_<-c("mean_wf_pm")
time_variable <- "study_time_period"
####################

## run nonlinear model using gam
####################
m2_adj<-gam(death_nw~s(mean_wf_pm,k=5)+gender+age_boq_s+mrtl_status+r_e+rqrs_intrprtr+
              smkn_status+pov_p+pop_den+year,data=dt_cln,method="REML",family="binomial", weights=sw_c)
gam.check(m2_adj)
summary(m2_adj)
m2_adj %>%
  glance()
plot(m2_adj)


predict<-predict(m2_adj, data=dt_nw,type="terms", se.fit=T)

wf_pm_pred <- cbind(dt_nw$mean_wf_pm, predict$fit[,"s(mean_wf_pm)"], predict$se.fit[,"s(mean_wf_pm)"]) 
wf_pm_pred <- data.frame(wf_pm_pred[order(wf_pm_pred[,1]),])
names(wf_pm_pred) <- c("mean_wf_pm","fit","sefit")
wf_pm_pred$rr <- exp(wf_pm_pred$fit) ##odds of mortality
wf_pm_pred$ll <- exp(wf_pm_pred$fit-(1.96*wf_pm_pred$sefit))
wf_pm_pred$ul <- exp(wf_pm_pred$fit+(1.96*wf_pm_pred$sefit))


## predict log odds of mortality at 0 exposure value
predict.dlnm <- crosspred(basis="mean_wf_pm", model = m2_adj, cen=0, at=wf_pm_pred$mean_wf_pm)
plot(predict.dlnm, ylab="Odds ratio", xlab="mean_wf_pm", col=2)

## Put predicted values from crosspred function into data frame and rename columsn
pred_dlnm<-as_data_frame(cbind(predict.dlnm$predvar,predict.dlnm$allRRfit,predict.dlnm$allRRlow,predict.dlnm$allRRhigh))
pred_dlnm<-pred_dlnm %>%
  rename(mean_wf_pm=V1,or=V2,or.ll=V3,or.ul=V4)

## Join to wf_pm_pred by exposure values
wf_pm_pred<-wf_pm_pred %>%
  left_join(pred_dlnm,by=join_by(mean_wf_pm))

wf_pm_pred %>%
  saveRDS(here(indir1, "results", "wf_pm_pred_mort_dlmn.rds"))
####################
