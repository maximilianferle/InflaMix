rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/000_functions_constants.R")

library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)
library(survcomp)
library(pROC)
library(cowplot)

# Select all patients not in the derivation cohort
df_all_chrt <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>%
  # To evaluate InflaMix validation with just six labs (albumin, hgb, crp, ldh, alp, ast) -
  # uncomment the line below, which ignores all other labs (This will give figures 3d-l)
  # mutate_at(vars(starts_with(focus_labs)), ~.*NA) %>%
  filter(cohort!="MSK Development") %>%
  mutate(record_id=record_id)

# Master data frame of all validation cohort data
val_df <- reformat_table_tp(df_all_chrt, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  rename(Cluster=cluster) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  mutate(preld_crpuln_ratio=noLCS_preld_crp/uln_preld_crp) %>%
  mutate(preld_ferruln_ratio=noLCS_preld_ferritin/uln_preld_ferritin) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# For each of the following validation analyses, define the cohort, calculate inference metrics
# Note that the functions will evaluate various inferences that are not relevant or of interest.
# Table 2 of the paper clarifies the relevant analyses.

#### II. MSK LBCL Validation
coh2_df <- val_df %>% filter(cohort=="MSK Validation")
coh2_metrics <- clust_metrics(df=coh2_df, dxres=lymphoma_dxres,
              covar=lbcl_covar,res_covar=lbcl_covar, wght=TRUE, aname="msk_validation")


##### III. SMC+HMH LBCL Validation
coh3_df <- val_df %>% filter(cohort=="Center Validation")
coh3_metrics <- clust_metrics(df=coh3_df, dxres=lymphoma_dxres,
                                 covar=lbcl_covar, res_covar=lbcl_covar, wght=TRUE, aname="center_validation")

#### Cohort IV. MCL + FL Validation
coh4_df <- val_df %>% filter(cohort=="Non-Hodgkin Lymphoma Validation")

coh4_metrics <- clust_metrics(df=coh4_df, dxres=lymphoma_dxres,
                                 covar=lymphoma_covar,res_covar=lymphoma_covar, wght=TRUE, aname="nhl_validation")

##############################
# Exploratory Analyses by costim domain
df_all_chrt_plusdev <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes")

all_chrt_df <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  rename(Cluster=cluster) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# All patients treated with CD28 CAR-T
cd28_df <- all_chrt_df %>% filter(costim=="CD28")
cd28_metrics <- clust_metrics(df=cd28_df, dxres=lymphoma_dxres,
                                 covar=costim_covar, res_covar=costim_covar, wght=TRUE, aname="cd28_validation")

# All patients treated with 41BB CAR-T
bb41_df <- all_chrt_df %>% filter(costim=="41BB")
bb41_metrics <- clust_metrics(df=bb41_df, dxres=lymphoma_dxres,
                                 covar=costim_covar, res_covar=costim_covar, wght=TRUE, aname="41BB_validation")


# Exploratory Analyses by MTV domain
df_all_chrt_plusdev <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes")

all_chrt_df <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  rename(Cluster=cluster) %>%
  filter(dx_simple=="LBCL") %>% 
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(bl_mtv)) %>%
  mutate(record_id=record_id)


# MVA with MTV
mtv_metrics <- clust_metrics(df=all_chrt_df, dxres=lymphoma_dxres,
                                covar="+age+costim+primary_ref+bl_mtv+bl_mtv*Cluster", 
                                res_covar="+age+costim+primary_ref+bl_mtv+bl_mtv*Cluster", wght=TRUE, aname="mtv_validation")

mtv_cut <- quantile( (df_all %>% filter(dx_simple=="LBCL" & !is.na(bl_mtv)))$bl_mtv, probs=c(0.33, 0.5, 0.66))[3]

# Subgroup High MTV
hightb_df <- all_chrt_df %>% filter(bl_mtv>mtv_cut)
hightb_metrics <- clust_metrics(df=hightb_df, dxres=lymphoma_dxres,
                              covar="+age+costim+primary_ref", 
                              res_covar="+age+costim+primary_ref", wght=TRUE, aname="hightb_validation")

# Subgroup Low MTV
lowtb_df <- all_chrt_df %>% filter(bl_mtv<mtv_cut)
lowtb_metrics <- clust_metrics(df=lowtb_df, dxres=lymphoma_dxres,
                                covar="+age+costim+primary_ref", 
                                res_covar="+age+costim+primary_ref", wght=TRUE, aname="lowtb_validation")


#### 
validation_inf <- bind_rows(
  coh2_metrics,
  coh3_metrics,
  coh4_metrics,
  cd28_metrics,
  bb41_metrics,
  mtv_metrics,
  hightb_metrics,
  lowtb_metrics,
) 


# Use the below code to plot PFS and OS KM survival estimate curves
# Fig 2c, d, f, g, i, j
# Fig 3d, e, g, h, j, k
df_analysis = coh4_df
anval = "nhl_validation"
tmax=25
qm = 5
sz1=20
pxm = 15
pym = 0.9
oxm = 7
oym = 0.13

pfs <- clus_ggsurvplot(df=df_analysis,wght=FALSE,sz=sz1,
                        xmi= pxm,
                        ymi= pym,
                        event=tolower("PFS"),
                        metric=validation_inf %>% filter(analysis==anval & outcome=="PFS"),
                        labadj=0.05,
                        labl="")

os <- clus_ggsurvplot(df=df_analysis,wght=FALSE,sz=sz1,
                xmi= oxm,
                ymi= oym,
                event=tolower("OS"),
                metric=validation_inf %>% filter(analysis==anval & outcome=="OS"),
                labadj=0.05,
                labl="")


# Use the below code to plot CR bar plots
# Fig 2e, h, k
# Fig 3f, i, l
m1 = validation_inf %>% filter(analysis==anval & outcome=="No CR")
sz1=10
xm = 1
ym = 0.98
labl = ""

cr <- plot_grid((clust_bar_plot(df=df_analysis, tp_pre="", list_res=c("CR"), sz=sz1/1.8, xmi=xm, ymi=ym, metric=m1) + ylim(0, 1.1) ), labels=c(labl), scale=1, label_size = 30)

plot_grid(ggsurvfit_build(pfs), ggsurvfit_build(os), cr, nrow=1, labels=c("", "", ""), rel_widths = c(2,2, 1)) 


