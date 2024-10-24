rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/000_functions_constants.R")

library(labelled)
library(survival)
library(gtsummary)
library(gt)
library(webshot2)

tbl1_df <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>%
  mutate(preld_crpuln_ratio=noLCS_preld_crp/uln_preld_crp) %>%
  select(
    record_id,
    center,
    cohort,
    age,
    kps90,
    dx_simple,
    car_t_product,
    stage_aph,
    trans_nhl,
    primary_ref,
    bin_preld_ldh,
    bridge_categories,
    num_rxline_3lvl,
    crs24,
    icans24,
    everCR_100,
    ev_os,
    ev_pfs,
    tt_os_m,
    tt_pfs_m,
    preld_crpuln_ratio
  ) %>%
  mutate(cohort=ifelse(cohort=="MSK Development", "I. MSK LBCL", ifelse(cohort=="MSK Validation", "II. MSK LBCL",
                                                                        ifelse(cohort=="Center Validation", "III. SMC and HMH LBCL", "IV. MCL and FL Lymphoma")))) %>%
  mutate(cohort=factor(as.character(cohort), levels=c(
    "I. MSK LBCL",
    "II. MSK LBCL",
    "III. SMC and HMH LBCL",
    "IV. MCL and FL Lymphoma"
  )))  %>%
  mutate(dx_simple=case_when(
    dx_simple=="LBCL" ~ "Large B-cell Lymphoma",
    dx_simple=="MCL" ~ "Mantle Cell Lymphoma",
    dx_simple=="FL" ~ "Follicular Lymphoma",
    .default=NA
    )) %>%
  mutate(center=factor(as.character(center), levels=c(
    "Memorial Sloan Kettering",
    "Sheba Medical Center",
    "Hackensack Meridian Health"
  ))) %>%
  mutate(dx_simple=factor(as.character(dx_simple), levels=c(
    "Large B-cell Lymphoma",
    "Follicular Lymphoma",
    "Mantle Cell Lymphoma"
  ))) %>%
  mutate(car_t_product=factor(as.character(car_t_product), levels=c(
    "Axicabtagene ciloleucel",
    "Point-Of-Care CD19-directed \n CD28-costimulatory domain CAR-T",
    "Tisagenlecleucel",
    "Lisocabtagene maraleucel",
    "Brexucabtagene autoleucel"
  )))

var_label(tbl1_df$center) <- "Treatment Center"
var_label(tbl1_df$age) <- "Age"
var_label(tbl1_df$kps90) <- "Karnofsky Performance Status (KPS)"
var_label(tbl1_df$car_t_product) <- "CAR-T Product"
var_label(tbl1_df$bridge_categories) <- "Bridging Therapy"
var_label(tbl1_df$dx_simple) <- "Disease"
var_label(tbl1_df$primary_ref) <- "Primary Refractory Disease"
var_label(tbl1_df$stage_aph) <- "Stage at Apheresis"
var_label(tbl1_df$num_rxline_3lvl) <- "Num. Prior Lines of Therapy"
var_label(tbl1_df$trans_nhl) <- "Transformed Lymphoma"
var_label(tbl1_df$crs24) <- "CRS Grade 2-4"
var_label(tbl1_df$icans24) <- "ICANS Grade 2-4"
var_label(tbl1_df$everCR_100) <- "CR by Day 100"
var_label(tbl1_df$cohort) <- "Cohort"
var_label(tbl1_df$tt_os_m) <- "Overall Survival (months)"
var_label(tbl1_df$tt_pfs_m) <- "Progression Free Survival (months)"
var_label(tbl1_df$bin_preld_ldh) <- "Prelymphodepletion LDH above ULN"
var_label(tbl1_df$preld_crpuln_ratio) <- "Prelymphodepletion CRP / ULN"

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

tbl1c1 <- tbl_survfit(
  tbl1_df %>% select(cohort, tt_os_m, tt_pfs_m, ev_os, ev_pfs, cohort),
  y=Surv(tt_os_m, ev_os==0),
  include=c(cohort),
  probs=0.5,
  label_header = "**Median Follow-up**"
)

tbl1b <- tbl_merge(list(tbl1c1), tab_spanner = c("")) %>%
  modify_footnote(
    all_stat_cols() ~ "Median (95% Confidence Interval)"
  ) %>% as_gt() %>% as_raw_html()

tbl1 <- tbl_summary(tbl1_df %>% select(!c(record_id, tt_pfs_m, tt_os_m, ev_os, ev_pfs)), by="cohort") %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1") ~ "**Model Derivation Cohort**") %>%
  modify_spanning_header(c("stat_2", "stat_3", "stat_4") ~ "**Validation Cohorts**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("****") %>%
  bold_labels() %>%
  as_gt() %>%
  #gt::tab_source_note(gt::html(tbl1b)) %>%
  gt::tab_source_note(gt::md("*FL, Follicular Lymphoma; LBCL, Large B-cell Lymphoma; MCL, Mantle Cell Lymphoma; MSK, Memorial Sloan Kettering Cancer Center; SMC, Sheba Medical Center; HMH, Hackensack Meridian Health; CRS, Cytokine Release Syndrome; ICANS, Immune Effector Cell-Associated Neurotoxicity Syndrome; CR, Complete Response; LDH, Lactate Dehydrogenase; ULN, Upper Limit of Normal*")) %>%
  tab_options(table.font.size = 9.75)


gt::gtsave(tbl1, file="tables_figures/tbl1.pdf")


tbl1 <- tbl_summary(tbl1_df %>% select(!c(record_id, tt_pfs_m, tt_os_m, ev_os, ev_pfs)), by="cohort") %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1") ~ "**Model Derivation Cohort**") %>%
  modify_spanning_header(c("stat_2", "stat_3", "stat_4") ~ "**Validation Cohorts**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  bold_labels()

tbl1b <- tbl_merge(list(tbl1c1), tab_spanner = c("")) %>%
  modify_footnote(
    all_stat_cols() ~ "Median (95% Confidence Interval)"
  )

survfit(Surv(tt_os_m, ev_os==0)~1, data=tbl1_df)


### Supplementary Table 1
df_dev <-  left_join(df_all, df_labs_all %>% select(!cohort) %>% select(record_id, starts_with("d0_")), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>% # Must have curated metadata
  filter(cohort=="MSK Development") %>%
  select(record_id, contains(cluster_labs))

inflamix_assignments <- reformat_table_tp(df=df_dev, gmm_mu = dev_mu, gmm_sig = dev_sig, gmm_pro = dev_pro, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp = "d0_") %>%
  group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>% left_join(df_all, by="record_id")

supptbl1_df <- inflamix_assignments %>%
  transmute(
    age, 
    kps90,
    primary_ref,
    trans_nhl,
    stage_aph,
    num_rxline_3lvl,
    costim,
    bridge_categories,
    cluster
  )

var_label(supptbl1_df$age) <- "Age"
var_label(supptbl1_df$kps90) <- "Karnofsky Performance Status (KPS)"
var_label(supptbl1_df$costim) <- "Costimulatory Domain"
var_label(supptbl1_df$bridge_categories) <- "Bridging Therapy"
var_label(supptbl1_df$primary_ref) <- "Primary Refractory Disease"
var_label(supptbl1_df$stage_aph) <- "Stage at Apheresis"
var_label(supptbl1_df$num_rxline_3lvl) <- "Num. Prior Lines of Therapy"
var_label(supptbl1_df$trans_nhl) <- "Transformed Lymphoma"
var_label(supptbl1_df$cluster) <- "Cluster"

supp_tbl1 <- tbl_summary(supptbl1_df, by="cluster") %>%
  add_overall() %>%
  add_p() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  #modify_caption("****") %>%
  bold_labels()


#### Supp Table 3
# Filter patients, analysis_type==1 signifies complete clinical metadata
df_all_chrt_plusdev <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id") %>%filter(analysis_type=="Lymphoma Outcomes") %>%
  mutate(record_id=record_id)

# Assign clusters using InflaMix at d0, pre-lymphodepletion, and pre-apheresis timepoints.
df_all_chrt_d0 <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  rename(d0_tau=tau)
df_all_chrt_preld <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preld_") %>%
  reframe(record_id, cluster=cluster, preld_tau=tau)
df_all_chrt_preaph <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preaph_") %>%
  reframe(record_id, cluster=cluster, preaph_tau=tau)

# Merge cluster assignment probabilities together in wide format
df <- df_all_chrt_d0 %>%
  left_join(df_all_chrt_preaph, by=c("record_id", "cluster")) %>%
  left_join(df_all_chrt_preld, by=c("record_id", "cluster")) %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(bridge_yn)) %>%
  filter(!is.na(primary_ref)) %>%
  filter(!is.na(car_t_product)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# Restructuring data frame and selecting only the most probable cluster assignment with its associated probability
prealluv <- df %>% group_by(record_id) %>% dplyr::slice(which.max(preaph_tau)) %>% ungroup() %>%
  select(!contains("tau")) %>% mutate(preaph_cluster=cluster) %>% select(!cluster) %>%
  left_join(df %>% select(record_id, cluster, preld_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(preld_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(preld_cluster=cluster) %>% select(!cluster), by="record_id") %>%
  left_join(df %>% select(record_id, cluster, d0_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(d0_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(d0_cluster=cluster) %>% select(!cluster), by="record_id") %>%
  column_to_rownames("record_id")

supptbl3_df <- prealluv %>%
  filter(!is.na(preaph_cluster)) %>%
  filter(!is.na(d0_cluster)) %>%
  mutate(
    h1_cluster=ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Inflammatory", "Infl. -> Infl.",
                      ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Non-Inflammatory","Infl. -> Non-Infl.",
                             ifelse(preaph_cluster=="Non-Inflammatory" & d0_cluster=="Inflammatory","Non-Infl. -> Infl.",
                                    "Non-Infl. -> Non-Infl.")))) %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Infl. -> Infl.","Infl. -> Non-Infl.", "Non-Infl. -> Non-Infl.", "Non-Infl. -> Infl."))) %>%
  mutate(h1_cluster=h1_cluster) %>% mutate(Cluster=h1_cluster) %>% rownames_to_column("record_id") %>%
  transmute(
    age, 
    kps90,
    dx_simple,
    primary_ref,
    trans_nhl,
    stage_aph,
    num_rxline_3lvl,
    car_t_product,
    bridge_categories,
    vvtime, 
    noLCS_preld_ldh,
    h1_cluster
  )

var_label(supptbl3_df$age) <- "Age"
var_label(supptbl3_df$kps90) <- "Disease"
var_label(supptbl3_df$dx_simple) <- "Karnofsky Performance Status (KPS)"
var_label(supptbl3_df$car_t_product) <- "CAR-T Product"
var_label(supptbl3_df$bridge_categories) <- "Bridging Therapy"
var_label(supptbl3_df$primary_ref) <- "Primary Refractory Disease"
var_label(supptbl3_df$stage_aph) <- "Stage at Apheresis"
var_label(supptbl3_df$num_rxline_3lvl) <- "Num. Prior Lines of Therapy"
var_label(supptbl3_df$trans_nhl) <- "Transformed Lymphoma"
var_label(supptbl3_df$vvtime) <- "Vein-to-vein Time (Days)"
var_label(supptbl3_df$noLCS_preld_ldh) <- "Pre-Lymphodepletion LDH"
var_label(supptbl3_df$h1_cluster) <- "Cluster Transition"

supp_tbl3 <- tbl_summary(supptbl3_df, by="h1_cluster") %>%
  add_overall() %>%
  #add_p() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  #modify_caption("****") %>%
  bold_labels()

supptbl3ld_df <- prealluv %>%
  filter(!is.na(preld_cluster)) %>%
  filter(!is.na(d0_cluster)) %>%
  mutate(
    h1_cluster=ifelse(preld_cluster=="Inflammatory" & d0_cluster=="Inflammatory", "Infl. -> Infl.",
                      ifelse(preld_cluster=="Inflammatory" & d0_cluster=="Non-Inflammatory","Infl. -> Non-Infl.",
                             ifelse(preld_cluster=="Non-Inflammatory" & d0_cluster=="Inflammatory","Non-Infl. -> Infl.",
                                    "Non-Infl. -> Non-Infl.")))) %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Infl. -> Infl.","Infl. -> Non-Infl.", "Non-Infl. -> Non-Infl.", "Non-Infl. -> Infl."))) %>%
  mutate(h1_cluster=h1_cluster) %>% mutate(Cluster=h1_cluster) %>% rownames_to_column("record_id") %>%
  transmute(
    age, 
    kps90,
    dx_simple,
    primary_ref,
    trans_nhl,
    stage_aph,
    num_rxline_3lvl,
    car_t_product,
    #bridge_categories,
    #vvtime, 
    noLCS_preld_ldh,
    h1_cluster
  )

var_label(supptbl3ld_df$age) <- "Age"
var_label(supptbl3ld_df$kps90) <- "Disease"
var_label(supptbl3ld_df$dx_simple) <- "Karnofsky Performance Status (KPS)"
var_label(supptbl3ld_df$car_t_product) <- "CAR-T Product"
#var_label(supptbl3ld_df$bridge_categories) <- "Bridging Therapy"
var_label(supptbl3ld_df$primary_ref) <- "Primary Refractory Disease"
var_label(supptbl3ld_df$stage_aph) <- "Stage at Apheresis"
var_label(supptbl3ld_df$num_rxline_3lvl) <- "Num. Prior Lines of Therapy"
var_label(supptbl3ld_df$trans_nhl) <- "Transformed Lymphoma"
#var_label(supptbl3ld_df$vvtime) <- "Vein-to-vein Time (Days)"
var_label(supptbl3ld_df$noLCS_preld_ldh) <- "Pre-Lymphodepletion LDH"
var_label(supptbl3ld_df$h1_cluster) <- "Cluster Transition"

supp_tbl3ld <- tbl_summary(supptbl3ld_df, by="h1_cluster") %>%
  add_overall() %>%
  #add_p() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  #modify_caption("****") %>%
  bold_labels()
