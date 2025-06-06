rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/000_functions_constants.R")

library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggsci)
library(pheatmap)

# Patients from all cohorts (metadata + lab data)
df_all_chrt_plusdev <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id")

# Apply InflaMix to predict clusters at day 0, pre-lymphodepletion, and pre-apheresis
df_all_chrt_d0 <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  rename(d0_tau=tau)
df_all_chrt_preld <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preld_") %>%
  reframe(record_id, cluster=cluster, preld_tau=tau)
df_all_chrt_preaph <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preaph_") %>%
  reframe(record_id, cluster=cluster, preaph_tau=tau)

# Combine cluster probabilities at all timepoints in wide format.
df_all_chrt <- df_all_chrt_d0 %>%
  left_join(df_all_chrt_preaph, by=c("record_id", "cluster")) %>%
  left_join(df_all_chrt_preld, by=c("record_id", "cluster")) %>%
  mutate(record_id=record_id)

df1 <- df_all_chrt

# Condense the data frame to only consider only cluster assignment and probability of assignment to that cluster across different timepoints.
df_hc <- df1 %>% group_by(record_id) %>% dplyr::slice(which.max(d0_tau)) %>% ungroup() %>%
  select(!contains("tau")) %>% mutate(d0_cluster=cluster) %>% select(!cluster) %>%
  left_join(df1 %>% select(record_id, cluster, preld_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(preld_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(preld_cluster=cluster) %>% select(!cluster), by="record_id") %>%
  left_join(df1 %>% select(record_id, cluster, preaph_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(preaph_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(preaph_cluster=cluster) %>% select(!cluster), by="record_id")


df_all_chrt1 <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id")



#####@ Total missing labs  - Figure 3a ########

# Get the number of patient of patients per cohort
numcohort <- df_hc %>% filter(analysis_type=="Lymphoma Outcomes") %>% select(starts_with("d0_"), cohort) %>% select(contains(cluster_labs), cohort) %>%
  group_by(cohort) %>% summarize(n=n()) %>% ungroup

# Calculate the percent of patients missing each lab by cohort.
missdat_cohort <- df_hc %>% filter(analysis_type=="Lymphoma Outcomes") %>% select(starts_with("d0_"), cohort) %>% select(contains(cluster_labs), cohort) %>%
  pivot_longer(!cohort, names_to = "labs", values_to = "value") %>%
  group_by(cohort, labs) %>% summarize(nacount=sum(is.na(value))) %>% ungroup() %>%
  filter(!grepl("Development", cohort)) %>%
  left_join(numcohort, by="cohort") %>% mutate(nacount=round(nacount/n,2)) %>% select(!n) %>%
  pivot_wider(names_from = "cohort", values_from = nacount) %>%
  column_to_rownames("labs") %>%
  as.matrix()

# Generate a heatmap of missing labs by cohort (colorscale is percent missing labs)
dfc1 <- missdat_cohort

rownames(dfc1)[rownames(dfc1) == 'd0_albumin'] <- "Albumin"
rownames(dfc1)[rownames(dfc1) == 'd0_wbc'] <- "WBC"
rownames(dfc1)[rownames(dfc1) == 'd0_alk'] <- "ALP"
rownames(dfc1)[rownames(dfc1) == 'd0_ast'] <- "AST"
rownames(dfc1)[rownames(dfc1) == 'd0_crp'] <- "CRP"
rownames(dfc1)[rownames(dfc1) == 'd0_ddimer'] <- "Ddimer"
rownames(dfc1)[rownames(dfc1) == 'd0_ferritin'] <- "Ferritin"
rownames(dfc1)[rownames(dfc1) == 'd0_hb'] <- "Hgb"
rownames(dfc1)[rownames(dfc1) == 'd0_il10'] <- "IL-10"
rownames(dfc1)[rownames(dfc1) == 'd0_il6'] <- "IL-6"
rownames(dfc1)[rownames(dfc1) == 'd0_ldh'] <- "LDH"
rownames(dfc1)[rownames(dfc1) == 'd0_plt'] <- "Plts"
rownames(dfc1)[rownames(dfc1) == 'd0_tbr'] <- "Tbili"
rownames(dfc1)[rownames(dfc1) == 'd0_tnfa'] <- "TNFa"

df_hm <- dfc1
clusheat <- function(x) {
  fig3a <- pheatmap::pheatmap(x, fontsize_row=18, fontsize_col = 18,
                              cluster_cols = FALSE,
                              legend=TRUE,

                              labels_col = c(paste0("MSK LBCL\nn=", numcohort$n[2]), paste0("Center LBCL\nn=", numcohort$n[3]), paste0("NHL\nn=", numcohort$n[4])),
                              treeheight_row = 0, fontsize = 20,
                              cluster_rows = TRUE,
                              clustering_distance_rows = "euclidean",
                              display_numbers = x,
                              angle_col = 0,
                              number_color="black"
  )
}

# Supp Fig 4b
sfig4b_main <- clusheat(df_hm)


# Calculate the number of patients missing n labs where n = 1:16
missdat <-  bind_rows(
  df_hc  %>% filter(analysis_type=="Lymphoma Outcomes") %>% select(record_id, starts_with("d0_"), cohort) %>% select(record_id, contains(cluster_labs)) %>%
    pivot_longer(!record_id, names_to = "labs", values_to = "value") %>%
    group_by(record_id) %>% summarize(nacount=sum(is.na(value))) %>% ungroup() %>% mutate(timepoint="d0"),
  df_hc %>% filter(analysis_type=="Lymphoma Outcomes")%>% select(record_id, starts_with("preld_")) %>% select(record_id, contains(cluster_labs)) %>%
    pivot_longer(!record_id, names_to = "labs", values_to = "value") %>%
    group_by(record_id) %>% summarize(nacount=sum(is.na(value))) %>% ungroup() %>% mutate(timepoint="preld"),
  df_hc %>% filter(analysis_type=="Lymphoma Outcomes")%>% select(record_id, starts_with("preaph_")) %>% select(record_id, contains(cluster_labs)) %>%
    pivot_longer(!record_id, names_to = "labs", values_to = "value") %>%
    group_by(record_id) %>% summarize(nacount=sum(is.na(value))) %>% ungroup() %>% mutate(timepoint="preaph")
) %>% left_join(df_hc %>% select(record_id, dx=dx_simple, cohort), by="record_id") %>%
  filter(nacount >=0 & nacount <= 16)

missdat2 <- missdat %>%
  mutate(timepoint=ifelse(timepoint=="d0", paste0("Peri-Infusion (n=",dim(unique(missdat[missdat$timepoint=="d0",1]))[1],")"), ifelse(
    timepoint=="preld", paste0("Pre-Lymphodepletion (n=",dim(unique(missdat[(missdat$timepoint=="preld"),1]))[1],")"), paste0("Pre-Apheresis (n=",dim(unique(missdat[missdat$timepoint=="preaph",1]))[1],")")
  )
  )
  ) %>%
  mutate(poptot=ifelse(dx=="LBCL", dim(unique(missdat[(missdat$dx=="LBCL"),1]))[1],
                       ifelse(dx=="MCL",  dim(unique(missdat[(missdat$dx=="MCL"),1]))[1],
                              dim(unique(missdat[(missdat$dx=="FL"),1]))[1]
                       )
  )) %>%
  mutate(dx=ifelse(dx=="LBCL", paste0("LBCL (n=",dim(unique(missdat[(missdat$dx=="LBCL"),1]))[1],")"),
                   ifelse(dx=="MCL",  paste0("MCL (n=",dim(unique(missdat[(missdat$dx=="MCL"),1]))[1],")"),
                          paste0("FL (n=",dim(unique(missdat[(missdat$dx=="FL"),1]))[1],")")
                   )
  )) %>%
  left_join(numcohort, by="cohort") %>% mutate(poptot=n) %>% filter(cohort!="MSK Development")

missdat_dx <- missdat2 %>% filter(grepl("Infusion", timepoint)) %>%
  select(!c(record_id, timepoint)) %>% group_by(nacount, cohort, poptot) %>%
  summarize(n=n()/poptot) %>% ungroup() %>% select(!poptot) %>% distinct() %>% pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  pivot_longer(!nacount, names_to = "cohort", values_to = "n") %>%
  mutate(cohort=ifelse(cohort=="MSK Validation", "MSK LBCL", ifelse(cohort=="Center Validation", "Center (SMC+HMH) LBCL", "NHL (MCL+FL)")))

# Supp Figure 4a
sfig4a <- ggplot(missdat_dx, aes(x=nacount, y=n, fill=cohort)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  labs(x = "No. of Pre-Infusion Missing Labs per Patient",
       y = "Proportion of Cohort",
       #subtitle = "",
       #caption=paste0("Adj. R2: ", summary(lm(ptau~ftau, data=fig2_df1))$adj.r.squared, "\n","p < 0.0001")
  )+
  scale_x_continuous(breaks=0:8, labels=0:8, limits=c(-1,8.5))+
  guides(fill = guide_legend(title = "Cohort")) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, size=17),
    #element_text(angle = 90, vjust=-0.001, hjust=1,  size=15, face="bold"),
    axis.text.y = element_text(angle = 0, size=17),
    axis.title.x = element_text(angle = 0, size=17, face="bold"),
    axis.title.y = element_text(size=17, face="bold"),
    title = element_text(angle = 0, size=15, face="bold"),
    legend.text = element_text(angle = 0, size=17),
    legend.position = c(.8, 0.8)
  ) +
  scale_fill_manual(values=pal_npg("nrc")(9)[3:5])


# Evaluate probability of inflammatory cluster assignment to patients who have no missing laboratory data at the
# d0 and preld timepoints (in long format). No patients have complete laboratory data at the pre-apheresis timepoint.
# For the below analyses, we do not need to apply the filter analysis_type==1 filter (complete clinical metadata)
# - because we only care about evaluating concordance of cluster assignment using laboratory data here.
partclus1 <- bind_rows(
  reformat_table_tp(df_all_chrt1 %>% column_to_rownames("record_id") %>% select(starts_with("d0_")) %>%
                      select(contains(cluster_labs)) %>% na.omit() %>%
                      rename_all(~stringr::str_replace(.,paste0("^", "d0_"),"")) %>%
                      rename_with(~ paste0("d0_", .)) %>% rownames_to_column("record_id"), dev_mu, dev_sig, dev_pro, dev_inflammClus, tp="d0_") %>%
    group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
    mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% 
    mutate(timepoint="d0"),
  reformat_table_tp(df_all_chrt1 %>% column_to_rownames("record_id") %>% select(starts_with("preld_")) %>%
                      select(contains(cluster_labs)) %>% na.omit() %>%
                      rename_all(~stringr::str_replace(.,paste0("^", "preld_"),"")) %>%
                      rename_with(~ paste0("d0_", .)) %>% rownames_to_column("record_id"), dev_mu, dev_sig, dev_pro, dev_inflammClus, tp="d0_") %>%
    group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
    mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% 
    mutate(timepoint="preld")
) %>% mutate(record_id=paste0(record_id, "_temp_", 1:n())) %>% select(!cluster)

# For these patients, evaluate inflammatory cluster probability when only albumin, hgb, crp, ldh, ast, and alp are used.
partclus2 <-reformat_table_tp(partclus1 %>% rename(fulldat_tau=tau) %>%
                                mutate_at(vars(starts_with(focus_labs)), ~.*NA),
                              dev_mu, dev_sig, dev_pro, dev_inflammClus, tp="d0_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
  mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% 
  select(!cluster) %>% rename(cytomissdat_tau=tau) %>%
  mutate(fulldat_tau = round(fulldat_tau, 3), cytomissdat_tau=round(cytomissdat_tau, 3))

# Link with metadata
mis6_df <- partclus2 %>% select(record_id, cytomissdat_tau, fulldat_tau, timepoint) %>%
  mutate(record_id=gsub("(.*)_temp.*", "\\1", record_id)) %>% left_join(
    df_hc %>% select(record_id, dx=dx_simple), by="record_id"
  ) %>% mutate(part_cluster=ifelse(cytomissdat_tau>=0.5, "Inflammatory", "Non-Inflammatory"), full_cluster=ifelse(fulldat_tau>=0.5, "Inflammatory", "Non-Inflammatory")) %>%
  mutate(dx=ifelse(dx=="FL" | dx =="MCL", "MCL+FL", dx))

# For LBCL at d0, LBCL at preld, MCL+FL at d0, and MCL+FL at preld, evaluate concordance in cluster assignment as well as a linear regression between
# cluster assignment probabilities when all labs are available vs when only albumin, hgb, ldh, crp, ast, and alp are available.
mis6_df_all <- bind_rows(
  mis6_df %>% filter(dx=="LBCL") %>%
    select(!c(contains("tau"), record_id, dx)) %>%
    mutate(concord=ifelse(part_cluster==full_cluster, 1, 0)) %>% #group_by(timepoint) %>%
    summarize(value=sum(concord)/n(), n=n()) %>% mutate(metric="ccd"), #%>% mutate(feat=timepoint) %>% select(!timepoint),
  mis6_df %>% filter(timepoint=="d0") %>%
    select(!c(contains("tau"), record_id, timepoint)) %>%
    mutate(concord=ifelse(part_cluster==full_cluster, 1, 0)) %>% #group_by(dx) %>%
    summarize(value=sum(concord)/n(), n=n()) %>% mutate(metric="ccd"),# %>% mutate(feat=dx) %>% select(!dx),
  mis6_df %>% filter(dx=="LBCL") %>%
    select(!c(contains("cluster"), record_id, dx)) %>%
    #group_by(timepoint) %>%
    reframe(value=ccc_vec(estimate = cytomissdat_tau, truth=fulldat_tau), n=n()) %>% mutate(metric="ccc"), #%>% mutate(feat=timepoint) %>% select(!timepoint),
  mis6_df %>% filter(timepoint=="d0") %>%
    select(!c(contains("cluster"), record_id, timepoint)) %>%
    #group_by(dx) %>%
    summarize(value=ccc_vec(estimate = cytomissdat_tau, truth=fulldat_tau), n=n()) %>% mutate(metric="ccc"), #%>% mutate(feat=dx) #%>% select(!dx)
) 

library(ggrepel)

###### Cluster assignment with randomly missing lab values Supp Figure 4c ######
# Pick the same group of patients missing no laboratory data at MSK.
set.seed(NULL)
pclusterval_PARTdat_randit <- partclus1 %>% column_to_rownames("record_id") %>%
  select(starts_with("d0_")) %>%
  select(contains(cluster_labs))
# Initialize data.frame - probabilities when no data is missing.
pred_clus <- data.frame(record_id="id_test", rartp1=0.5, misdat=0, iter=0)
for(m in 1:dim(pclusterval_PARTdat_randit)[1]){
  pred_clus <- bind_rows(pred_clus, data.frame(
    record_id=rownames(pclusterval_PARTdat_randit)[m],
    rartp1=assign_clust(x=pclusterval_PARTdat_randit[m,], dev_mu, dev_sig, dev_pro)[dev_inflammClus],
    misdat=0,
    iter=1
  )
  )
}
pred_clus <- pred_clus[-1,]

# Evaluate probability of inflammatory cluster assignment -
# for 100 different iterations of n missing lab values per patients (where n=1 to 7)
m <- 0
for(l in 1:100) {
  old_time = Sys.time()
  for(i in 1:7){
    pclusterval_PARTdat_randit <- partclus1 %>% column_to_rownames("record_id") %>%
      select(starts_with("d0_")) %>%
      select(contains(cluster_labs))

    for(k in 1:dim(pclusterval_PARTdat_randit)[1]){
      m <- m+1
      set.seed(m)
      pclusterval_PARTdat_randit[k,c(sample(colnames(pclusterval_PARTdat_randit %>%
                                                       select(starts_with("d0_")) %>%
                                                       select(contains(cluster_labs))), i, replace=FALSE))] = rep(NA, i)
    }
    for(j in 1:dim(pclusterval_PARTdat_randit)[1]){
      pred_clus <- bind_rows(pred_clus, data.frame(
        record_id=rownames(pclusterval_PARTdat_randit)[j],
        rartp1=assign_clust(x=pclusterval_PARTdat_randit[j,], dev_mu, dev_sig, dev_pro)[dev_inflammClus],
        misdat=i,
        iter=l
      )
      )
    }
    print(i)
  }
  print(paste0("l=", l))
  print(Sys.time() - old_time)
}
pred_clus1 <- pred_clus %>% left_join(partclus1 %>% select(record_id, ftau=tau, timepoint), by="record_id")

# Link cluster predictions with with metadata
pred_clus2 <- pred_clus1 %>% mutate(record_id=gsub("(.*)_temp.*", "\\1", record_id)) %>%
  left_join(df_hc %>% select(record_id, dx=dx_simple), by="record_id")

# Calculate cluster assignment concordance and adjusted R2 from linear regression between cluster assignment -
# - probabilities with and without missing lab data.
df_tp <- data.frame(ccc=0, concord = 0,   timepoint="first", misdat=0, iter=0)
df_dx <- data.frame(ccc=0, concord = 0,  timepoint="temp", misdat=0, iter=0)
df_lbcld0 <- data.frame(ccc=0, concord = 0, pval=1, timepoint="ref", misdat=0, iter=0)
for(i in 1:7){
  old_time = Sys.time()
  tdf <- pred_clus2 %>% filter(misdat==i)
  for(j in 1:100){
    tdfit1 <- tdf %>% filter(iter==j)
    for(k in 1:length(unique(tdfit1$timepoint))){
      tdfit <- tdfit1 %>% filter(dx=="LBCL" & timepoint=="preld") %>%
        mutate(pclus=ifelse(rartp1>=0.5, 1,-1), fclus=ifelse(ftau>=0.5, 1, -1)) %>%
        mutate(ccd=ifelse(pclus==fclus, 1, 0))
      tfitfit_summ <- ccc(estimate=rartp1,  truth=ftau, data=tdfit)
      df_tp <- bind_rows(df_tp, data.frame(
        ccc = tfitfit_summ$.estimate,
        concord = sum(tdfit$ccd)/length(tdfit$ccd),
        timepoint = "Pre-Lymphodepletion\nLBCL",
        misdat=i,
        iter=j
      ))
    }
    for(k in 1:length(unique(tdfit1$dx))){
      tdfit <- tdfit1 %>% filter(timepoint=="d0"&dx!="LBCL") %>%
        mutate(pclus=ifelse(rartp1>=0.5, 1,-1), fclus=ifelse(ftau>=0.5, 1, -1)) %>%
        mutate(ccd=ifelse(pclus==fclus, 1, 0))
      tfitfit_summ <- ccc(estimate=rartp1,  truth=ftau, data=tdfit)
      df_dx <- bind_rows(df_dx, data.frame(
        ccc = tfitfit_summ$.estimate,
        concord = sum(tdfit$ccd)/length(tdfit$ccd),
        timepoint = "Pre-Infusion\nMCL+FL",
        misdat=i,
        iter=j
      ))
    }
    for(k in 1:1){
      tdfit <- tdfit1 %>% filter(timepoint=="d0"&dx=="LBCL")%>% #filter(dx==unique(tdfit1$dx)[1])%>%
        mutate(pclus=ifelse(rartp1>=0.5, 1,-1), fclus=ifelse(ftau>=0.5, 1, -1)) %>%
        mutate(ccd=ifelse(pclus==fclus, 1, 0))
      tfitfit_summ <- ccc(estimate=rartp1,  truth=ftau, data=tdfit)
      df_lbcld0 <- bind_rows(df_lbcld0, data.frame(
        ccc = tfitfit_summ$.estimate,
        concord = sum(tdfit$ccd)/length(tdfit$ccd),
        timepoint = "Pre-Infusion\nLBCL",
        misdat=i,
        iter=j
      ))
    }
    print(paste0("j=", j))
  }
  print(paste0("i=", i))
  print(Sys.time() - old_time)
}

# Link with metadata
partclus3 <- partclus1 %>% mutate(record_id=gsub("(.*)_temp.*", "\\1", record_id)) %>%
  left_join(df_hc %>% select(record_id, dx=dx_simple), by="record_id")

# Generate summary statistics for mean concordance and adj. R2 as well as +/- 2 standard deviation range.
partclus_df1 <- bind_rows(df_tp[-1,], df_dx[-1,], df_lbcld0[-1,]) %>%
  filter(!grepl("Lymphodep", timepoint)) %>% 
  group_by(misdat) %>%
  summarize(mean=mean(ccc), mean_ccd=mean(concord), low_ci=mean(ccc)-(2*sd(ccc)), high_ci=mean(ccc)+(2*sd(ccc)), low_ci_ccd=mean(concord)-(2*sd(concord)), high_ci_ccd=mean(concord)+(2*sd(concord))) %>%
  ungroup() %>%
  #mutate(timepoint1=ifelse(timepoint=="Pre-Lymphodepletion\nLBCL", paste0("Pre-Lymphodepletion\nLBCL (n=",dim(partclus3[(partclus3$timepoint=="preld"&partclus3$dx=="LBCL"),])[1],")"), ifelse(
  #  timepoint=="Pre-Infusion\nLBCL", paste0("Pre-Infusion\nLBCL (n=",dim(partclus3[(partclus3$timepoint=="d0"&partclus3$dx=="LBCL"),])[1],")"), ifelse(
  #    timepoint=="Pre-Infusion\nMCL+FL", paste0("Pre-Infusion\nMCL+FL (n=",dim(partclus3[((partclus3$dx %in% c("MCL", "FL"))&partclus3$timepoint=="d0"),])[1],")"), "Other"
  #  )))) %>%
  mutate(low_ci=ifelse(low_ci < 0, 0, low_ci), high_ci=ifelse(high_ci>1, 1, high_ci), low_ci_ccd=ifelse(low_ci_ccd < 0, 0, low_ci_ccd), high_ci_ccd=ifelse(high_ci_ccd>1, 1, high_ci_ccd))

partclus_df2 <- partclus_df1 %>% mutate(
  mean=mean_ccd,
  low_ci=low_ci_ccd,
  high_ci=high_ci_ccd
)

df1_tp <- partclus_df1 #%>% mutate(timepoint1=factor(timepoint1, levels=(c("Pre-Infusion\nLBCL (n=249)","Pre-Infusion\nMCL+FL (n=39)", "Pre-Lymphodepletion\nLBCL (n=26)" )))) #%>% filter(!grepl("Lymphodep", timepoint1)) %>% drop_levels()
df2_tp <- partclus_df2 #%>% mutate(timepoint1=factor(timepoint1, levels=(c("Pre-Infusion\nLBCL (n=249)","Pre-Infusion\nMCL+FL (n=39)", "Pre-Lymphodepletion\nLBCL (n=26)" )))) #%>% filter(!grepl("Lymphodep", timepoint1)) %>% drop_levels()

# Generate barchart of mean adjusted R2 values for each n randomly missing laboratory values (left half of figure 2b)
sfig4c_p1 <- ggplot(df1_tp, aes(x=misdat, y=mean)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=low_ci, ymax=high_ci), width=.2,
                position=position_dodge(.9)) +
  labs(x = "",
       y = "CCC",
  )+
  scale_x_continuous(breaks=1:8, limits=c(0.05,7.5))+
  guides(fill = guide_legend(title = ""), ncol=2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, size=17),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, size=17, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    title = element_text(angle = 0, size=15, face="bold"),
    legend.text = element_text(angle = 0, size=17),
    legend.position = "bottom"
  ) +
  scale_fill_npg() + coord_flip() +scale_y_reverse()

# Just to grab the legend
legend <- ggpubr::get_legend(sfig4c_p1)

# Re-generate barchart of mean adjusted R2 values for each n randomly missing laboratory values (left half of figure 2b)/
# This time without the legend.
sfig4c_p1 <- ggplot(df1_tp, aes(x=misdat, y=mean)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=low_ci, ymax=high_ci), width=.2,
                position=position_dodge(.9)) +
  labs(x = "",
       y = "CCC",
  )+
  scale_x_continuous(breaks=8:1, labels=8:1, transform = scales::transform_reverse())+
  scale_y_continuous(breaks=c(1, 0.75, 0.5, 0.25, 0), labels=c(1.00, 0.75, 0.50, 0.25, 0.00), transform = scales::transform_reverse())+
  guides(fill = guide_legend(title = "Timepoint"), ncol=2) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, size=17),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, size=17, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    title = element_text(angle = 0, size=15, face="bold"),
    legend.text = element_text(angle = 0, size=17),
    legend.position = "none"
  ) +
  scale_fill_npg() + coord_flip() 

# Generate barchart of mean concordance values for each n randomly missing laboratory values (right half of figure 2b)
sfig4c_p2 <- ggplot(df2_tp, aes(x=misdat, y=mean)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=low_ci, ymax=high_ci), width=.2,
                position=position_dodge(.9)) +
  labs(x = "",
       y = "Hard Label Concordance",
  )+
  scale_x_continuous(breaks=8:1, labels=8:1, transform = scales::transform_reverse())+
  guides(fill = guide_legend(title = "Timepoint")) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, size=17),
    axis.text.y = element_text(angle = 0, size=17, face="bold"),
    axis.title.x = element_text(angle = 0, size=17, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    title = element_text(angle = 0, size=15, face="bold"),
    legend.text = element_text(angle = 0, size=17),
    legend.position = "none"
  ) +
  scale_fill_npg()+ coord_flip()

# Figure supp4c 
sfig4c <- plot_grid(plot_grid(sfig4c_p1, NULL, sfig4c_p2, labels=c("", "", ""), nrow=1, rel_widths=c(1, 0, 1),label_size = 30, scale=1),
                   legend, labels=c("", ""), ncol=1, rel_heights = c(1, 0.1), scale=0.9, label_size=30
)

