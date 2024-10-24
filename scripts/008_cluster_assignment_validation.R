rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
library(tidyverse)
library(ggsci)

source("scripts/000_functions_constants.R")

## Generate unit vectors of inflammatory and non-inflammatory parameter-mean vectors
ref_inflam_vec <- dev_mu[,which.max(dev_mu["d0_il6",])]
unit_ref_inflam_vec <- ref_inflam_vec / sqrt(sum(ref_inflam_vec^2))

ref_LEASTinflam_vec <- dev_mu[,which.min(dev_mu["d0_il6",])]
unit_ref_LEASTinflam_vec <- ref_LEASTinflam_vec / sqrt(sum(ref_LEASTinflam_vec^2))

######## 
coh_df1 <- df_labs_all %>% select(record_id, cohort, starts_with("noLCS_d0_"), starts_with("uln_d0_")) %>% 
  filter(cohort=="MSK Development") %>% select(!cohort) %>% 
  pivot_longer(!c(record_id), 
               names_to=c(".value", "lab"),
               names_sep="_d0_"
  ) %>% 
  mutate(truln=noLCS/uln) %>%
  mutate(truln=ifelse(lab %in% c("crp", "ferritin"), log10(truln), truln)) %>% 
  mutate(lab=paste0("d0_", lab))  %>% select(record_id, lab, truln) 

# Get InflaMix assignment for cohort. 
coh1_metrics <- coh_df1 %>%
  group_by(lab) %>% dplyr::reframe(truln_mean=mean(truln), truln_sd=sd(truln)) %>% ungroup()

coh1_preim <- coh_df1 %>% left_join(coh1_metrics, by="lab") %>% 
  mutate(scld=(truln-truln_mean)/truln_sd) %>%
  select(record_id, lab, scld) %>% 
  pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
  select(record_id, contains(cluster_labs))

d0_labs <- df_labs_all %>% 
  select(starts_with("d0"))

misrate <- as.data.frame(t(colSums(is.na(d0_labs))/dim(d0_labs)[1]))

coh1_df_variants <- list(alldat=coh1_preim)
for(z in 1:10){
  set.seed(z)
  rmat <- as.data.frame(matrix(runif(dim(coh1_preim)[1]*6), dim(coh1_preim)[1]))
  colnames(rmat) <- c("lfts", "cytokines", "crp", "ldh", "ddimer", "ferritin")
  rmat$lfts[rmat$lfts < misrate$d0_ast] <- NA
  rmat$lfts[rmat$lfts > misrate$d0_ast] <- 1
  
  rmat$cytokines[rmat$cytokines < misrate$d0_il6] <- NA
  rmat$cytokines[rmat$cytokines > misrate$d0_il6] <- 1
  
  rmat$cytokines[rmat$cytokines < misrate$d0_il6] <- NA
  rmat$cytokines[rmat$cytokines > misrate$d0_il6] <- 1
  
  rmat$crp[rmat$crp < misrate$d0_crp] <- NA
  rmat$crp[rmat$crp > misrate$d0_crp] <- 1
  
  rmat$ldh[rmat$ldh < misrate$d0_ldh] <- NA
  rmat$ldh[rmat$ldh > misrate$d0_ldh] <- 1
  
  rmat$ddimer[rmat$ddimer < misrate$d0_ddimer] <- NA
  rmat$ddimer[rmat$ddimer > misrate$d0_ddimer] <- 1
  
  rmat$ferritin[rmat$ferritin < misrate$d0_ferritin] <- NA
  rmat$ferritin[rmat$ferritin > misrate$d0_ferritin] <- 1
  
  dfcalv <- data.frame(
    record_id=coh1_preim$record_id,
    d0_albumin=rmat$lfts * coh1_preim$d0_albumin,
    d0_alk=rmat$lfts * coh1_preim$d0_alk,
    d0_ast=rmat$lfts * coh1_preim$d0_ast,
    d0_ferritin=rmat$ferritin * coh1_preim$d0_ferritin,
    d0_hb=coh1_preim$d0_hb,
    d0_ldh=rmat$ldh * coh1_preim$d0_ldh,
    d0_plt=coh1_preim$d0_plt,
    d0_tbr=rmat$lfts * coh1_preim$d0_tbr,
    d0_il10=rmat$cytokines * coh1_preim$d0_il10,
    d0_il6=rmat$cytokines * coh1_preim$d0_il6,
    d0_tnfa=rmat$cytokines * coh1_preim$d0_tnfa,
    d0_crp=rmat$crp * coh1_preim$d0_crp,
    d0_ddimer=rmat$ddimer * coh1_preim$d0_ddimer,
    d0_wbc=coh1_preim$d0_wbc
  ) 
  coh1_df_variants[[z+1]] <- dfcalv
}

coh1_df_variants[[12]] <- coh1_preim %>% mutate_at(vars(starts_with(focus_labs)), ~.*NA) 
coh1_preim <- coh1_preim %>% column_to_rownames("record_id")

allits <- data.frame(record_id="NA", imi_tau=0, im_tau=0, z=0, i=0)
k <- 1
for(z in 1:12){
  dfcalv <- coh1_df_variants[[z]]
  
  inflamix_assign_full <- reformat_table_tp(df=dfcalv,
                                            gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig,
                                            inflammClus = dev_inflammClus,
                                            least_inflammClus = dev_LEASTinflammClus,
                                            tp="d0_") %>%
    select(record_id, cluster, tau) %>% 
    #filter(cluster=="Inflammatory") %>% select(record_id, im_tau=tau)
    group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
    mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% select(record_id, im_tau=tau)
  
  i <- 1
  n = 100
  concord <- data.frame(ccd=1:n)
  lmcof <- data.frame(r2=1:n)
  oob_all <- data.frame(record_id="NA", imi_tau=0, im_tau=0, i=0)
  
  while(i <= n){
    set.seed(k)
    
    # Resampling with replacement
    resamples <- coh_df1 %>% pivot_wider(names_from = "lab", values_from = "truln", values_fill = NA)
    resamples <- resamples[sample(1:dim(resamples)[1], size=dim(resamples)[1], replace = TRUE),] %>%
      mutate(id=1:n()) %>%
      pivot_longer(!c(record_id, id), names_to = "lab", values_to = "truln") 
    
    resample_metrics <- resamples %>%
      group_by(lab) %>% dplyr::reframe(truln_mean=mean(truln), truln_sd=sd(truln)) %>% ungroup()
    
    resamples <- resamples %>% left_join(resample_metrics, by="lab") %>% 
      mutate(scld=(truln-truln_mean)/truln_sd) %>%
      select(record_id, id, lab, scld) %>% 
      pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
      column_to_rownames("id") %>% select(!record_id) %>%
      select(contains(cluster_labs))
    
    inflamixi <- Mclust(resamples, modelNames = "VVV", G = 2, verbose = FALSE)
    
    if(is.null(inflamixi)){
      k <- k+1
      next
    }

    # When we generate new clusters for each bootstrap, the cluster labels are not consistently carried over.
    # Identify which cluster is inflammatory by finding the cluster with a unit vector of mean laboratory values that is most parallel -
    # by calculating the dot product between the two - with the reference unit vector of our original model.
    dotprods <-  as.data.frame(inflamixi$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
      rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
      dplyr::group_by(cluster) %>% dplyr::summarize(dotprod = unit_ref_inflam_vec %*% meanvalues) %>% ungroup()
    LEASTdotprods <-  as.data.frame(inflamixi$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
      rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
      group_by(cluster) %>% summarize(dotprod = unit_ref_LEASTinflam_vec %*% meanvalues)
    inflamm_index <- which.max(as.double(unlist(dotprods[,2])))
    LEASTinflamm_index <- which.max(as.double(unlist(LEASTdotprods[,2])))
    
    coh1_preimi <- coh_df1 %>% left_join(resample_metrics, by="lab") %>% 
      mutate(scld=(truln-truln_mean)/truln_sd) %>%
      select(record_id, lab, scld) %>% 
      pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
      select(record_id, contains(cluster_labs))
    
    oob_imi <- reformat_table_tp(df=coh1_preimi ,
                                 gmm_mu=inflamixi$parameters$mean, 
                                 gmm_pro=inflamixi$parameters$pro, 
                                 gmm_sig=inflamixi$parameters$variance$sigma, 
                                 inflammClus = inflamm_index, 
                                 least_inflammClus = LEASTinflamm_index, 
                                 tp="d0_") %>%
      select(record_id, cluster, tau) %>% 
      group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% select(record_id, imi_tau=tau) %>%
      left_join(inflamix_assign_full, by="record_id")
    
    oob_all <- bind_rows(oob_all, oob_imi %>% mutate(i=i))
    
    i <- i+1
    k <- k+1
    print(paste0(z, "_", i))
  }
  allits <- bind_rows(allits, oob_all[-1,] %>% mutate(z=z)) 
}

oob_all1 <- allits[-1,] 

mypal <- pal_npg("nrc", alpha = 0.7)(9)[c(3,9,8)]

oob_all1 %>% 
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  group_by(mv_min, record_id) %>% summarize(im_tau=mean(im_tau), imi_tau=mean(imi_tau)) %>%
  #filter(mv_min=="No missing") %>%
  #filter(mv_min=="1-7 missing labs") %>%
  filter(mv_min=="8 missing labs (focus set)") %>%
  ggplot(aes(x=im_tau, y=imi_tau, color=mv_min, fill=mv_min)) + #, color=factor(mv_min))) + 
  xlim(0, 1)+
  theme_bw() +
  theme(
    text=element_text(size=22)
  ) + 
  geom_point(alpha=0.3, size=8)+
  stat_smooth(method = "lm", se = TRUE, fullrange = T, alpha=0.4) +
  geom_abline(slope=1, intercept=0, size=1.5) + 
  annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1, fill= "lightgray", alpha=0.4)  + 
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5 , fill= "lightgray", alpha=0.4) + 
  #annotate("rect", xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1, fill= "yellow") + 
  #annotate("rect", xmin = 0.5, xmax = 1, ymin = 0, ymax = 0.5, fill= "green") +
  xlab("") + ylab("") +
  scale_color_manual(values=mypal[3])+
  scale_fill_manual(values=mypal[3]) + guides(color="none", fill="none")


cccpre_val <- oob_all1 %>% 
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  group_by(mv_min, record_id) %>% summarize(im_tau=mean(im_tau), imi_tau=mean(imi_tau)) %>%
  filter(mv_min=="8 missing labs (focus set)") %>%
  #filter(mv_min=="1-7 missing labs") %>%
  #filter(mv_min=="No missing") %>%
  mutate(im_tau=im_tau)

sqrt(summary(lm(imi_tau~im_tau, data=cccpre_val))$r.squared)
cor(cccpre_val$imi_tau, cccpre_val$im_tau, method="pearson")
ccc(data=cccpre_val, truth=imi_tau, estimate=im_tau)

dev_boxplots <- oob_all1 %>%
  mutate_at(c("imi_tau", "im_tau"), function(x) {ifelse(x>=0.5, 1, 0)}) %>%
  mutate(ccd=im_tau+imi_tau) %>% mutate(ccd=ifelse(ccd==1, 0, 1)) %>%
  group_by(z, i) %>% mutate(ccds=sum(ccd), n=n()) %>% summarize(ccd=ccds/n) %>% ungroup() %>% distinct() %>%
  mutate(mv_min = ifelse(z==12, "8 missing labs (focus set)", ifelse(z==1, "No missing", "1-7 missing labs"))) %>%
  mutate(mv_min=factor(mv_min, levels=c("No missing", "1-7 missing labs", "8 missing labs (focus set)"))) %>% 
  ggplot(aes(x=mv_min, y=ccd, color=mv_min, fill=mv_min)) + 
  geom_boxplot(alpha=0.2) + ylab("Accuracy") + xlab("") +
  scale_color_manual(values=mypal)+
  scale_fill_manual(values=mypal) +theme_bw() + 
  theme(
    text=element_text(size=22),
    legend.position = "none"
  ) 

test <- oob_all1 %>%
  mutate_at(c("imi_tau", "im_tau"), function(x) {ifelse(x>=0.5, 1, 0)}) %>%
  mutate(ccd=im_tau+imi_tau) %>% mutate(ccd=ifelse(ccd==1, 0, 1)) %>%
  group_by(z, i) %>% mutate(ccds=sum(ccd), n=n()) %>% summarize(ccd=ccds/n) %>% ungroup() %>% distinct() %>%
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  mutate(mv_min=factor(mv_min, levels=c("No missing", "1-7 missing labs", "8 missing labs (focus set)"))) %>%
  filter(mv_min=="8 missing labs (focus set)") 

median(test$ccd)

####
######## 
coh_df2 <- df_labs_all %>% select(record_id, cohort, starts_with("noLCS_d0_"), starts_with("uln_d0_")) %>% 
  filter(cohort != "MSK Development") %>% select(!cohort) %>% 
  pivot_longer(!c(record_id), 
               names_to=c(".value", "lab"),
               names_sep="_d0_"
  ) %>% 
  mutate(truln=noLCS/uln) %>%
  mutate(truln=ifelse(lab %in% c("crp", "ferritin"), log10(truln), truln)) %>% 
  mutate(lab=paste0("d0_", lab))  %>% select(record_id, lab, truln)  %>%
  pivot_wider(names_from = "lab", values_from = "truln", values_fill = NA) %>%
  drop_na() %>%
  pivot_longer(!record_id, names_to = "lab", values_to = "truln")

# Get InflaMix assignment for cohort. 
coh2_preim <- coh_df2 %>% left_join(coh1_metrics, by="lab") %>% 
  mutate(scld=(truln-truln_mean)/truln_sd) %>%
  select(record_id, lab, scld) %>% 
  pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
  select(record_id, contains(cluster_labs)) 

coh2_df_variants <- list(alldat=coh2_preim)
for(z in 1:10){
  set.seed(z)
  rmat <- as.data.frame(matrix(runif(dim(coh2_preim)[1]*6), dim(coh2_preim)[1]))
  colnames(rmat) <- c("lfts", "cytokines", "crp", "ldh", "ddimer", "ferritin")
  rmat$lfts[rmat$lfts < misrate$d0_ast] <- NA
  rmat$lfts[rmat$lfts > misrate$d0_ast] <- 1
  
  rmat$cytokines[rmat$cytokines < misrate$d0_il6] <- NA
  rmat$cytokines[rmat$cytokines > misrate$d0_il6] <- 1
  
  rmat$cytokines[rmat$cytokines < misrate$d0_il6] <- NA
  rmat$cytokines[rmat$cytokines > misrate$d0_il6] <- 1
  
  rmat$crp[rmat$crp < misrate$d0_crp] <- NA
  rmat$crp[rmat$crp > misrate$d0_crp] <- 1
  
  rmat$ldh[rmat$ldh < misrate$d0_ldh] <- NA
  rmat$ldh[rmat$ldh > misrate$d0_ldh] <- 1
  
  rmat$ddimer[rmat$ddimer < misrate$d0_ddimer] <- NA
  rmat$ddimer[rmat$ddimer > misrate$d0_ddimer] <- 1
  
  rmat$ferritin[rmat$ferritin < misrate$d0_ferritin] <- NA
  rmat$ferritin[rmat$ferritin > misrate$d0_ferritin] <- 1
  
  dfcalv <- data.frame(
    record_id=coh2_preim$record_id,
    d0_albumin=rmat$lfts * coh2_preim$d0_albumin,
    d0_alk=rmat$lfts * coh2_preim$d0_alk,
    d0_ast=rmat$lfts * coh2_preim$d0_ast,
    d0_ferritin=rmat$ferritin * coh2_preim$d0_ferritin,
    d0_hb=coh2_preim$d0_hb,
    d0_ldh=rmat$ldh * coh2_preim$d0_ldh,
    d0_plt=coh2_preim$d0_plt,
    d0_tbr=rmat$lfts * coh2_preim$d0_tbr,
    d0_il10=rmat$cytokines * coh2_preim$d0_il10,
    d0_il6=rmat$cytokines * coh2_preim$d0_il6,
    d0_tnfa=rmat$cytokines * coh2_preim$d0_tnfa,
    d0_crp=rmat$crp * coh2_preim$d0_crp,
    d0_ddimer=rmat$ddimer * coh2_preim$d0_ddimer,
    d0_wbc=coh2_preim$d0_wbc
  ) 
  coh2_df_variants[[z+1]] <- dfcalv
}

coh2_df_variants[[12]] <- coh2_preim %>% mutate_at(vars(starts_with(focus_labs)), ~.*NA) 
coh2_preim <- coh2_preim %>% column_to_rownames("record_id")

allits <- data.frame(record_id="NA", im2_tau=0, im_tau=0, z=0, i=0)
k <- 1
for(z in 1:12){
  dfcalv <- coh2_df_variants[[z]]
  
  inflamix_assign_full <- reformat_table_tp(df=dfcalv,
                                            gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig,
                                            inflammClus = dev_inflammClus,
                                            least_inflammClus = dev_LEASTinflammClus,
                                            tp="d0_") %>%
    select(record_id, cluster, tau) %>% 
    #filter(cluster=="Inflammatory") %>% select(record_id, im_tau=tau)
    group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
    mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% select(record_id, im_tau=tau)
  
  i <- 1
  n = 100
  oob_all <- data.frame(record_id="NA", im2_tau=0, im_tau=0, i=0)
  
  while(i <= n){
    set.seed(k)
    
    # Resampling with replacement
    resamples <- coh_df2 %>% pivot_wider(names_from = "lab", values_from = "truln", values_fill = NA)
    resamples <- resamples[sample(1:dim(resamples)[1], size=dim(resamples)[1], replace = TRUE),] %>%
      mutate(id=1:n()) %>%
      pivot_longer(!c(record_id, id), names_to = "lab", values_to = "truln") 
    
    resample_metrics <- resamples %>%
      group_by(lab) %>% dplyr::reframe(truln_mean=mean(truln), truln_sd=sd(truln)) %>% ungroup()
    
    resamples <- resamples %>% left_join(resample_metrics, by="lab") %>% 
      mutate(scld=(truln-truln_mean)/truln_sd) %>%
      select(record_id, id, lab, scld) %>% 
      pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
      column_to_rownames("id") %>% select(!record_id) %>%
      select(contains(cluster_labs))
    
    inflamixi <- Mclust(resamples, modelNames = "VVV", G = 2, verbose = FALSE)
    
    if(is.null(inflamixi)){
      k <- k+1
      next
    }
    
    # When we generate new clusters for each bootstrap, the cluster labels are not consistently carried over.
    # Identify which cluster is inflammatory by finding the cluster with a unit vector of mean laboratory values that is most parallel -
    # by calculating the dot product between the two - with the reference unit vector of our original model.
    dotprods <-  as.data.frame(inflamixi$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
      rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
      dplyr::group_by(cluster) %>% dplyr::summarize(dotprod = unit_ref_inflam_vec %*% meanvalues) %>% ungroup()
    LEASTdotprods <-  as.data.frame(inflamixi$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
      rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
      group_by(cluster) %>% summarize(dotprod = unit_ref_LEASTinflam_vec %*% meanvalues)
    inflamm_index <- which.max(as.double(unlist(dotprods[,2])))
    LEASTinflamm_index <- which.max(as.double(unlist(LEASTdotprods[,2])))
    
    coh2_preimi <- coh_df2 %>% left_join(resample_metrics, by="lab") %>% 
      mutate(scld=(truln-truln_mean)/truln_sd) %>%
      select(record_id, lab, scld) %>% 
      pivot_wider(names_from = "lab", values_from = "scld", values_fill = NA) %>%
      select(record_id, contains(cluster_labs))
    
    oob_imi <- reformat_table_tp(df=coh2_preimi ,
                                 gmm_mu=inflamixi$parameters$mean, 
                                 gmm_pro=inflamixi$parameters$pro, 
                                 gmm_sig=inflamixi$parameters$variance$sigma, 
                                 inflammClus = inflamm_index, 
                                 least_inflammClus = LEASTinflamm_index, 
                                 tp="d0_") %>%
      select(record_id, cluster, tau) %>% 
      group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% select(record_id, im2_tau=tau) %>%
      left_join(inflamix_assign_full, by="record_id")
    
    oob_all <- bind_rows(oob_all, oob_imi %>% mutate(i=i))
    
    i <- i+1
    k <- k+1
    print(paste0(z, "_", i))
  }
  allits <- bind_rows(allits, oob_all[-1,] %>% mutate(z=z)) 
}

oob_all1 <- allits[-1,] 

mypal <- pal_npg("nrc", alpha = 0.7)(9)[c(3,9,8)]

oob_all1 %>% 
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  group_by(mv_min, record_id) %>% summarize(im_tau=mean(im_tau), im2_tau=mean(im2_tau)) %>%
  #filter(mv_min=="No missing") %>%
  #filter(mv_min=="1-7 missing labs") %>%
  filter(mv_min=="8 missing labs (focus set)") %>%
  ggplot(aes(x=im_tau, y=im2_tau, color=mv_min, fill=mv_min)) + #, color=factor(mv_min))) + 
  xlim(0, 1)+
  theme_bw() +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.00), limits=c(0, 1.0)) + 
  geom_point(alpha=0.3, size=8)+
  stat_smooth(method = "lm", se = TRUE, fullrange = T, alpha=0.4) +
  geom_abline(slope=1, intercept=0, size=1.5) + 
  annotate("rect", xmin = 0.5, xmax = 1, ymin = 0.5, ymax = 1, fill= "lightgray", alpha=0.4)  + 
  annotate("rect", xmin = 0, xmax = 0.5, ymin = 0, ymax = 0.5 , fill= "lightgray", alpha=0.4) + 
  theme(
    text=element_text(size=22)
  ) + 
  #annotate("rect", xmin = 0, xmax = 0.5, ymin = 0.5, ymax = 1, fill= "yellow") + 
  #annotate("rect", xmin = 0.5, xmax = 1, ymin = 0, ymax = 0.5, fill= "green") +
  xlab("") + ylab("") +
  scale_color_manual(values=mypal[3])+
  scale_fill_manual(values=mypal[3]) + guides(color="none", fill="none")


cccpre_val <- oob_all1 %>% 
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  group_by(mv_min, record_id) %>% summarize(im_tau=mean(im_tau), im2_tau=mean(im2_tau)) %>%
  filter(mv_min=="8 missing labs (focus set)") %>%
  #filter(mv_min=="1-7 missing labs") %>%
  #filter(mv_min=="No missing") %>%
  mutate(im_tau=im_tau)

sqrt(summary(lm(im2_tau~im_tau, data=cccpre_val))$r.squared)
cor(cccpre_val$im2_tau, cccpre_val$im_tau, method="pearson")
ccc(data=cccpre_val, truth=im2_tau, estimate=im_tau)

oob_all1 %>%
  mutate_at(c("im2_tau", "im_tau"), function(x) {ifelse(x>=0.5, 1, 0)}) %>%
  mutate(ccd=im_tau+im2_tau) %>% mutate(ccd=ifelse(ccd==1, 0, 1)) %>%
  group_by(z, i) %>% mutate(ccds=sum(ccd), n=n()) %>% summarize(ccd=ccds/n) %>% ungroup() %>% distinct() %>%
  mutate(mv_min = ifelse(z==12, "8 missing labs (focus set)", ifelse(z==1, "No missing", "1-7 missing labs"))) %>%
  mutate(mv_min=factor(mv_min, levels=c("No missing", "1-7 missing labs", "8 missing labs (focus set)"))) %>% 
  ggplot(aes(x=mv_min, y=ccd, color=mv_min, fill=mv_min)) + 
  geom_boxplot(alpha=0.2) + ylab("Accuracy") + xlab("") +
  scale_color_manual(values=mypal)+
  scale_fill_manual(values=mypal) +theme_bw()

test <- oob_all1 %>%
  mutate_at(c("im2_tau", "im_tau"), function(x) {ifelse(x>=0.5, 1, 0)}) %>%
  mutate(ccd=im_tau+im2_tau) %>% mutate(ccd=ifelse(ccd==1, 0, 1)) %>%
  group_by(z, i) %>% mutate(ccds=sum(ccd), n=n()) %>% summarize(ccd=ccds/n) %>% ungroup() %>% distinct() %>%
  mutate(mv_min = ifelse(z==1, "No missing", ifelse(z==12, "8 missing labs (focus set)", "1-7 missing labs"))) %>%
  mutate(mv_min=factor(mv_min, levels=c("No missing", "1-7 missing labs", "8 missing labs (focus set)"))) %>%
  filter(mv_min=="8 missing labs (focus set)") 

median(test$ccd)

















