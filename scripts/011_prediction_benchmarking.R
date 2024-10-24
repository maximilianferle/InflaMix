rm(list=ls())
library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(ggsurvfit)
library(ggpubr)
library(broom)
library(pROC)
library(glmnet)
library(riskRegression)
library(compareC)
library(ggsci)
library(dcurves)
library(caret)
library(rms)
library(yardstick)
library(MESS)
library(cowplot)
library(data.table)

load("output/model.RData")
load("output/scaledata.RData")

inflamix <- d0clusters
im_mu=dev_mu
im_sig=dev_sig
im_pro=dev_pro
ind_inflclus=dev_inflammClus
ind_ninflclus=dev_LEASTinflammClus
source("scripts/000_functions_constants.R")

# No cytokine mixture modeling --------------------------------------------
cluster_labs_nocytos <- cluster_labs[!(cluster_labs %in% c("il6", "il10", "tnfa"))]

reformat_table_tp_nocytos <- function(df, gmm_mu, gmm_sig, gmm_pro, inflammClus, least_inflammClus, tp){
  df_2clust_d0labs <- df %>% select(record_id, starts_with(tp)) %>%
    select(record_id, contains(cluster_labs_nocytos)) %>% column_to_rownames("record_id") %>%
    rename_all(~stringr::str_replace(.,paste0("^", tp),"")) %>%
    rename_with(~ paste0("d0_", .))
  reformat.dat <- data.frame(cluster=1, tau=0.5, record_id="id_test")
  for(k in 1:dim(df)[1]){
    if(num_na <- sum(is.na(df_2clust_d0labs[k,]))==dim(df_2clust_d0labs)[2]){
      reformat.dat <- reformat.dat %>% bind_rows(data.frame(
        cluster = 1:length(gmm_pro),
        tau = NA,
        record_id=rownames(df_2clust_d0labs)[k]
      ))
    } else(
      reformat.dat <- reformat.dat %>% bind_rows(data.frame(
        cluster = 1:length(gmm_pro),
        tau = assign_clust(x=df_2clust_d0labs[k,], gmm_mu, gmm_sig, gmm_pro),
        record_id=rownames(df_2clust_d0labs)[k]
      ))
    )
  }
  if(length(gmm_pro)>2){
    clus_recode <- reformat.dat %>% select(cluster) %>% filter(!(cluster %in% c(inflammClus, least_inflammClus))) %>%
      distinct() %>% mutate(newclus=1:n())
    reformat.dat <- reformat.dat[-1,] %>% left_join(clus_recode, by="cluster") %>%
      inner_join(df, by="record_id") %>%
      mutate(cluster=factor(ifelse(cluster==inflammClus, "Inflammatory",
                                   ifelse(cluster==least_inflammClus,
                                          "Non-Inflammatory",
                                          ifelse(as.integer(newclus)==1, "Neutral Cluster",
                                                 paste0("Neutral Cluster-", as.integer(newclus) ))))))%>%
      mutate(cluster=factor(cluster, levels=rev(levels(cluster)))) %>%
      mutate(cluster=relevel(cluster, "Non-Inflammatory")) %>%
      filter(tau!=0) %>% filter(!is.nan(tau))
  } else{
    reformat.dat <- reformat.dat[-1,] %>%
      inner_join(df, by="record_id") %>%
      mutate(cluster=factor(ifelse(cluster==inflammClus, "Inflammatory","Non-Inflammatory"))) %>%
      mutate(cluster=factor(cluster, levels=rev(levels(cluster)))) %>%
      mutate(cluster=relevel(cluster, "Non-Inflammatory")) %>%
      filter(tau!=0) %>% filter(!is.nan(tau))
  }
  return(reformat.dat)
}
df_dev_nocytos <-  left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>% # Must have curated metadata
  filter(cohort=="MSK Development") %>%
  select(record_id, starts_with("d0_")) %>%
  select(record_id, contains(cluster_labs_nocytos)) %>%
  column_to_rownames("record_id")

# Define model parameters
d0clusters_nocytos <- Mclust(df_dev_nocytos %>% select(starts_with("d0_")), modelNames = "VVV", G = 2)

# Define the mean laboratory vector for the inflammatory cluster (the one that has the higher IL-6 value)
ref_inflam_vec_nocytos <- d0clusters_nocytos$parameters$mean[,which.max(d0clusters_nocytos$parameters$mean["d0_crp",])]
unit_ref_inflam_vec_nocytos <- ref_inflam_vec_nocytos / sqrt(sum(ref_inflam_vec_nocytos^2))

# Define the model paramaters
dev_mu_nocytos <- d0clusters_nocytos$parameters$mean
dev_sig_nocytos <- d0clusters_nocytos$parameters$variance$sigma
dev_pro_nocytos <- d0clusters_nocytos$parameters$pro
dev_inflammClus_nocytos <- which.max(d0clusters_nocytos$parameters$mean["d0_crp",])
dev_LEASTinflammClus_nocytos <- which.min(d0clusters_nocytos$parameters$mean["d0_crp",])

# Functions ---------------------------------------------------------------
coxsplity=function(y, nfolds, method, prop){
  N=nrow(y)
  tem=data.frame(y, i=seq(N), foldid=0)
  tem=tem[order(y[, "time"], y[, "status"]), ]
  n1=sum(y[, "status"]);n2=N-n1
  
  if(!missing(prop)) {invprop=round(1/(1-prop)); tip=floor(prop*invprop); vip=invprop-tip}
  
  if(method=="crossfold") {
    tem$foldid[tem[, "status"]==1]=sample(rep(seq(nfolds), length=n1))
    tem$foldid[tem[, "status"]==0]=sample(rep(seq(nfolds), length=n2))
  }
  
  if(method=="simple") {
    tem$foldid[tem[, "status"]==1]=sample( c(rep(1, length=floor(n1*prop)), rep(2, length=ceiling(n1-(n1*prop))) ) ) 
    tem$foldid[tem[, "status"]==0]=sample( c(rep(1, length=floor(n2*prop)), rep(2, length=ceiling(n2-(n2*prop))) ) )
  }
  
  if(method=="bootstrap") {
    tem$foldid[tem[, "status"]==1]=sample( c(rep(1, length=floor(n1*prop)), rep(2, length=ceiling(n1-(n1*prop))) ), replace = TRUE) 
    tem$foldid[tem[, "status"]==0]=sample( c(rep(1, length=floor(n2*prop)), rep(2, length=ceiling(n2-(n2*prop))) ), replace = TRUE )
  }
  
  foldid=tem$foldid[order(tem$i)]
  return(foldid)
}
collect_riskReg <- function(score, pobject, analyses) {
  pobjtable <- pobject$control$addtable2plot$table
  nanum=length(score$Calibration$plotframe$risk[score$Calibration$plotframe$model==analyses[1]]) - length(pobject$plotFrames[[1]]$Pred)
  namat=as.data.frame(matrix(data=rep(NA, nanum*2), nrow=nanum, ncol=2))
  colnames(namat) <- c("Pred", "Obs")
  i=2
  preploplot <- bind_rows(pobject$plotFrames[[analyses[1]]], namat) %>% mutate(analysis=analyses[1])
  while(i<=length(analyses)){
    nanum=length(score$Calibration$plotframe$risk[score$Calibration$plotframe$model==analyses[i]]) - length(pobject$plotFrames[[i]]$Pred)
    namat=as.data.frame(matrix(data=rep(NA, nanum*2), nrow=nanum, ncol=2))
    colnames(namat) <- c("Pred", "Obs")
    preploplot <- bind_rows(preploplot, bind_rows(pobject$plotFrames[[analyses[i]]], namat) %>% mutate(analysis=analyses[i]))
    i=i+1
  }
  if(score$response.type=="survival") {
    preploplot <- preploplot %>% mutate(id=score$Calibration$plotframe$ID, rug=score$Calibration$plotframe$risk, time=score$Calibration$plotframe$time, status=score$Calibration$plotframe$status, pseudovalue=score$Calibration$plotframe$pseudovalue)
  } else {
    preploplot <- preploplot %>% mutate(id=score$Calibration$plotframe$ID, rug=score$Calibration$plotframe$risk, time=NA, status=score$Calibration$plotframe$ReSpOnSe, pseudovalue=score$Calibration$plotframe$ReSpOnSe)
  }
  return(preploplot)
}
fit_riskReg_spline <- function(score, pobject, analyses, knots, lintype="spline") {
  pobjtable <- pobject$control$addtable2plot$table
  nanum=length(score$Calibration$plotframe$risk[score$Calibration$plotframe$model==analyses[1]]) - length(pobject$plotFrames[[1]]$Pred)
  namat=as.data.frame(matrix(data=rep(NA, nanum*2), nrow=nanum, ncol=2))
  colnames(namat) <- c("Pred", "Obs")
  i=2
  preploplot <- bind_rows(pobject$plotFrames[[analyses[1]]], namat) %>% mutate(analysis=analyses[1])
  while(i<=length(analyses)){
    nanum=length(score$Calibration$plotframe$risk[score$Calibration$plotframe$model==analyses[i]]) - length(pobject$plotFrames[[i]]$Pred)
    namat=as.data.frame(matrix(data=rep(NA, nanum*2), nrow=nanum, ncol=2))
    colnames(namat) <- c("Pred", "Obs")
    preploplot <- bind_rows(preploplot, bind_rows(pobject$plotFrames[[analyses[i]]], namat) %>% mutate(analysis=analyses[i]))
    i=i+1
  }
  preploplot <- preploplot %>% mutate(rug=score$Calibration$plotframe$risk) %>% mutate(analysis=factor(analysis, levels=model.names))
  
  if(lintype=="spline"){
    ploplot <- preploplot %>% ggplot() +#ggplot(aes(x=Pred, y=Obs, color=analysis)) + 
      geom_abline(slope=1, intercept=0, size=1.5, color="gray") +
      xlim(0, 1) + ylim(0,1) +
      #geom_rug(aes(x=rug, y=0, color=analysis), sides = "b", size=0.5)+
      #geom_point(size=1, alpha=0.5) + 
      #geom_line(size=1.1) + 
      geom_smooth(aes(x=Pred, y=Obs, color=analysis), method="lm", se=FALSE, formula=y~rcs(x, knots)) + 
      geom_density(aes(x=rug, color=analysis, fill=analysis, y=after_stat(scaled/(10))), alpha=0.4) +
      #scale_color_manual(values=npgcols) +
      #annotation_custom(tableGrob(pobjtable, theme=ttheme_minimal(rowhead=list(fg_params = list(col = c(npgcols[3],npgcols[1], npgcols[2])))), rows=analyses), xmin=0.15, xmax=0.2, ymin=0.8, ymax=0.9) + 
      theme_classic() +
      theme(legend.position=legend_on, 
            plot.title = element_text(size=20))
  } else {
    ploplot <- preploplot %>% ggplot(aes(x=Pred, y=Obs, color=analysis)) + 
      geom_abline(slope=1, intercept=0, size=1.5, color="gray") +
      xlim(0, 1) + ylim(0,1) +
      geom_rug(aes(x=rug, y=0, color=analysis), sides = "b")+
      geom_line(size=1.1) + 
      geom_point(aes(x=Pred, y=Obs, color=analysis), fill="white", shape=21, size=2, alpha=0.8) + 
      #geom_smooth(method="lm", se=FALSE, formula=y~poly(x, knots)) + 
      #scale_color_manual(values=npgcols) +
      #annotation_custom(tableGrob(pobjtable, theme=ttheme_minimal(rowhead=list(fg_params = list(col = c(npgcols[3],npgcols[1], npgcols[2])))), rows=analyses), xmin=0.15, xmax=0.2, ymin=0.8, ymax=0.9) + 
      theme_classic() +
      theme(legend.position=legend_on, 
            plot.title = element_text(size=20))
  }
  
  ploplot <- ploplot + xlab("Predicted Event Risk") + ylab("Observed Event Rate")
  
  return(ploplot)
}
model_train <- function(tdf, outcome, outcome_mtype, mvars, engines, penalties, alphas, mnames){
  if(outcome_mtype=="binomial"){
    models=list()
    i=1
    while(length(mvars)-i >= 0){
      if(engines[i]=="glm"){
        models[[i]] <- glm(as.formula(paste0(outcome,"~", paste(unique(unlist(mvars[[i]])), collapse="+"))), data=as.data.frame(tdf), family="binomial", x=T)
      }
      if(engines[i]=="cv.glmnet"){
        models_temp <- cv.glmnet(x=tdf[,mvars[[i]]], y=tdf[,outcome], family="binomial", keep=TRUE, relax=FALSE, penalty.factor=penalties[[i]], alpha=alphas[i])
        models[[i]] <- list()
        models[[i]][[1]] <- models_temp
        models[[i]][[2]] <- GLMnet(formula=as.formula(paste0(outcome,"~", paste(unique(unlist(mvars[[i]])), collapse="+"))), data=as.data.frame(tdf[,c(outcome, mvars[[i]])]), family="binomial", keep=TRUE, relax=FALSE, penalty.factor=c(0, penalties[[i]]), alpha=alphas[i], foldid=models_temp$foldid)
      }
      i=i+1
    }
  }
  if(outcome_mtype=="cox"){
    i=1
    models=list()
    while(length(mvars)-i >= 0){
      if(engines[i]=="coxph"){
        models[[i]] <- coxph(as.formula(paste0(outcome,"~", paste(unique(unlist(mvars[[i]])), collapse="+"))), data=as.data.frame(tdf), x=T)
      }
      if(engines[i]=="cv.glmnet"){
        models_temp <- cv.glmnet(x=tdf[,mvars[[i]]], y=cbind(time=tdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcome), "\\s+"))[1]], status=tdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcome), "\\s+"))[2]]), family="cox", keep=TRUE, relax=FALSE, penalty.factor=penalties[[i]], alpha=alphas[i])
        models[[i]] <- list()
        models[[i]][[1]] <- models_temp
        models[[i]][[2]] <- GLMnet(formula=as.formula(paste0(outcome,"~", paste(unique(unlist(mvars[[i]])), collapse="+"))), data=as.data.frame(tdf[,c(unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcome), "\\s+")), mvars[[i]])]),  family="cox", keep=TRUE, relax=FALSE, penalty.factor=penalties[[i]], alpha=alphas[i], foldid=models_temp$foldid)
      }
      i=i+1
    }
  }
  names(models) <- mnames
  return(models)
}
recal <- function(models, engines, mvars, vdf, outcome, method, split.recal, folds, vfold, seed){
  if(length(unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+")))==1) {
    otype="non-tte"
    time=rep(1, dim(vdf)[1])
    status=vdf[,outcome]
  } else{
    otype="tte"
    time=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[1]]
    status=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]]
  }
  models1 <- lapply(models, function(x) if("list" %in% class(x)){x[[1]]} else{x})
  vdf=as.data.frame(vdf)
  set.seed(seed)
  if(method=="simple" & otype=="non-tte") {
    vdf$split <- rep(2, dim(vdf)[1])
    vdf[caret::createDataPartition(y = factor(vdf[,outcome]), times=1, p = split.recal)[[1]], "split"] <- 1
    vdf_cal <- vdf[vdf$split==1,] ; vdf_val <- vdf[vdf$split==2,]
    recmodels <- list()
    i=1
    while(length(models1)-i>=0){
      if(engines[i]=="glm"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_cal[,mvars[[i]]])  %*% coef(models1[[i]])[-1])
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_val[,mvars[[i]]])  %*% coef(models1[[i]])[-1])
        recmodels[[i]] <- glm(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", outcome))), data=vdf_cal, family="binomial", x=T)
      }
      if(engines[i]=="cv.glmnet"){
        vdf_cal[,paste0("pred", i, "_", names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_cal[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        vdf_val[,paste0("pred", i, "_", names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_val[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        recmodels[[i]] <- glm(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", outcome))), data=vdf_cal, family="binomial", x=T)
      }
      i=i+1
    }
    
  }
  if(method=="simple" & otype=="tte") {
    vdf$split <- rep(2, dim(vdf)[1])
    vdf$split <- coxsplity(y=as.matrix(Surv(time, status)), method = "simple", prop = split.recal)
    vdf_cal <- vdf[vdf$split==1,] ; vdf_val <- vdf[vdf$split==2,]
    recmodels <- list()
    i=1
    while(length(models1)-i>=0){
      if(engines[i]=="coxph"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_cal[,mvars[[i]]])  %*% coef(models1[[i]]))
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_val[,mvars[[i]]])  %*% coef(models1[[i]]))
        recmodels[[i]] <- coxph(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]))), data=vdf_cal, x=T)
      }
      if(engines[i]=="cv.glmnet"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_cal[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_val[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        recmodels[[i]] <- coxph(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]))), data=vdf_cal, x=T)
      }
      i=i+1
    }
  }
  if(method=="crossfold" & otype=="non-tte") {
    vdf$split <- caret::createFolds(factor(vdf$everCR_100),k=folds,list=FALSE)
    vdf_cal <- vdf[vdf$split!=vfold,] ; vdf_val <- vdf[vdf$split==vfold,]
    recmodels <- list()
    i=1
    while(length(models1)-i>=0){
      if(engines[i]=="glm"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_cal[,mvars[[i]]])  %*% coef(models1[[i]])[-1])
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_val[,mvars[[i]]])  %*% coef(models1[[i]])[-1])
        recmodels[[i]] <- glm(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", outcome))), data=vdf_cal, family="binomial", x=T)
      }
      if(engines[i]=="cv.glmnet"){
        vdf_cal[,paste0("pred", i, "_", names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_cal[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        vdf_val[,paste0("pred", i, "_", names(models1)[i], "_", outcome)] <- as.numeric(as.matrix(vdf_val[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        recmodels[[i]] <- glm(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", outcome))), data=vdf_cal, family="binomial", x=T)
      }
      i=i+1
    }
    
  }
  if(method=="crossfold" & otype=="tte") {
    vdf$split <- coxsplity(y=as.matrix(Surv(time, status)), method = "crossfold", nfolds = folds)
    vdf_cal <- vdf[vdf$split!=vfold,] ; vdf_val <- vdf[vdf$split==vfold,]
    recmodels <- list()
    i=1
    while(length(models1)-i>=0){
      if(engines[i]=="coxph"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_cal[,mvars[[i]]])  %*% coef(models1[[i]]))
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_val[,mvars[[i]]])  %*% coef(models1[[i]]))
        recmodels[[i]] <- coxph(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]))), data=vdf_cal, x=T)
      }
      if(engines[i]=="cv.glmnet"){
        vdf_cal[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_cal[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        vdf_val[,paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2])] <- as.numeric(as.matrix(vdf_val[,names(models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])])  %*% models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]][models1[[i]]$glmnet.fit$beta[,models1[[i]]$index[1]]!=0])
        recmodels[[i]] <- coxph(formula=as.formula(paste0(outcome,"~",paste0("pred", i, "_",names(models1)[i], "_", unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]))), data=vdf_cal, x=T)
      }
      i=i+1
    }
  }
  
  names(recmodels) <- names(models1)
  return(list(recalibrated_models=recmodels, vdf_cal=vdf_cal, vdf_val=vdf_val, recal_predvars=colnames(vdf_cal)[grep("pred", colnames(vdf_cal))]))
}
discrim <- function(models, engines,  mvars, vdf, outcome, refmodel, teval){
  if(length(unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+")))==1) {
    otype="non-tte"
    time=rep(1, dim(vdf)[1])
    status=vdf[,outcome]
  } else{
    otype="tte"
    time=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[1]]
    status=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]]
  }
  if(missing(teval) & otype=="tte"){
    predictions <- list()
    i=1
    while(length(models)-i>=0){
      if(engines[i]=="glm"){
        predictions[[i]] <- predict.glm(models[[i]], newdata=as.data.frame(vdf), type="response")
        compareC(timeX=rep(1, dim(vdf)[1]), statusX = vdf[,outcome], scoreY=predict.glm(models[[1]], newdata=as.data.frame(vdf), type="response"), scoreZ=predict.glm(models[[3]], newdata=as.data.frame(vdf), type="response"))
      }
      if(engines[i]=="cv.glmnet"){
        predictions[[i]] <- predict(models[[i]][[1]], newx=vdf[,mvars[[i]]], s="lambda.min")
      }
      if(engines[i]=="coxph"){
        predictions[[i]] <- predict(models[[i]], newdata=as.data.frame(vdf))
      }
      i=i+1
    }
    return(data.frame(
      models=names(models),
      predictors=sapply(mvars, function(x) paste(x, collapse="+")),
      engine=engines,
      cstat=sapply(predictions, function(x) 1-compareC(timeX=time, statusX = status, scoreY=x, scoreZ=predictions[[refmodel]])$est.c["Cxy"]),
      cstat.se=sapply(predictions, function(x) 1.96*sqrt(compareC(timeX=time, statusX = status, scoreY=x, scoreZ=predictions[[refmodel]])$est.varCxy)),
      cstat.pval=sapply(predictions, function(x) compareC(timeX=time, statusX = status, scoreY=x, scoreZ=predictions[[refmodel]])$pval)
    ) %>% mutate(outcome=outcome) %>% reframe(models, outcome, predictors, engine, cstat, cstat.lower=cstat-cstat.se, cstat.upper=cstat+cstat.se, cstat.pval))
    
  }
  
  if(missing(teval) & otype=="non-tte"){
    models1 <- lapply(models, function(x) if("list" %in% class(x)){x[[2]]} else{x})
    names(models1) <- names(models)
    sc <- Score(models1,
                data=vdf,
                formula=as.formula(paste0(outcome,"~1")),
                null.model=FALSE,
                x=TRUE, keep=TRUE, contrasts=lapply(1:length(models), function(x) c(refmodel, x))
    )
    
    return(
      data.frame(
        models=names(models1),
        predictors=sapply(mvars, function(x) paste(x, collapse="+")),
        engine=engines,
        mname=names(models1),
        cstat=sc$AUC$score$AUC,
        cstat.lower=sc$AUC$score$lower,
        cstat.upper=sc$AUC$score$upper,
        cstat.pval=ifelse(length(sc$AUC$score$AUC)==1, NA, sc$AUC$contrasts$p),
        brier=sc$Brier$score$Brier,
        brier.lower=sc$Brier$score$lower,
        brier.upper=sc$Brier$score$upper,
        brier.pval=ifelse(length(sc$AUC$score$AUC)==1, NA, sc$Brier$contrasts$p)
      ) %>% mutate(outcome=outcome) %>% reframe(models, outcome, predictors, engine, cstat, cstat.lower, cstat.upper, cstat.pval, brier, brier.lower, brier.upper, brier.pval)
    )
  }
  
  if(!missing(teval) & otype=="tte"){
    models1 <- lapply(models, function(x) if("list" %in% class(x)){x[[2]]} else{x})
    names(models1) <- names(models)
    sc <- Score(models1,
                data=vdf,
                formula=as.formula(paste0(outcome,"~1")),
                null.model=FALSE,
                x=TRUE, keep=TRUE, contrasts=lapply(1:length(models), function(x) c(refmodel, x)),
                times=teval
    )
    
    return(data.frame(
      models=names(models1),
      predictors=sapply(mvars, function(x) paste(x, collapse="+")),
      engine=engines,
      mname=names(models1),
      cstat=sc$AUC$score$AUC,
      cstat.lower=sc$AUC$score$lower,
      cstat.upper=sc$AUC$score$upper,
      cstat.pval=ifelse(length(sc$AUC$score$AUC)==1, NA, sc$AUC$contrasts$p),
      brier=sc$Brier$score$Brier,
      brier.lower=sc$Brier$score$lower,
      brier.upper=sc$Brier$score$upper,
      brier.pval=ifelse(length(sc$AUC$score$AUC)==1, NA, sc$Brier$contrasts$p)
    ) %>% mutate(outcome=outcome, time=teval) %>% reframe(models, outcome, predictors, engine, time, cstat, cstat.lower, cstat.upper, cstat.pval, brier, brier.lower, brier.upper, brier.pval))
  }
  
}
calib <- function(models, engines,  mvars, vdf, outcome, teval, binmeth, quant, knots, mcolors){
  if(length(unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+")))==1) {
    otype="non-tte"
    time=rep(1, dim(vdf)[1])
    status=vdf[,outcome]
  } else{
    otype="tte"
    time=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[1]]
    status=vdf[,unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+"))[2]]
  }
  
  if(missing(teval) & otype=="non-tte"){
    models1 <- lapply(models, function(x) if("list" %in% class(x)){x[[2]]} else{x})
    names(models1) <- names(models)
    cr_preds <- lapply(models1, function(x) predict.glm(x, newdata=as.data.frame(vdf), type="response"))
    id_link <- bind_rows(map2(cr_preds, names(cr_preds), function(df, name) {
      return(data.frame(id=names(df), risk=df, model=rep(name, length(df))))
    })) %>% rownames_to_column("rm") %>%  left_join(vdf %>% select(!contains("pred")) %>% rownames_to_column("id")) %>%
      rename(record_id=id) %>% select(!c(rm))
    id_link$status=id_link[,outcome]
    
    sc <- Score(models1,
                data=vdf,
                formula=as.formula(paste0(outcome,"~1")),
                null.model=FALSE, plots="cal",
                x=TRUE, keep=TRUE
    )
    p_sc <- plotCalibration(sc, method=binmeth, q=quant, cens.method="local", pseudo=F, rug=T, plot=F)
    caldfsc <- sc$Calibration$plotframe %>% 
      left_join(id_link %>% rename(ReSpOnSe=status), by=c("model", "risk", "ReSpOnSe")) %>% distinct()
    caldf <- collect_riskReg(sc, p_sc, analyses=names(models1)) %>% 
      reframe(model=analysis, binpred=Pred, binobs=Obs, risk=rug, time, status, pseudovalue)
    calconc <- caldf[,c("model", "binpred", "binobs")] %>% group_by(model) %>% ccc(binobs, binpred) %>% reframe(model, ccc=.estimate, outcome=outcome, time=NA)
    caldf <- caldf %>% select(!c(binpred, binobs)) %>% 
      left_join(id_link, by=c("model", "risk", "status")) %>% distinct()
    calplot <- fit_riskReg_spline(score=sc, pobject = p_sc, analyses=names(models1), knots=knots) + scale_color_manual(values=mcolors)
    return(list(caldf=caldf, calplot=calplot, calconc=calconc, caldfsc=caldfsc, sc=sc))
  }
  
  if(!missing(teval) & otype=="tte"){
    models1 <- lapply(models, function(x) if("list" %in% class(x)){x[[2]]} else{x})
    names(models1) <- names(models)
    vdft <- vdf; vdft$time <- time; vdft$status <- status
    id_link <- bind_rows(map2(models1, names(models1), function(mx, mnamex){
      df <- as.data.frame(1-t(summary(survfit(mx, newdata=vdf), times=teval)$surv)) %>% rownames_to_column("record_id") %>% rename(risk=V1) %>% mutate(model=mnamex)
      return(df)
    })) %>% left_join(vdft %>% rownames_to_column("record_id"), by="record_id")
    
    sc <- Score(models1,
                data=vdf,
                formula=as.formula(paste0(outcome,"~1")),
                null.model=FALSE, plots="cal",
                x=TRUE, keep=TRUE,
                times=teval
    )
    p_sc <- plotCalibration(sc, method=binmeth, q=quant, cens.method="local",pseudo = T, rug=T, plot=F)
    caldf <- collect_riskReg(sc, p_sc, analyses=names(models1)) %>% reframe(model=analysis, binpred=Pred, binobs=Obs, risk=rug, risk_link=signif(risk, 7), time, status, pseudovalue)
    calconc <- caldf[,c("model", "binpred", "binobs")] %>% group_by(model) %>% ccc(binobs, binpred) %>% reframe(model, ccc=ifelse(is.na(.estimate), 0, .estimate), outcome=outcome, time=teval)
    caldfsc <- sc$Calibration$plotframe %>% 
      mutate(risk_link=signif(risk, 7)) %>% 
      left_join(id_link %>% mutate(risk_link=signif(risk, 7)), by=c("model", "risk_link", "time", "status")) %>% mutate(risk=coalesce(risk.x, risk.y)) %>% select(!c(risk.x, risk.y, risk_link)) %>% distinct()
    
    caldf <- caldf %>% select(!c(binpred, binobs)) %>% 
      left_join(id_link %>% mutate(risk_link=signif(risk, 7)), by=c("model", "risk_link", "time", "status")) %>% mutate(risk=coalesce(risk.x, risk.y)) %>% select(!c(risk.x, risk.y, risk_link)) %>% distinct()
    
    calplot <- fit_riskReg_spline(score=sc, pobject = p_sc, analyses=names(models1), knots=knots) + scale_color_manual(values=mcolors)
    return(list(caldf=caldf, calplot=calplot, calconc=calconc, caldfsc=caldfsc, sc=sc))
  }
  
}
decurves <- function(caldf, outcome, teval, dca_range){
  if(length(unlist(str_split(gsub("Surv|\\(|\\)|\\,","", outcome), "\\s+")))==1) {
    otype="non-tte"
  } else{
    otype="tte"
  }
  if(missing(teval) & otype=="non-tte"){
    dcap <- dca(as.formula(paste0("status~", paste(unique(caldf$model), collapse="+"))),
                data=caldf %>% select(record_id, time, status, pseudovalue, risk, model) %>% pivot_wider(id_cols=c(record_id, time, status, pseudovalue), names_from = "model", values_from="risk"),
                thresholds=1:100/100) %>% plot(smooth=TRUE) 
    nb <- list()
    dcadf <- dcap$data
    i=1
    while(length(dca_range)-i >=0){
      #nb[[i]] <- dcap$data %>% filter(threshold >= dca_range[[i]][1] & threshold <= dca_range[[i]][2]) %>% group_by(label) %>% dplyr::summarize(menb=mean(net_benefit)) %>% ungroup() %>% mutate(outcome=outcome, time=NA, dcrange=paste(dca_range[[i]], collapse="-"))
      nb[[i]] <- dcap$data %>%  group_by(label) %>% dplyr::summarize(menb=MESS::auc(x=threshold, y=net_benefit, from=dca_range[[i]][1], to=dca_range[[i]][2])) %>% ungroup() %>% mutate(outcome=outcome, time=NA, dcrange=paste(dca_range[[i]], collapse="-"))
      i=i+1
    }
    nb <- bind_rows(nb)
    return(list(dcap=dcap, nb=nb, dcadf=dcadf))
  }
  
  if(!missing(teval) & otype=="tte"){
    dcap <- dca(as.formula(paste0("Surv(time, status)~", paste(unique(caldf$model), collapse="+"))),
                data=caldf %>% select(record_id, time, status, pseudovalue, risk, model) %>% pivot_wider(id_cols=c(record_id, time, status, pseudovalue), names_from = "model", values_from="risk"),
                time=teval,
                thresholds=1:100/100) %>% plot(smooth=TRUE) 
    nb <- list()
    dcadf <- dcap$data
    i=1
    while(length(dca_range)-i >=0){
      #nb[[i]] <- dcap$data %>% filter(threshold >= dca_range[[i]][1] & threshold <= dca_range[[i]][2]) %>% group_by(label) %>% dplyr::summarize(menb=mean(net_benefit)) %>% ungroup() %>% mutate(outcome=outcome, time=teval, dcrange=paste(dca_range[[i]], collapse="-"))
      nb[[i]] <- dcap$data %>%  group_by(label) %>% dplyr::summarize(menb=MESS::auc(x=threshold, y=net_benefit, from=dca_range[[i]][1], to=dca_range[[i]][2])) %>% ungroup() %>% mutate(outcome=outcome, time=NA, dcrange=paste(dca_range[[i]], collapse="-"))
      i=i+1
    }
    nb <- bind_rows(nb)
    return(list(dcap=dcap, nb=nb, dcadf=dcadf))
  }
  
}

# Setup -------------------------------------------------------------------

df_all <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id") %>%
  mutate(bl_crp=noLCS_preld_crp/uln_preld_crp) 

df_train <- df_all %>%
  filter(cohort %in% c("MSK Development"))

df_val <- df_all %>%
  filter(dx_simple=="LBCL") %>%
  filter(!(cohort %in% c("MSK Development"))) 

inflamix_assignments <- reformat_table_tp(df_all, gmm_mu=im_mu, gmm_sig=im_sig, gmm_pro = im_pro, inflammClus = ind_inflclus, least_inflammClus = ind_ninflclus, tp="d0_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
  mutate(tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% 
  select(record_id, tau) %>%
  mutate(logtau=log10(tau))

nocytos_assignments <- reformat_table_tp_nocytos(df=df_all, gmm_mu=dev_mu_nocytos, gmm_sig=dev_sig_nocytos, gmm_pro = dev_pro_nocytos, inflammClus = dev_inflammClus_nocytos, least_inflammClus = dev_LEASTinflammClus_nocytos, tp="d0_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
  mutate(nocytos_tau=ifelse(cluster=="Inflammatory", tau, 1-tau)) %>% 
  select(record_id, nocytos_tau) %>%
  mutate(nocytos_logtau=log10(nocytos_tau))

# Analysis_Vars -----------------------------------------------------------
# Constant Vars
basevars <-  c("age", "primary_ref", "bin_preld_ldh", "costim")
imvars <- c("logtau")
lab14vars <- c("d0_albumin", 
               "d0_alk", 
               "d0_ast", 
               "d0_ferritin", 
               "d0_hb", 
               "d0_ldh", 
               "d0_plt", "d0_tbr",
               "d0_il10", "d0_il6", "d0_tnfa",
               "d0_crp",
               "d0_ddimer", "d0_wbc"
)
lab11vars <- c("d0_albumin", 
               "d0_alk", 
               "d0_ast", 
               "d0_ferritin", 
               "d0_hb", 
               "d0_ldh", 
               "d0_plt", "d0_tbr",
               "d0_crp",
               "d0_ddimer", "d0_wbc"
)
outcomes=c("everCR_100", "Surv(tt_pfs_m, ev_pfs)", "Surv(tt_os_m, ev_os)")
crpvars <- c("bl_crp")
nocytos_vars <- c("nocytos_logtau")

training_pre <- df_train 
validation_pre <- df_val %>% 
  filter( (center %in% c("Memorial Sloan Kettering", "Hackensack Meridian Health", "Sheba Medical Center"))) %>% filter(costim !="Other")

splitmethod="crossfold"
nfolds=2
splitprop=0.5
evaltime=6
knots=3
dcrange= list(c(0.15, 0.30), c(0.30, 0.45))
specvar="response_28_4lvl"
svarfilt=list(c("PR", "SD"), c("CR", "PR"))
iterations=100
legend_on="none"


# Define model comparisons ------------------------------------------------
### For this section of the code - only run the lines for the model comparisons of interest. 

#### InflaMix Alone
modelcols <- c("#E64B35FF")
model.vars=list(
  c(basevars, imvars)
)
model.names=c("InflaMix")
reference.model.index=1
resp_engines = c("glm")
resp_recal_engines = c("glm")
surv_engines = c("coxph")
surv_recal_engines = c("coxph")
modelreg.penalties=list(NA)
modelreg.alpha=c(NA)
trainregseed = 1

#### CRP Benchmark
modelcols <- c("#4DBBD5FF", "#00A087FF", "#E64B35FF")
model.vars=list(
  c(basevars),
  c(basevars, crpvars),
  c(basevars, imvars)
)
model.names=c("Base", "CRP", "InflaMix")
reference.model.index=3
resp_engines = c("glm", "glm", "glm")
resp_recal_engines = c("glm", "glm", "glm")
surv_engines = c("coxph", "coxph", "coxph")
surv_recal_engines = c("coxph", "coxph", "coxph")
modelreg.penalties=list(NA, NA,  NA)
modelreg.alpha=c(NA, NA, NA)
trainregseed = 1

#### CRP + NoCytos Benchmark
modelcols <- c("#4DBBD5FF", "#00A087FF",  "#F39B7FFF", "#E64B35FF")
model.vars=list(
  c(basevars),
  c(basevars, crpvars),
  c(basevars, nocytos_vars),
  c(basevars, imvars)
)
model.names=c("Base", "CRP", "NoCytosMM", "InflaMix")
reference.model.index=4
resp_engines = c("glm", "glm", "glm", "glm")
resp_recal_engines = c("glm", "glm", "glm", "glm")
surv_engines = c("coxph", "coxph", "coxph", "coxph")
surv_recal_engines = c("coxph", "coxph", "coxph", "coxph")
modelreg.penalties=list(NA, NA,  NA, NA)
modelreg.alpha=c(NA, NA, NA, NA)
trainregseed = 1

#### Lab 14Reg Only Benchmark
modelcols <- c("#8491B4FF", "#E64B35FF")
model.vars=list(
  c(basevars, lab14vars),
  c(basevars, imvars)
)
model.names=c("Lab14Reg", "InflaMix")
reference.model.index=2
resp_engines = c("cv.glmnet", "glm")
resp_recal_engines = c("glm", "glm")
surv_engines = c("cv.glmnet", "coxph")
surv_recal_engines = c("coxph", "coxph")
modelreg.penalties=list(c(rep(0, 4), rep(1, 14)),  NA)
modelreg.alpha=c(1, NA)
trainregseed = 1

#### Lab 11Reg Benchmark
modelcols <- c("#8491B4FF", "#F39B7FFF" , "#E64B35FF")
model.vars=list(
  c(basevars, lab11vars),
  c(basevars, nocytos_vars),
  c(basevars, imvars)
)
model.names=c("Lab11Reg", "NoCytosMM", "InflaMix")
reference.model.index=3
resp_engines = c("cv.glmnet", "glm", "glm")
resp_recal_engines = c("glm", "glm", "glm")
surv_engines = c("cv.glmnet", "coxph", "coxph")
surv_recal_engines = c("coxph","coxph", "coxph")
modelreg.penalties=list(c(rep(0, 4), rep(1, 11)), NA,  NA)
modelreg.alpha=c(1, NA, NA)
trainregseed = 1

# Common Code to Run: Loop ------------------------------------------------------
totalcohort <- bind_rows(training_pre, validation_pre) %>% left_join(inflamix_assignments) %>% left_join(nocytos_assignments) %>% 
  mutate(svar=!!sym(specvar)) %>% 
  select(record_id, everCR_100, tt_pfs_m, ev_pfs, tt_os_m, ev_os, starts_with(unique(unlist(model.vars))), svar) %>%
  drop_na(!svar) %>%
  mutate(svar1=ifelse(is.na(svar), NA, ifelse(svar %in% svarfilt[[1]], 1, 0))) %>%
  mutate(svar2=ifelse(is.na(svar), NA, ifelse(svar %in% svarfilt[[2]], 1, 0))) %>%
  select(!svar)

training <- totalcohort %>% filter(record_id %in% training_pre$record_id) 
rownames(training) <- training$record_id
validation <- totalcohort %>% filter(record_id %in% validation_pre$record_id) 
rownames(validation) <- validation$record_id

tmat <- model.matrix(as.formula(paste0(
  "~",
  paste(unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcomes), "\\s+")), collapse = "+"),
  "+",
  paste(unique(unlist(model.vars)), collapse="+")
)),
training
) [,-1]

colnames(tmat) <- c(unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcomes), "\\s+")), unique(unlist(model.vars)))
tmat <- as.matrix(cbind(tmat,
                        (totalcohort %>% filter(record_id %in% rownames(tmat)) %>% column_to_rownames("record_id") %>% select(contains("svar")))
))

vmat <- model.matrix(as.formula(paste0(
  "~",
  paste(unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcomes), "\\s+")), collapse = "+"),
  "+",
  paste(unique(unlist(model.vars)), collapse="+")
)),
validation
) [,-1]
colnames(vmat) <- c(unlist(str_split(gsub("Surv|\\(|\\)|\\,", "", outcomes), "\\s+")), unique(unlist(model.vars)))
vmat <- as.matrix(cbind(vmat,
                        (totalcohort %>% filter(record_id %in% rownames(vmat)) %>% column_to_rownames("record_id") %>% select(contains("svar")))
))

# Train Models
set.seed(trainregseed)
cr_models <- model_train(
  tdf=tmat,
  outcome="everCR_100",
  outcome_mtype="binomial",
  mvars=model.vars,
  engines=resp_engines,
  penalties=modelreg.penalties,
  alpha=modelreg.alpha,
  mnames=model.names
)

net_cr_discrim <- discrim(models=cr_models, engines = resp_engines, mvars = model.vars, vdf = vmat, outcome = "everCR_100", refmodel = reference.model.index)


pfs_models <- model_train(
  tdf=tmat,
  outcome="Surv(tt_pfs_m, ev_pfs)",
  outcome_mtype="cox",
  mvars=model.vars,
  engines=surv_engines,
  penalties=modelreg.penalties,
  alpha=modelreg.alpha,
  mnames=model.names
)

net_pfs_discrim_cindex <- discrim(models=pfs_models, engines = surv_engines, mvars = model.vars, vdf = vmat, outcome = "Surv(tt_pfs_m, ev_pfs)", refmodel = reference.model.index)
net_pfs_discrim <- discrim(models=pfs_models, engines = surv_engines, mvars = model.vars, vdf = vmat, outcome = "Surv(tt_pfs_m, ev_pfs)", refmodel = reference.model.index, teval=evaltime)

os_models <- model_train(
  tdf=tmat,
  outcome="Surv(tt_os_m, ev_os)",
  outcome_mtype="cox",
  mvars=model.vars,
  engines=surv_engines,
  penalties=modelreg.penalties,
  alpha=modelreg.alpha,
  mnames=model.names
)

net_os_discrim_cindex <- discrim(models=os_models, engines = surv_engines, mvars = model.vars, vdf = vmat, outcome = "Surv(tt_os_m, ev_os)", refmodel = reference.model.index)
net_os_discrim <- discrim(models=os_models, engines = surv_engines, mvars = model.vars, vdf = vmat, outcome = "Surv(tt_os_m, ev_os)", refmodel = reference.model.index, teval=evaltime*2)

net_discrim <- bind_rows(
  net_cr_discrim,
  net_pfs_discrim_cindex,
  net_pfs_discrim,
  net_os_discrim_cindex,
  net_os_discrim
) %>% select(!c(predictors, engine))

# Recalibration 
n_performance<- list()
n_calibration<- list()
s=1
while(iterations-sum(unlist(lapply(n_performance,function(x) length(x)==nfolds))) > 0){
  try({
    n_performance[[s]] <- list()
    n_calibration[[s]] <- list()
    for(i in 1:nfolds){
      recal_crmodels <- recal(
        models=cr_models,
        engines=resp_engines,
        mvars=model.vars,
        vdf=vmat,
        outcome="everCR_100",
        method=splitmethod,
        folds=nfolds,
        vfold=i,
        seed=s
      )
      
      recal_pfsmodels <- recal(
        models=pfs_models,
        engines=surv_engines,
        mvars=model.vars,
        vdf=vmat,
        outcome="Surv(tt_pfs_m, ev_pfs)",
        method=splitmethod,
        folds=nfolds,
        vfold=i,
        seed=s
      )
      
      recal_osmodels <- recal(
        models=os_models,
        engines=surv_engines,
        mvars=model.vars,
        vdf=vmat,
        outcome="Surv(tt_os_m, ev_os)",
        method=splitmethod,
        folds=nfolds,
        vfold=i,
        seed=s
      )
      
      pfscalibs <- calib(models=recal_pfsmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_pfsmodels$recal_predvars, vdf = recal_pfsmodels$vdf_val, outcome = "Surv(tt_pfs_m, ev_pfs)", binmeth = "nne", quant = 3, knots = knots, teval=evaltime, mcolors = modelcols)
      pfsdca <- decurves(caldf = pfscalibs$caldf, outcome = "Surv(tt_pfs_m, ev_pfs)", teval=evaltime, dca_range = dcrange)
      pfsdca_spec1 <- decurves(caldf = pfscalibs$caldf %>% filter(svar1==1), outcome = "Surv(tt_pfs_m, ev_pfs)", teval=evaltime, dca_range = dcrange)
      pfsdca_spec2 <- decurves(caldf = pfscalibs$caldf %>% filter(svar2==1), outcome = "Surv(tt_pfs_m, ev_pfs)", teval=evaltime, dca_range = dcrange)
      
      
      crcalibs <- calib(models=recal_crmodels$recalibrated_models, engines = resp_recal_engines, mvars = recal_crmodels$recal_predvars, vdf = recal_crmodels$vdf_val, outcome = "everCR_100", binmeth = "nne", quant = 3, knots = knots, mcolors = modelcols)
      crdca <- decurves(caldf = crcalibs$caldf, outcome = "everCR_100", dca_range = dcrange)
      crdca_spec1 <- decurves(caldf = crcalibs$caldf %>% filter(svar1==1), outcome = "everCR_100", dca_range = dcrange)
      crdca_spec2 <- decurves(caldf = crcalibs$caldf %>% filter(svar2==1), outcome = "everCR_100", dca_range = dcrange)
      
      
      oscalibs <- calib(models=recal_osmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_osmodels$recal_predvars, vdf = recal_osmodels$vdf_val, outcome = "Surv(tt_os_m, ev_os)", binmeth = "nne", quant = 3, knots = knots, teval=evaltime*2, mcolors = modelcols)
      #osdca <- decurves(caldf = oscalibs$caldf, outcome = "Surv(tt_os_m, ev_os)", teval=evaltime*2, dca_range = dcrange)
      #osdca_spec1 <- decurves(caldf = oscalibs$caldf %>% filter(svar1==1), outcome = "Surv(tt_os_m, ev_os)", teval=evaltime*2, dca_range = dcrange)
      #osdca_spec2 <- decurves(caldf = oscalibs$caldf %>% filter(svar2==1), outcome = "Surv(tt_os_m, ev_os)", teval=evaltime*2, dca_range = dcrange)
      
      n_calibration[[s]][[i]] <- bind_rows(
        crcalibs$caldfsc %>% mutate(outcome="cr", fold=i, iteration=s),
        pfscalibs$caldfsc %>% mutate(outcome="pfs", fold=i, iteration=s),
        oscalibs$caldfsc %>% mutate(outcome="os", fold=i, iteration=s)
      )
      
      n_performance[[s]][[i]] <- bind_rows(
        discrim(models=recal_crmodels$recalibrated_models, engines = resp_recal_engines, mvars = recal_crmodels$recal_predvars, vdf = recal_crmodels$vdf_val, outcome = "everCR_100", refmodel = reference.model.index),
        discrim(models=recal_pfsmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_pfsmodels$recal_predvars, vdf = recal_pfsmodels$vdf_val, outcome = "Surv(tt_pfs_m, ev_pfs)", refmodel = reference.model.index),
        discrim(models=recal_pfsmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_pfsmodels$recal_predvars, vdf = recal_pfsmodels$vdf_val, outcome = "Surv(tt_pfs_m, ev_pfs)", refmodel = reference.model.index, teval=evaltime),
        discrim(models=recal_osmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_osmodels$recal_predvars, vdf = recal_osmodels$vdf_val, outcome = "Surv(tt_os_m, ev_os)", refmodel = reference.model.index),
        discrim(models=recal_osmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_osmodels$recal_predvars, vdf = recal_osmodels$vdf_val, outcome = "Surv(tt_os_m, ev_os)", refmodel = reference.model.index, teval=evaltime*2),
      ) %>% reframe(model=models, outcome, time, cstat, brier, fold=i, iteration=s) %>%
        left_join(
          bind_rows(
            crcalibs$calconc,
            pfscalibs$calconc,
            oscalibs$calconc
          ) 
        ) %>% full_join(
          bind_rows(
            crdca$nb,
            pfsdca$nb
          )  %>% rename(model=label) %>% mutate(fold=i, iteration=s)
          
        ) %>% left_join(
          bind_rows(
            crdca_spec1$nb %>% rename(menbf1=menb),
            pfsdca_spec1$nb %>% rename(menbf1=menb)
          )  %>% rename(model=label) %>% mutate(fold=i, iteration=s)
        ) %>% left_join(
          bind_rows(
            crdca_spec2$nb %>% rename(menbf2=menb),
            pfsdca_spec2$nb %>% rename(menbf2=menb)
          )  %>% rename(model=label) %>% mutate(fold=i, iteration=s)
        )
    }
  })
  s=s+1
  print(s)
}


# Downstream analysis -----------------------------------------------------

legend_on <- "none"

master_metrics <- bind_rows(lapply(n_performance[c(1:101)], function(x) {bind_rows(x)}))
master_calibration <- bind_rows(lapply(n_calibration[c(1:101)], function(x) {bind_rows(x)}))

master <- master_metrics %>%
  mutate(outcome=case_when(
    outcome=="everCR_100" ~ "No CR by Day 100",
    outcome=="Surv(tt_pfs_m, ev_pfs)" ~ "PFS",
    outcome=="Surv(tt_os_m, ev_os)" ~ "OS",
  )) %>%
  mutate(outcome=paste0(outcome, ifelse(is.na(time), "", paste0(" at ", time, " months")))) %>% select(!time) %>%
  pivot_longer(!c(model, outcome, iteration, dcrange, fold), names_to = "metric", values_to = "value") %>%
  mutate(metric=case_when(
    metric=="menb" ~ "Decision at Pre-Infusion",
    metric=="menbf1" ~ "Decision if PR at Day +30",
    metric=="menbf2" ~ "Decision if CR or PR at Day +30",
    .default = metric
  )) %>% filter(!grepl("OS", outcome)) %>% filter(outcome!="PFS") %>%
  mutate(model=factor(model, levels=model.names))

crsc <- calib(models=recal_crmodels$recalibrated_models, engines = resp_recal_engines, mvars = recal_crmodels$recal_predvars, vdf = recal_crmodels$vdf_val, outcome = "everCR_100", binmeth = "nne", quant = 3, knots = knots, mcolors = modelcols)$sc
pfssc <- calib(models=recal_pfsmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_pfsmodels$recal_predvars, vdf = recal_pfsmodels$vdf_val, outcome = "Surv(tt_pfs_m, ev_pfs)", binmeth = "nne", quant = 3, knots = knots, teval=evaltime, mcolors = modelcols)$sc
ossc <- calib(models=recal_osmodels$recalibrated_models, engines = surv_recal_engines, mvars = recal_osmodels$recal_predvars, vdf = recal_osmodels$vdf_val, outcome = "Surv(tt_os_m, ev_os)", binmeth = "nne", quant = 3, knots = knots, teval=evaltime*2, mcolors = modelcols)$sc

crsc$Calibration$plotframe <- as.data.table(master_calibration %>% filter(outcome=="cr") %>% 
                                              #mutate(ID=record_id) %>% 
                                              #group_by(ID, model) %>% dplyr::summarize(ID, ReSpOnSe=ReSpOnSe, model, risk=median(risk)) %>% ungroup() %>% distinct()%>%
                                              reframe(ID, ReSpOnSe, model, risk))
pcrsc <- plotCalibration(crsc, method="quantile", q=20, cens.method="local", pseudo = T, rug=F, plot=T)
crcal <- fit_riskReg_spline(score = crsc, pobject=pcrsc, analyses=model.names, knots=3, lintype="spline") + scale_color_manual(values = modelcols) + scale_fill_manual(values = modelcols)

bind_rows(map2(pcrsc$plotFrames, names(pcrsc$plotFrames), function(x, y) {x %>% mutate(model=y)})) %>% 
  group_by(model) %>% ccc(truth=Obs, estimate=Pred) %>% ungroup()

pfssc$Calibration$plotframe <- as.data.table(master_calibration %>% filter(outcome=="pfs") %>% 
                                               reframe(ID, times, pseudovalue, time, status, WTi, model, risk, Wt))
ppfssc <- plotCalibration(pfssc, method="quantile", q=20, cens.method="local", pseudo = T, rug=F, plot=T)

bind_rows(map2(ppfssc$plotFrames, names(ppfssc$plotFrames), function(x, y) {x %>% mutate(model=y)})) %>% 
  group_by(model) %>% ccc(truth=Obs, estimate=Pred) %>% ungroup()

legend_on <- "none"
pfscal <- fit_riskReg_spline(score = pfssc, pobject=ppfssc, analyses=model.names, knots=3, lintype="spline") + scale_color_manual(values = modelcols) + scale_fill_manual(values = modelcols)

dcaplot <- dca(as.formula(paste0("Surv(time, status)~", paste(unique((master_calibration %>% filter(outcome=="pfs"))$model), collapse="+"))),
               data= master_calibration %>% filter(outcome=="pfs") %>%
                 filter(svar1==1) %>%
                 select(record_id, time, status, pseudovalue, risk, model, fold, iteration) %>% pivot_wider(id_cols=c(record_id, time, status, pseudovalue, fold, iteration), names_from = "model", values_from="risk"),
               time=6,
               thresholds=1:100/100) %>% plot(smooth=TRUE) +
  scale_color_manual(values=c("violetred4","black", (modelcols))) +
  scale_fill_manual(values=c("violetred4","black", (modelcols))) +
  geom_vline(xintercept = 0.3, size=1.5, color="gray") +
  coord_cartesian(xlim=c(0.10, 0.50), ylim=c(0, 0.35), clip="on") +
  ggtitle("") +  
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        plot.title=element_text(size=14),
        plot.subtitle = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)
        
  )

cstat_median <- master %>% 
  filter(metric=="cstat") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange)) %>% distinct() %>% 
  group_by(model, outcome) %>% dplyr::summarize(med.cstat=median(value)) %>% ungroup() 

cstat_rank <- master %>% 
  filter(metric=="cstat") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange))  %>% distinct() %>% 
  group_by(outcome, iteration, fold) %>% mutate(rank=rank(-as.numeric(value))) %>% ungroup() %>%
  group_by(outcome, model) %>% dplyr::summarize(avg.rank=mean(rank)) %>% ungroup()

plot1 <- 
  plot_grid(
    master %>% 
      filter(metric=="cstat") %>% 
      filter(!is.na(value)) %>% select(!c(dcrange, metric)) %>%
      ggplot(aes(x=model, y=value,color=model, fill=model)) + 
      geom_boxplot(alpha=0.5, outliers = FALSE) + #stat_compare_means(paired=TRUE) +  
      geom_label(data=cstat_median, aes(x=model, y=med.cstat, label=round(med.cstat, 2)), position=position_dodge(width=0.75), vjust=-0.3, fill="white")+
      #geom_text(data=cstat_rank %>% mutate(y=0.88), aes(x=model, y=y, label=round(avg.rank, 2)), position=position_dodge(width=0.75), size=4.5)+
      #annotate(geom="text", x=3, y=0.47, label="Average Rank", size=4.5) + 
      facet_grid(.~outcome) +
      scale_color_manual(values=c(modelcols)) + 
      scale_fill_manual(values=c(modelcols)) + 
      theme_classic() + xlab("") + ylab("C-Statistic") + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=0.01),
            axis.title.y = element_text(size=14),
            axis.text = element_text(size=12)
      ),
    ggplot() + geom_label(data=cstat_rank %>% mutate(y="Mean\nRank"), aes(x=model, y=y, label=round(avg.rank, 2), color=model), position=position_dodge(width=0.75), size=4.5) +
      scale_color_manual(values=modelcols) + 
      facet_grid(.~outcome) +
      scale_color_manual(values=c(modelcols)) + 
      scale_fill_manual(values=c(modelcols)) + 
      theme_classic() + xlab("") + ylab("") + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=0.01),
            axis.title.y = element_text(size=1),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=14),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text  = element_blank()
      ),
    ncol = 1,align = "v", rel_heights = c(5, 1)
  )

brier_median <- master %>% 
  filter(metric=="brier") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange)) %>% distinct() %>% 
  group_by(model, outcome) %>% dplyr::summarize(med.brier=median(value)) %>% ungroup() 

brier_rank <- master %>% 
  filter(metric=="brier") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange))  %>% distinct() %>% 
  group_by(outcome, iteration, fold) %>% mutate(rank=rank(as.numeric(value))) %>% ungroup() %>%
  group_by(outcome, model) %>% dplyr::summarize(avg.rank=mean(rank)) %>% ungroup()

plot2 <- 
  plot_grid(master %>% 
              filter(metric=="brier") %>% 
              filter(!is.na(value)) %>% select(!c(dcrange, metric)) %>% distinct() %>%
              ggplot(aes(x=model, y=value, color=model, fill=model)) + 
              geom_boxplot(alpha=0.5, outliers = FALSE) + #stat_compare_means(paired=TRUE, label.y=0.25) + 
              geom_label(data=brier_median, aes(x=model, y=med.brier, label=round(med.brier, 2)), position=position_dodge(width=0.75), vjust=-0.3, fill="white")+
              #geom_text(data=cstat_rank %>% mutate(y=0.88), aes(x=model, y=y, label=round(avg.rank, 2)), position=position_dodge(width=0.75), size=4.5)+
              #annotate(geom="text", x=3, y=0.47, label="Average Rank", size=4.5) + 
              facet_grid(.~outcome) +
              scale_color_manual(values=c(modelcols)) + 
              scale_fill_manual(values=c(modelcols)) + 
              theme_classic() + xlab("") + ylab("Brier Score") + 
              theme(legend.position = "none",
                    legend.title = element_blank(),
                    legend.text = element_text(size=12),
                    axis.title.x = element_text(size=0.01),
                    axis.title.y = element_text(size=14),
                    axis.text = element_text(size=12)
              ),
            ggplot() + geom_label(data=brier_rank %>% mutate(y="Mean\nRank"), aes(x=model, y=y, label=round(avg.rank, 2), color=model), position=position_dodge(width=0.75), size=4.5) +
              scale_color_manual(values=modelcols) + 
              facet_grid(.~outcome) +
              scale_color_manual(values=c(modelcols)) + 
              scale_fill_manual(values=c(modelcols)) + 
              theme_classic() + xlab("") + ylab("") + 
              theme(legend.position = "none",
                    legend.title = element_blank(),
                    legend.text = element_text(size=12),
                    axis.title.x = element_text(size=0.01),
                    axis.title.y = element_text(size=1),
                    axis.text.x = element_blank(),
                    axis.text.y = element_text(size=14),
                    axis.line.x = element_blank(),
                    axis.line.y = element_blank(),
                    axis.ticks = element_blank(),
                    strip.background = element_blank(),
                    strip.text  = element_blank()
              ),
            ncol = 1, align = "v", rel_heights = c(5, 1)
  )

ccc_median <- master %>% 
  filter(metric=="ccc") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange)) %>% distinct() %>% 
  group_by(model, outcome) %>% dplyr::summarize(med.ccc=median(value)) %>% ungroup() 

master %>% 
  filter(metric=="cstat") %>% 
  filter(outcome=="PFS at 6 months") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange)) %>% distinct() %>% 
  group_by(model, outcome) %>%
  dplyr::summarize(
    lqr=quantile(value, probs=c(0.25, 0.5, 0.75))[1],
    med=quantile(value, probs=c(0.25, 0.5, 0.75))[2],
    uqr=quantile(value, probs=c(0.25, 0.5, 0.75))[3]
  ) %>%
  ungroup() 

ccc_rank <- master %>% 
  filter(metric=="ccc") %>% 
  filter(!is.na(value)) %>% select(!c(metric, dcrange))  %>% distinct() %>% 
  group_by(outcome, iteration, fold) %>% mutate(rank=rank(-as.numeric(value))) %>% ungroup() %>%
  group_by(outcome, model) %>% dplyr::summarize(avg.rank=mean(rank)) %>% ungroup()

plot3 <- 
  plot_grid(
    master %>% 
      filter(metric=="ccc") %>% 
      filter(!is.na(value)) %>% select(!c(dcrange, metric)) %>% distinct() %>%
      ggplot(aes(x=model, y=value, color=model, fill=model)) + 
      geom_boxplot(alpha=0.5, outliers = FALSE) + #stat_compare_means(paired=TRUE, label.y=1.1) + 
      geom_label(data=ccc_median, aes(x=model, y=med.ccc, label=round(med.ccc, 2)), position=position_dodge(width=0.75), vjust=-0.3, fill="white")+
      #geom_text(data=cstat_rank %>% mutate(y=0.88), aes(x=model, y=y, label=round(avg.rank, 2)), position=position_dodge(width=0.75), size=4.5)+
      #annotate(geom="text", x=3, y=0.47, label="Average Rank", size=4.5) + 
      facet_grid(.~outcome) +
      scale_color_manual(values=c(modelcols)) + 
      scale_fill_manual(values=c(modelcols)) + 
      theme_classic() + xlab("") + ylab("C-Statistic") + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=0.01),
            axis.title.y = element_text(size=14),
            axis.text = element_text(size=12)
      ),
    ggplot() + geom_label(data=ccc_rank %>% mutate(y="Mean\nRank"), aes(x=model, y=y, label=round(avg.rank, 2), color=model), position=position_dodge(width=0.75), size=4.5) +
      scale_color_manual(values=modelcols) + 
      facet_grid(.~outcome) +
      scale_color_manual(values=c(modelcols)) + 
      scale_fill_manual(values=c(modelcols)) + 
      theme_classic() + xlab("") + ylab("") + 
      theme(legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=0.01),
            axis.title.y = element_text(size=1),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=14),
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text  = element_blank()
      ),
    ncol = 1,align = "v", rel_heights = c(5, 1)
  )

models=pfs_models
engines=surv_engines
mvars = model.vars
vdf = as.data.frame(vmat)
outcome = "Surv(tt_pfs_m, ev_pfs)"
teval = 6
binmeth = "nne"
quant=3
knots = 3
mcolors=modelcols


#pfs_cal <- calib(models=pfs_models, engines=surv_engines, mvars = model.vars, vdf = as.data.frame(vmat), outcome = "Surv(tt_pfs_m, ev_pfs)", teval = 6, binmeth = "nne", quant = 3, knots = 3, mcolors = modelcols)

set.seed(1)
vdf <- as.data.frame(vmat)
vdf$tt_pfs_m <-vdf$tt_pfs_m+runif(n=length(vdf$tt_pfs_m), min=-1e-6, max=1e-6)
models1 <- lapply(pfs_models, function(x) if("list" %in% class(x)){x[[2]]} else{x})
names(models1) <- names(models)
vdft <- vdf; vdft$time <- vdf$tt_pfs_m; vdft$status <- vdf$ev_pfs

sc <- Score(models1,
            data=vdf,
            formula=as.formula(paste0(outcome,"~1")),
            null.model=FALSE, plots="cal",
            x=TRUE, keep=TRUE,
            times=teval
)
p_sc <- plotCalibration(sc, method=binmeth, q=quant, cens.method="local",pseudo = T, rug=T, plot=T)
caldfsc <- sc$Calibration$plotframe %>% 
  left_join(vdft, by=c("time", "status")) %>% distinct()

dcaplot_nocal <- dca(as.formula(paste0("Surv(time, status)~", paste(unique(caldfsc$model), collapse="+"))),
                     data= caldfsc %>% 
                       filter(svar1==1) %>%
                       select(ID, time, status, pseudovalue, risk, model) %>% pivot_wider(id_cols=c(ID, time, status, pseudovalue), names_from = "model", values_from="risk"),
                     time=6,
                     thresholds=1:100/100) %>% plot(smooth=TRUE) +
  scale_color_manual(values=c("violetred4","black", (modelcols))) +
  geom_vline(xintercept = 0.3, size=1.5, color="gray") +
  coord_cartesian(xlim=c(0.10, 0.50), ylim=c(0, 0.35), clip="on") +
  ggtitle("") +  
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        plot.title=element_text(size=14),
        plot.subtitle = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12)
        
  )

crcal
pfscal
dcaplot
plot1
plot2
plot3
dcaplot_nocal
net_discrim


# End Code ----------------------------------------------------------------





