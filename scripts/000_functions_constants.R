library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)

# Obtain color scales
colorclusters <- function(nclus){
  fullspectrum=c("#9467BDFF", "#2CA02CFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF")
  colorscale=c("#1F77B4FF",  "#FF7F0EFF")
  if(nclus>2){
    colorscale=c("#1F77B4FF", fullspectrum[1:(nclus-2)], "#FF7F0EFF")
  }
  return(colorscale)
}

# Function to assign a cluster from partial labs
assign_clust <- function(x, mu, sig, pro){
  misfeats <- colnames(x)[!c(colnames(x) %in% rownames(as.data.frame(t(x)) %>% filter(!is.na(.))))]
  x1 <- x %>% select(!contains(misfeats))
  mu1 <- mu[!(rownames(mu) %in% misfeats),]
  sig1 <- sig[!(rownames(sig) %in% misfeats),!(rownames(sig) %in% misfeats),]
  cnum <- length(pro)
  densfun <- 1:cnum
  if(sum(dim(x1))==2){
    mu1 <- as.matrix(data.frame(t(mu1)))
    colnames(mu1) <- NULL
    rownames(mu1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    sig1 <- array(sig1, dim=c(1,1,2))
    colnames(sig1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    rownames(sig1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    return(NaN)
  }
  for(i in 1:cnum){
    densfun[i] <- dmvnorm(x1, mean = mu1[,i], sigma = sig1[,,i], log = FALSE)
  }
  return( (pro * densfun) / sum(pro * densfun) )
}

# Function to reformat data table to include cluster assignments with tau probabilities for any timepoint
reformat_table_tp <- function(df, gmm_mu, gmm_sig, gmm_pro, inflammClus, least_inflammClus, tp){
  df_2clust_d0labs <- df %>% select(record_id, starts_with(tp)) %>%
    select(record_id, contains(cluster_labs)) %>% column_to_rownames("record_id") %>%
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

# Plot barplots for binary outcomes
clust_bar_plot <- function(df, tp_pre, list_res, metric, xmi, ymi, sz) {
  colorscale <- colorclusters(length(unique(df$Cluster)))
  df_restox <- bind_rows(
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, everCR_100) %>% reframe(n=n()) %>%
      filter(!is.na(everCR_100)) %>% group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(everCR_100=="CR") %>%
      rename(feat=everCR_100) %>% mutate(feat_name="CR"),
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, crs24) %>% reframe(n=n()) %>%
      filter(!is.na(crs24)) %>%group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(crs24=="CRS >1") %>%
      rename(feat=crs24) %>% mutate(feat_name="CRS"),
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, icans24) %>% reframe(n=n()) %>%
      filter(!is.na(icans24)) %>%group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(icans24=="ICANS >1") %>%
      rename(feat=icans24) %>% mutate(feat_name="ICANS")
  ) %>% filter(feat_name %in% list_res)

  trp <- df_restox %>%
    ggplot(aes(x=feat_name, y=prop, fill=Cluster)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge())+
    geom_label(aes(y=prop/2, group=Cluster, label=paste0(n,"\n",tot)), color="black", fill="white",
               fontface="bold", position = position_dodge(width = 0.9), size=sz*1.5) +
    guides(fill = FALSE) +
    theme_minimal()+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, size=15),
      axis.title.x = element_text(angle = 0, size=15, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
    ) +
    labs(x = "Non-Infl.       Infl.",
         y = paste0("Proportion ", list_res),
         title="",
    )+
    scale_fill_manual(values=colorscale) +
    annotate(geom = "text", x=xmi, y=ymi, label = paste0(
      "Adj. OR (95% CI)\n", metric$expEstimate, " (", metric$low_ci, " - ", metric$high_ci, ")", "\np ",  cpval(metric$pvalue)
    ), size=sz)
  return(trp)
}

# Plot KM curves for survival outcomes
clus_ggsurvplot <- function(df, event, metric, wght, sz, xmi, ymi, labadj, labl){

  if(wght==FALSE){
    df <- df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(tau=1)
  }

  if(event=="pfs"){
    event <- "PFS"
    sp <- survfit(Surv(tt_pfs_m, ev_pfs) ~ Cluster, data = df, weights=tau)
  } else {
    if(event=="os") {
      event <- "OS"
      sp <- survfit(Surv(tt_os_m, ev_os) ~ Cluster, data = df, weights=tau)
    } else {event <- ""}
  }

  timedat <- sp$time
  strat1 <- sp$strata[1]
  strat2 <- sp$strata[2]
  timind1 <- which.min(abs(timedat[1:strat1]-20))
  timind2 <- which.min(abs(timedat[-(1:strat1)]-20))+strat1
  survind1 <- sp$surv[timind1]+0.1
  survind2 <- sp$surv[timind2]+0.1

  survpl <- sp %>% ggsurvfit(linewidth=0.75) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent,
      expand = c(0.01, 0)
    ) +
    coord_cartesian(xlim=c(0, 24), clip="on")+
    scale_x_continuous(
      #limits = c(0, 25),
      breaks = c(0,6,12,18,24),
      labels = c(0, 6,12,18,24)
    ) +
    scale_color_manual(values = colorclusters(2)) +
    theme_minimal() +
    theme(legend.position = FALSE) +
    add_risktable(risktable_height=0.2, size=5,
                  risktable_stats = "n.risk",
                  theme = list( theme_risktable_default(axis.text.y.size = 13, plot.title.size = 13))
    )+
    add_risktable_strata_symbol(symbol = "\U25CF", size = 25)+
    add_censor_mark() +
    labs(
      title = "",
      y = paste0(labl, " Probability"),
      x="Months from CAR-T Infusion"
    ) +
    guides(color=FALSE, linetype=FALSE)+
    theme(
      axis.text.x = element_text(angle = 0, size=15),
      axis.text.y = element_text(angle = 0, size=15),
      axis.title.x = element_text(angle = 0, size=18),
      axis.title.y = element_text(size=18),
      legend.text = element_text(size=16),
      legend.position="bottom"
    ) +
    annotate("text", x=20, y=survind1+labadj, label="Non-Inflammatory", size=sz/3, color=colorclusters(2)[1])+
    annotate("text", x=20, y=survind2+labadj, label="Inflammatory", size=(sz/3)+0.5, color=colorclusters(2)[2])+
    labs(
      title ="",
      y = paste0(event, " Probability"),
      x="Months from Infusion"
    ) +
    annotate(geom = "text", x=xmi, y=ymi, label = paste0(
      "Adj. HR (95% CI)\n", metric$expEstimate, " (", metric$low_ci, " - ", metric$high_ci, ")", "\np ",  cpval(metric$pvalue)
    ), size=6)

  return(survpl)
}

# Calculate Inferences for Validation
clust_metrics <- function(df, dxres, covar, res_covar, aname, wght){

  if(wght==FALSE){
      df <- df %>% group_by(record_id) %>%
        dplyr::slice(which.max(tau)) %>%ungroup() %>% mutate(tau=1)
    df <- df %>% mutate(tau=1)
  } else {
    df <- df
  }

  df_metrics <- bind_rows(
    tidy(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial())) %>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="No CR") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial()))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial()))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=res_covar),

    tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))%>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="PFS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=covar),

    tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))%>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="OS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=covar) 

  )  %>% mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2)) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2)) %>% mutate(analysis=aname)
  return(df_metrics)
}

# Generate Table Grob
metric_grob_gen <- function(df, metricsi, oh, sz) {
  if(oh=="OR"){
    ohtext = "Adj. OR (95% CI)"
  } else {ohtext ="Adj. HR (95% CI)"}

  dfi <- data.frame(cluster=c("Non-Infl.", "Infl."),
                    Adj_OR = c("Reference", paste0(metricsi$expEstimate, " (", metricsi$low_ci, " - ", metricsi$high_ci, ")")),
                    pvalue= c("", cpval(metricsi$pvalue))) %>%
    tableGrob(rows=NULL, cols=c("", ohtext, "p value"), theme=ttheme_minimal(base_size=sz,
                                                                             core = list(padding=unit(c(2, 2), "mm"))))
  gi <- gtable_add_grob(dfi,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(dfi), l = 1, r = ncol(dfi)) %>%
    gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 1, r = ncol(dfi)) %>%
    gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                    t = 1, b=3, l = 2, r = ncol(dfi))
}

cpval <- function(x) {
  y <- as.character(x)
  if(x >  0.1){y <- "> 0.1"}
  if(x <  0.01){y <- "< 0.01"}
  if(x <  0.001){y <- "< 0.001"}
  return(y)
}

cluster_labs = c(
  "albumin",
  "alk",
  "ast",
  "ferritin",
  "hb",
  "ldh",
  "plt",
  "tbr",
  "il10",
  "il6",
  "tnfa",
  "crp",
  "ddimer",
  "wbc"
)

focus_labs <- paste0("d0_", cluster_labs)[!c(paste0("d0_", cluster_labs) %in% c(
  "d0_albumin",
  "d0_alk",
  "d0_ast",
  "d0_hb",
  "d0_ldh",
  "d0_crp",
  ""
))]

preld_trans_covar <-  "+age + bin_preld_ldh + dx_simple + costim + primary_ref"
trans_covar <- "+age + bridge_categories + bin_preld_ldh + dx_simple + costim + primary_ref"
lbcl_covar="+age + costim + bin_preld_ldh + primary_ref"
lymphoma_covar="+ age + bin_preld_ldh + dx_simple + primary_ref + costim"
costim_covar="+ age + primary_ref + bin_preld_ldh  + car_t_product + dx_simple"
lymphoma_dxres = "everCR_100"
