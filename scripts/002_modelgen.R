library(tidyverse)
library(mclust)
rm(list=ls())
load("output/scaledata.RData")
source("scripts/000_functions_constants.R")

# Preprocess and Generate Model from Development Cohort
df_dev <-  left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>% # Must have curated metadata
  filter(cohort=="MSK Development") %>%
  select(record_id, starts_with("d0_")) %>%
  select(record_id, contains(cluster_labs)) %>%
  drop_na() %>% column_to_rownames("record_id")

# Define model parameters
dev_model <- "VVV"
dev_modelG <- 2 # 2 clusters
d0clusters <- Mclust(df_dev %>% select(starts_with("d0_")), modelNames = dev_model, G = dev_modelG)

entropy = 1 + sum((d0clusters$z+1e-8)*log(d0clusters$z+1e-8))/(nrow(df_dev)*log(2))

# Define the mean laboratory vector for the inflammatory cluster (the one that has the higher IL-6 value)
ref_inflam_vec <- d0clusters$parameters$mean[,which.max(d0clusters$parameters$mean["d0_il6",])]
unit_ref_inflam_vec <- ref_inflam_vec / sqrt(sum(ref_inflam_vec^2))

# Define the model paramaters
dev_mu <- d0clusters$parameters$mean
dev_sig <- d0clusters$parameters$variance$sigma
dev_pro <- d0clusters$parameters$pro
dev_inflammClus <- which.max(d0clusters$parameters$mean["d0_il6",])
dev_LEASTinflammClus <- which.min(d0clusters$parameters$mean["d0_il6",])

save(d0clusters, dev_mu, dev_sig, dev_pro, dev_inflammClus, dev_LEASTinflammClus, file="output/model.RData")

dr <- MclustDR(d0clusters, lambda = 0.9)
plot(dr, what = "evalues")
plot(dr, what = "contour")
plot(dr, what = "boundaries", ngrid = 200)


suppfig5 <- data.frame(d1=dr$dir[,1], d2=dr$dir[,2], class=dr$classification) %>%
  mutate(class=ifelse(class==1, "Non-Inflammatory", "Inflammatory")) %>%
  ggplot(aes(x=d1, y=d2, color=class, fill=class)) + geom_point(size=2, alpha=0.7) + 
  theme_bw() + 
  xlab("GMMDR Direction 1") + 
  ylab("GMMDR Direction 2") + 
  scale_color_manual(values = c("#FF7F0EFF", "#1F77B4FF"))+
  guides(color=FALSE, fill=FALSE)

# Extended Figure 2h
feats <- c(3, 12)
d0clusters1 <- d0clusters
d0clusters1$parameters$mean <- d0clusters1$parameters$mean[feats,]
d0clusters1$parameters$variance$sigma <- d0clusters1$parameters$variance$sigma[feats,feats,]
d0clusters1$parameters$variance$shape <- d0clusters1$parameters$variance$shape[feats,]
d0clusters1$data <- d0clusters1$data[,feats]
library(factoextra)

efig2h <- fviz_mclust(d0clusters1, "classification", geom = "point",
            ylab="Log-Scaled C-Reactive Protein (CRP)", xlab="Scaled Aspartate Transaminase (AST)",
            ggtheme=theme_minimal) +
  labs(
       title = "",
       subtitle="INFLAMIX Clustering",
       color="Cluster"
  ) +
  scale_fill_manual(values=colorclusters(2))+
scale_color_manual(values=colorclusters(2), labels=c('Non-Inflammatory', 'Inflammatory'))+
  guides(fill=FALSE, shape=FALSE) + xlim(-1, 2.5)



# Extended Figure 2i
cluster_labs = c( # Cluster with or without considering WBC, Plt, Tbili (tbr)
    "albumin",
    "alk",
    "ast",
    "ferritin",
    "hb",
    "ldh",
    #"plt",
    #"tbr",
    "il10",
    "il6",
    "tnfa",
    "crp",
    "ddimer"
    #"wbc"
  )

df_dev <-  left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>% # Must have curated metadata
  filter(cohort=="MSK Development") %>%
  select(record_id, starts_with("d0_")) %>%
  select(record_id, contains(cluster_labs)) %>%
  drop_na() %>% column_to_rownames("record_id")

# Define model parameters
feats <- c(3, 10)
dev_model <- "VVV"
dev_modelG <- 2 # 2 clusters
d0clusters <- Mclust(df_dev %>% select(starts_with("d0_")), modelNames = dev_model, G = dev_modelG)

d0clusters1 <- d0clusters
d0clusters1$parameters$mean <- d0clusters1$parameters$mean[feats,]
d0clusters1$parameters$variance$sigma <- d0clusters1$parameters$variance$sigma[feats,feats,]
d0clusters1$parameters$variance$shape <- d0clusters1$parameters$variance$shape[feats,]
d0clusters1$data <- d0clusters1$data[,feats]

efig2hi <- fviz_mclust(d0clusters1, "classification", geom = "point",
                      ylab="Log-Scaled C-Reactive Protein (CRP)", xlab="Scaled Aspartate Transaminase (AST)",
                      ggtheme=theme_minimal) +
  labs(
    title = "",
    subtitle="Model derived without WBC, Plt, Tbili",
    color="Cluster"
  ) +
  scale_fill_manual(values=colorclusters(2))+
  scale_color_manual(values=colorclusters(2), labels=c('Non-Inflammatory', 'Inflammatory'))+
  guides(fill=FALSE, shape=FALSE) + xlim(-1, 2.5)

