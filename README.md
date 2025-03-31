---
title: "InflaMix"
author: "Sandeep Raj"
date: "2025-02-25"
---

## InflaMix Manuscript Code

This github repository contains code and models that accompany the InflaMix manuscript while it is under review. 

After publication, this git will be updated with:
 - a link to the paper
 - guidance on how to use InflaMix
 - additional code to further facilitate direct implementation of InflaMix for research with new data sets

A description of the two data inputs for this project and a description of their column values is given in the comments at the beginning of 001_base_scaling.R. 

Scripts should be run in the following order which contains code for the corresponding figures:

For additional details, please email rajs@mskcc.org

## Online Calculator
The link for the InflaMix online calculator is provided here: https://ssraj017.github.io/inflamix_app_prod/

Legacy Link: https://sraj17.shinyapps.io/inflamix_app/

For optimal results, users enter as many of the 14 labs as are available, giving precedence to the 6 laboratory measurements from the limited panel.

## Set-up

1. Download this entire directory from GitHub. Open the "InflaMix.Rproj" file in a suitable IDE like R-Studio.

2. Install all required packages, either manually by reviewing required libraries at the top of every script in the "scripts/" directory, or use the "renv" package "https://rstudio.github.io/renv/articles/renv.html"
 - Note that the "survcomp" and "ComplexHeatmap" packages require installation through BioConductor.
  - https://bioconductor.org/packages/release/bioc/html/survcomp.html
  - https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html 
 - Download and Installation should not more than a few minutes to an hour. 

3. Follow the steps below. 

### 0 . Please place Dataset 1 (deriv_cohort_d0_labs_v2_dataset1.csv) into the data folder. This was provided as a supplementary material. 

Every script here should run within a minute. There are two exceptions (see steps 5 and 7)

1. 001_base_scaling.R 
 - Supplementary Figures 2 and 8

2. 002_modelgen.R 
 - Extended Figure 3h and 3i, Supplementary Figure 9
  
3. 003_icl_bic.R
 - Extended Figure 3a
  
4. 004_cluster_lab_distributions.R
 - Extended Figure 3b-g
  
5. 005_deriv_cohort_properties.R 
  - Figure 1a-g, Supplementary Figure 3
  - Note that Figures 1d, e, f cannot be plotted with provided data, code should be run line-by-line.
  - In generating Figure 1g, the run time for generating each independent iteration of random forest will take 3-5 minutes. For 100 iterations this will take several hours. 
  
6. 006_deriv_outcomes.R
 - Figure 1h-l
 - Note that this script cannot be run with provided data. 
 
7. 007_partial_lab_clustering.R
 - Supplementary Figure 4
 - Note that this script cannot be run with provided data. Though the code and comments describe the analysis to evaluate the quality of clustering with partially available laboratory data 
 - The run time for this script will take ~2 hours. 
 
8. 008_cluster_assignment_validation.R
 - Extended Figure 4
 - Note that this script cannot be run with provided data. 
 
9. 009_outcome_association_validation.R
 - Figure 2, Figure 4, Extended Figure 6, and Extended Figure 7
 - Note that this script cannot be run with provided data. 
 
10. 010_cohort_heatmaps.R
 - Extended Figure 5
 - Note that only Extended Figure 5a be plotted with provided data

11. 011_prediction_benchmarking.R
 - Figure 3 and Supplementary Figure 5
 - Note that this script cannot be run with provided data. 
 
 12. 012_cluster_transitions.R
 - Figure 5, Supplementary Figure 6, and Supplementary Figure 7
 - Note that this script cannot be run with provided data. 
 
 13. 100_tables.R
  - Table 1, Supplementary Table 3
  - Note that this script cannot be run with provided data. 

All required libraries are listed within each script. 

The InflaMix model features are in the "model" directory. 

## Data availability

Data requests for patient-related laboratory measurements or clinical outcomes will be reviewed by the corresponding author in consultation with coauthors from Hackensack University Medical Center and Sheba Medical Center. Any data and materials that can be shared will be released via data transfer agreement. Laboratory values and their corresponding upper limits of normal for the model-derivation cohort of MSK patients are provided in the Supplementary Data (Data Set 1). With this data, the following figures can be reproduced:
 - Figure 1 a-c, g
 - Extended Figure 3a-i
 - Extended Figure 5a
 - Supplementary Figure 2
 - Supplementary Figure 3
 - Supplementary Figure 8
 - Supplementary Figure 9
