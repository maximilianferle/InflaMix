---
title: "InflaMix"
author: "Sandeep Raj"
date: "2025-03-31"
---

## InflaMix

This github repository contains code for deriving InflaMix as described in https://www.nature.com/articles/s41591-025-03532-x 

InflaMix (INFLAmmation MIXture Model) is a Gaussian mixture model trained on pre-infusion laboratory and cytokine data from patients with large B-cell lymphoma treated with CD19-CAR-T. It defines two clusters, the inflammatory and non-inflammatory clusters. The inflammatory cluster is enriched for patients with a serological signature that is consistent with systemic inflammation (e.g., elevated CRP, elevated IL6, low albumin).

The model can be applied to estimate the probability of cluster assignment even with several missing laboratory or cytokine values. Enrichment for the inflammatory signature reproducibly stratifies risk of treatment failure in multiple non-Hodgkin lymphoma cohorts and is predictive of CD19-CAR-T treatment failure in large B-cell lymphoma beyond standard biomarkers and risk factors. 

For additional details, please email rajs@mskcc.org

Terms of Use as described in the LICENSE file in this repository - please review for details. Briefly, InflaMix is intended for educational and research use only and is not yet intended for clinical use or as a medical device. 

## Online Calculator
The link for the InflaMix online calculator is provided here: https://ssraj017.github.io/inflamix_app_prod/

Legacy Link: https://sraj17.shinyapps.io/inflamix_app/

For optimal results, users enter as many of the 14 labs as are available, giving precedence to the 6 laboratory measurements from the limited panel.

## Set-up for analyses included in the InflaMix paper

1. Download this entire directory from GitHub. Open the "InflaMix.Rproj" file in a suitable IDE like R-Studio.

2. Install all required packages manually by reviewing required libraries at the top of every script in the "scripts/" directory
 - Note that the "survcomp" and "ComplexHeatmap" packages require installation through BioConductor.
  - https://bioconductor.org/packages/release/bioc/html/survcomp.html
  - https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html 
 - Download and Installation should not more than a few minutes to an hour. 

3. Follow the steps below.

A description of the data inputs for this project and a description of their column values is given in the comments at the beginning of 001_base_scaling.R. 

Scripts should be run in the following order which contains code for the corresponding figures:

### Script Descriptions

The following 

0. Please place Dataset 1 (deriv_cohort_d0_labs_v2_dataset1.csv) into the data folder. This was provided as source data with the paper. 

Every script here should run within a minute. There are two exceptions (see steps 5 and 7)

1. 001_base_scaling.R 
 - Supplementary Figures 2 and 7

2. 002_modelgen.R 
 - Extended Figure 3h and 3i, Supplementary Figure 8
  
3. 003_icl_bic.R
 - Extended Figure 3a
  
4. 004_cluster_lab_distributions.R
 - Extended Figure 3b-g
  
5. 005_deriv_cohort_properties.R 
  - Figure 1a-g, Supplementary Figure 3
  - Note that Figures 1d, e, f cannot be plotted with provided data, code should be run line-by-line.
  - In generating Figure 1g, the run time for generating each independent iteration of random forest will take 3-5 minutes. For 100 iterations this will take several hours on a local environment. 
  
6. 006_deriv_outcomes.R
 - Figure 1h-l
 - Note that this script cannot be run with provided data. 
 
7. 007_partial_lab_clustering.R
 - Supplementary Figure 4
 - Note that this script cannot be run with provided data. Though the code and comments describe the analysis to evaluate the quality of clustering with partially available laboratory data 
 - The run time for this script will take ~2 hours on a local environment. 
 
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
 - Figure 3 and Extended Figure 7
 - Note that this script cannot be run with provided data. 
 
 12. 012_cluster_transitions.R
 - Figure 5, Supplementary Figure 5, and Supplementary Figure 6
 - Note that this script cannot be run with provided data. 
 
 13. 100_tables.R
  - Table 1, Extended Tables 1 and 2
  - Note that this script cannot be run with provided data. 

All required libraries are listed within each script. 

The InflaMix model features and parameters will be in the "model" directory saved as an R data object. 

## Data availability

Data requests for patient-related laboratory measurements or clinical outcomes will be reviewed by the corresponding author in consultation with coauthors from Hackensack University Medical Center and Sheba Medical Center. Any data and materials that can be shared will be released via data transfer agreement. Laboratory values and their corresponding upper limits of normal for the model-derivation cohort of MSK patients are provided in the Source Data provided with the paper (Data Set 1). With this data, the following figures can be reproduced:
 - Figure 1 a-c, g
 - Extended Figure 3a-i
 - Extended Figure 5a
 - Supplementary Figure 2
 - Supplementary Figure 3
 - Supplementary Figure 7
 - Supplementary Figure 8
