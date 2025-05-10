library(mgcv)
library(dplyr)
library(stringr)
library(gratia)

setwd("C:\\HCP_processing\\HCP_significance")


# Loops that check interaction between BIS and other covariates like age and gender

datasets <- list(
  relational_relational = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\relational_relational_cleaned.csv"),
  WM_Body_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Body_2bk_cleaned.csv"),
  WM_Face_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Face_2bk_cleaned.csv"),
  WM_Place_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Place_2bk_cleaned.csv"),
  WM_Tool_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Tool_2bk_cleaned.csv"),
  math = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\math_cleaned.csv"),
  story = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\story_cleaned.csv"),
  ET_Shape = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\ET_Shape_cleaned.csv"),
  ET_Face = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\ET_Face_cleaned.csv")
)

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Initialize results data frame
  results <- data.frame()
  
  # Initialize p-value vectors
  gam.int.F <- c()
  gam.int.pvalue <- c()
  
  # Loop through each ROI brain region
  for (roi in roi_columns) {
    # Run gam.fit.smooth function
    result <- gam.factorsmooth.interaction(dataset, region = roi, 
                             smooth_var = smooth_var, int_var = int_var, covariates = covariates, 
                             knots = 3, set_fx = FALSE)
    
    # Add results to the results data frame
    results <- rbind(results, result)
    
    # Collect F-values and p-values
    gam.int.F <- c(gam.int.F, result[1, "gam.int.F"])
    gam.int.pvalue <- c(gam.int.pvalue, result[1, "gam.int.pvalue"])
    
  }
  
  # Calculate FDR-adjusted p-values
  fdr_adjusted_pvalues <- p.adjust(gam.int.pvalue, method = 'fdr')
  
  # Add FDR-adjusted p-values to the results data frame
  results$fdr_adjusted_pvalue <- fdr_adjusted_pvalues
  
  # Save results to a CSV file
  write.csv(results, file = paste0(dataset_name, "_BIS_age_results.csv"), row.names = FALSE)
}



# Loops that check the significance of general additive model significance

setwd("C:\\HCP_processing\\HCP_significance\\smooth_var\\control_for_age_gender")

roi_columns = colnames(WM_Body_2bk)[10:369] # Get the column names of ROI


covariates = c('Gender', 'Age')
smooth_var = 'BIS'
knots = 3


datasets <- list(
  WM_Face_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Face_2bk_cleaned.csv"),
  WM_Place_2bk = read.csv("C:\\HCP_processing\\BIS_outlier_removed\\WM_Place_2bk_cleaned.csv")
)

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  
  # Initialize results data frame
  results <- data.frame()
  
  # Initialize p-value vectors
  var_p_values <- c()
  anova_p_values <- c()
  
  # Loop through each ROI brain region
  for (roi in roi_columns) {
    # Run gam.fit.smooth function
    result <- gam.fit.smooth(dataset, region = roi, 
                             smooth_var = smooth_var, covariates = covariates, 
                             knots = 3, set_fx = FALSE, stats_only = TRUE)
    
    # Add results to the results data frame
    results <- rbind(results, result)
    
    # Collect p-values
    var_p_values <- c(var_p_values, result[1, "gam.smooth.pvalue"])
    anova_p_values <- c(anova_p_values, result[1, "anova.smooth.pvalue"])
    
  }
  
  # Calculate FDR-adjusted p-values
  var_fdr_adjusted_pvalues <- p.adjust(var_p_values, method = "fdr")
  anova_fdr_adjusted_pvalues <- p.adjust(anova_p_values, method = 'fdr')
  
  # Add FDR-adjusted p-values to the results data frame
  results$var_fdr_adjusted_pvalue <- var_fdr_adjusted_pvalues
  results$anova_fdr_adjusted_pvalue <- anova_fdr_adjusted_pvalues
  
  # Save results to a CSV file
  write.csv(results, file = paste0(dataset_name, "_gam_results.csv"), row.names = FALSE)
}
