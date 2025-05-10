# set working directory and load libraries
setwd("~/Desktop/hcp_data/model_comparison")
library(mgcv)
library(dplyr)
library(tidyr)

# Read data
face  <- read.csv("~/Desktop/hcp_data/WM_Face_2bk_cleaned.csv")
place <- read.csv("~/Desktop/hcp_data/WM_Place_2bk_cleaned.csv")

# function to perform k & m selection for one dataset with progress bar
run_km_selection <- function(dat, k.values, m.values, covariates=c("Gender","Age")) {
  # identify ROI columns by suffix _ROI
  roi.cols <- grep("_ROI$", names(dat), value=TRUE)
  n.rois <- length(roi.cols)
  total  <- n.rois * length(k.values) * length(m.values)
  
  # prepare storage
  results <- data.frame()
  
  # setup progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  counter <- 0
  
  # loop through ROIs
  for (roi in roi.cols) {
    y <- dat[[roi]]
    if (all(is.na(y))) {
      counter <- counter + length(k.values) * length(m.values)
      setTxtProgressBar(pb, counter)
      next
    }
    
    # test each (k, m) pair
    for (K in k.values) {
      for (M in m.values) {
        form <- as.formula(paste(roi, "~ s(BIS, k=", K, ", m=", M,
                                 ") +", paste(covariates, collapse=" + ")))
        fit <- try(gam(form, data=dat, method="REML"), silent=TRUE)
        aic_val <- if (inherits(fit, "try-error")) NA else AIC(fit)
        results <- rbind(results, data.frame(ROI=roi, k=K, m=M, AIC=aic_val))
        
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
      }
    }
  }
  close(pb)
  
  # Find best (k, m) per ROI
  best_km <- results %>%
    group_by(ROI) %>%
    filter(!is.na(AIC)) %>%
    slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Count how many times each (k, m) combination wins
  win.count <- best_km %>%
    count(k, m, name = "wins") %>%
    arrange(k, m)
  
  return(list(results = results, wins = win.count))
}

# define k and m values to test
k.grid <- c(3,4,5,6)
m.grid <- c(4,5,6)

# run for Face task with progress
cat("Running k&m selection on Face data:\n")
eres_face <- run_km_selection(face, k.grid, m.grid)

# run for Place task with progress
cat("\nRunning k&m selection on Place data:\n")
eres_place <- run_km_selection(place, k.grid, m.grid)

# summarize results
wins_face  <- eres_face$wins  %>% rename(wins_face=wins)
wins_place <- eres_place$wins %>% rename(wins_place=wins)

wins_all <- full_join(wins_face, wins_place, by=c("k", "m")) %>%
  replace_na(list(wins_face=0, wins_place=0))

# display and save
print(wins_all)
write.csv(wins_all, "km_selection_wins.csv", row.names=FALSE)
