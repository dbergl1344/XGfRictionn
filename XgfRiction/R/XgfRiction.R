library(roxygen2)
library(devtools)
library(REddyProc)
library(caret)
library(Metrics)
library(stats)
library(tidymodels)
library(ggplot2)
library(parsnip)
library(tune)
library(base)

#' This class represents a processor for XGFriction data.
#' XGFrictionProcessor Class
#'
#' This class represents a processor for XGFriction data.
#'
#' @field site The name of the site associated with the processor.
#' @import REddyProc
#' @export
XGFrictionProcessor <- setRefClass("XGFrictionProcessor",
                                   fields = list(site = "character"),
                                   methods = list(
                                     initialize = function(site) {
                                       site <<- site
                                     },
                                     initialize_xgfriction_processing = function(dataframe, LatDeg, LongDeg, TimeZoneHour) {
                                       require(REddyProc)
                                       # Check if the data is in 30-minute timestamp format
                                       if (any(diff(dataframe$DateTime) != 30*60)) {
                                         mins <- 15*60
                                         print("Adding 15 minutes to DateTime")
                                         print(dataframe$DateTime + mins)
                                         dataframe$DateTime <- (dataframe$DateTime + mins)
                                       }

                                       # Initialize the processing
                                       xgfriction_proc <- REddyProc::sEddyProc$new(site, dataframe$DateTime,
                                                                                   c('NEE', 'Rg', 'Tair', 'VPD', 'Ustar'))
                                       xgfriction_proc$sSetLocationInfo(LatDeg = LatDeg, LongDeg = LongDeg, TimeZoneHour = TimeZoneHour)
                                       return(xgfriction_proc)
                                     }
                                   )
)




#' Initialize XGFriction Processor
#'
#' This function initializes the XGFriction processor object.
#'
#' @param site Site name.
#' @return XGFrictionProcessor object.
#' @export
initialize_xgfriction_processor <- function(site) {
  return(XGFrictionProcessor$new(site))
}


#' Gap Fill MDS Variables
#'
#' This function performs gap filling for specified MDS variables using REddyProc.
#' @param processor XGFriction processing object
#' @param variables Vector of variables to gap fill
#' @param fill_all Logical indicating whether to fill all gaps
#' @return Updated XGFriction processing object
#' @importFrom REddyProc sEddyProc
#' @export
gap_fill_mds_variables <- function(processor, variables, fill_all = FALSE) {
  require(REddyProc)
  # Initialize REddyProc object for gap filling
  reddy_proc <- REddyProc::sEddyProc$new()

  # Loop through variables
  for (variable in variables) {
    # Perform gap filling using REddyProc's sMDSGapFill method
    gap_filled_data <- reddy_proc$sMDSGapFill(processor$data[[variable]], FillAll = fill_all)

    # Update processor with gap-filled data
    processor$data[[variable]] <- gap_filled_data
  }

  return(processor)
}





#' Perform IQR Filtering
#'
#' This function performs IQR filtering on specified variables.
#' @param processor XGFriction processing object
#' @param dataframe Data frame containing the variables to filter
#' @param variables Vector of variables to filter
#' @param threshold_multiplier Multiplier for IQR threshold
#' @return Updated XGFriction processing object
#' @importFrom stats IQR
#' @export
perform_iqr_filtering <- function(xgfriction_proc, dataframe, variables, threshold_multiplier = 6) {
  require(stats)
  for (variable in variables) {
    residual <- dataframe[[paste0(variable, "_orig")]] - dataframe[[paste0(variable, "_fall")]]
    IQR_val <- IQR(residual, na.rm = TRUE)
    outlier <- abs(residual) > (IQR_val * threshold_multiplier)
    xgfriction_proc$data[[variable]] <- ifelse(outlier == 0, dataframe[[paste0(variable, "_orig")]], NA)
  }
  return(xgfriction_proc)
}


#' Perform Ustar Threshold Distribution Estimation and Gap Filling
#'
#' This function estimates the u* threshold distribution and performs gap filling.
#' @param processor XGFriction processing object
#' @param dataframe Data frame containing the variable to perform gap filling on
#' @param variable Variable to perform gap filling on
#' @param nSample Number of samples for estimation
#' @param probs Probabilities for estimation
#' @param fill_all Logical indicating whether to fill all gaps
#' @return Updated XGFriction processing object
#' @importFrom base set.seed
#' @importFrom REddyProc usGetAnnualSeasonUStarMap
#' @export
perform_ustar_gap_fill <- function(processor, dataframe, variable, nSample = 1000L, probs = c(0.05, 0.5, 0.95), fill_all = TRUE) {
  require(stats)
  require(REddyProc)

  base::set.seed(2000)
  processor$sEstUstarThresholdDistribution(nSample = nSample, probs = probs)
  uStarTh <- processor$sGetEstimatedUstarThresholdDistribution()
  uStarThAnnual <- usGetAnnualSeasonUStarMap(uStarTh)[-2]
  uStarSuffixes <- colnames(uStarThAnnual)[-1]
  processor$sMDSGapFillAfterUStarDistr(variable, uStarTh = uStarThAnnual, uStarSuffixes = uStarSuffixes, FillAll = fill_all)
  return(processor)
}



#' Train XGBoost model with custom grid
#'
#' This function trains an XGBoost model using the provided data and a custom grid for hyperparameter tuning.
#'
#' @param data The dataset
#' @param formula The formula for the model
#' @param grid Custom grid for hyperparameter tuning
#' @param folds Number of folds for cross-validation
#' @param repeats Number of repeats for cross-validation
#' @param strata Variable for stratified sampling
#' @param n_trees Number of trees
#' @param tree_depth Maximum tree depth
#' @param min_n Minimum observations in terminal nodes
#' @param loss_reduction Minimum loss reduction to make further splits
#' @param sample_size Fraction of the training set to sample
#' @param mtry Number of variables to sample at each split
#' @param learn_rate Learning rate
#' @return Trained XGBoost model
#' @importFrom parsnip boost_tree
#' @importFrom tune control_grid finalize_model select_best tune_grid
#' @importFrom rsample vfold_cv
#' @importFrom workflows workflow
#' @importFrom doParallel registerDoParallel
#' @import xgboost
#' @export
train_xgboost_model <- function(data, formula, grid, folds = 5, repeats = 2, strata = NULL, n_trees = 2000,
                                tree_depth = tune(), min_n = tune(), loss_reduction = tune(),
                                sample_size = sample_prop(), mtry = tune(), learn_rate = tune()) {
  require(parsnip)
  require(tune)
  require(rsample)
  require(workflows)
  require(doParallel)
  require(xgboost)
  require(base)

  # Prepare data
  base::set.seed(60723)
  split <- rsample::initial_split(data)
  train_data <- rsample::training(split)

  # Set up model specification
  xgb_spec <- parsnip::boost_tree(trees = n_trees,
                         tree_depth = tree_depth,
                         min_n = min_n,
                         loss_reduction = loss_reduction,
                         sample_size = sample_size,
                         mtry = mtry,
                         learn_rate = learn_rate) %>%
    parsnip::set_engine("xgboost") %>%
    parsnip::set_mode("regression")

  # Create workflow
  xgb_workflow <- workflows::workflow() %>%
    workflows::add_formula(formula) %>%
    workflows::add_model(xgb_spec)

  # Create folds
  base::set.seed(60723)
  folds <- rsample::vfold_cv(train_data, v = folds, repeats = repeats, strata = strata)

  # Tune grid
  doParallel::registerDoParallel()
  base::set.seed(60723)
  res <- tune::tune_grid(xgb_workflow, resamples = folds, grid = grid, control = control_grid(save_pred = TRUE))

  # Train model with best parameters
  best_params <- tune::select_best(res, "rmse")
  trained_model <- tune::finalize_model(xgb_spec, best_params)

  return(trained_model)
}




