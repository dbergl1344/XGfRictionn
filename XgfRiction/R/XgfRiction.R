library(roxygen2)
library(devtools)
library(REddyProc)
library(Metrics)
library(stats)
library(tidymodels)
library(ggplot2)
library(parsnip)
library(tune)
library(base)
library(rsample)
library(caret)
library(xgboost)


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
#' @param seed Seed value for random number generation
#' @return Updated XGFriction processing object
#' @importFrom REddyProc usGetAnnualSeasonUStarMap
#' @export
perform_ustar_gap_fill <- function(processor, dataframe, variable, nSample = 1000L, probs = c(0.05, 0.5, 0.95), fill_all = TRUE, seed = NULL) {
  require(stats)
  require(REddyProc)

  if (!is.null(seed)) {
    set.seed(seed)
  }

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
#' @param seed_count Number of seed values to use
#' @param seed_for_data Seed value for data preparation
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
                                sample_size = sample_prop(), mtry = tune(), learn_rate = tune(),
                                seed_count = 1, seed_for_data = 60723) {
  require(parsnip)
  require(tune)
  require(rsample)
  require(workflows)
  require(doParallel)
  require(xgboost)
  require(base)

  # Prepare data
  set.seed(seed_for_data)
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
  seeds <- sample.int(1e6, seed_count)
  folds <- lapply(seeds, function(seed) {
    set.seed(seed)
    rsample::vfold_cv(train_data, v = folds, repeats = repeats, strata = strata)
  })

  # Tune grid
  doParallel::registerDoParallel()
  res <- tune::tune_grid(xgb_workflow, resamples = folds, grid = grid, control = control_grid(save_pred = TRUE))

  # Train model with best parameters
  best_params <- tune::select_best(res, "rmse")
  trained_model <- tune::finalize_model(xgb_spec, best_params)

  return(trained_model)
}





#' Select Best Model
#'
#' This function selects the best model based on specified evaluation metrics.
#'
#' @param results A data frame containing the results of model training (e.g., cross-validation results).
#' @param metric The evaluation metric used for selecting the best model (e.g., "rmse", "rsq").
#' @param higher_is_better Logical indicating whether a higher value of the evaluation metric is better.
#' @return The best model based on the specified evaluation metric.
#' @importFrom tune select_best
#' @importFrom tune select_best
#' @export
select_best_model <- function(results, metric, higher_is_better = TRUE) {
  # Add code here to select the best model based on the specified evaluation metric
}



#' Plot Variables of Importance
#'
#' This function plots the variables of importance of a trained XGBoost model.
#'
#' @param model The trained XGBoost model.
#' @param n_top Number of top variables to plot.
#' @return A plot showing the variables of importance.
#' @import ggplot2
#' @export
plot_variables_of_importance <- function(model, n_top = 10) {
  ggplot2::ggplot(ggplot2::aes(x = Importance, y = reorder(Variable, Importance))) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Top Variables of Importance",
                  x = "Importance",
                  y = "Variable")
}




#' Plot Model Metrics
#'
#' This function plots the metrics of a trained XGBoost model.
#'
#' @param metrics The metrics data frame.
#' @param metric_name The name of the metric to plot.
#' @return A plot showing the specified metric.
#' @import ggplot2
#' @export
plot_model_metrics <- function(metrics, metric_name) {
  ggplot2::ggplot(metrics, ggplot2::aes(x = .pred, y = .resid)) +
    ggplot2::geom_point(color = "purple") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "yellow") +
    ggplot2::labs(x = "Predicted NEE",
                  y = "Residuals",
                  title = "Residuals vs. Predicted Values (Test Set)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}



#' Use trained XGBoost model to predict missing data
#'
#' This function uses a trained XGBoost model to predict missing data in a dataset.
#'
#' @param model The trained XGBoost model.
#' @param data The dataset with missing values.
#' @return The dataset with missing values filled using the XGBoost model.
#' @import xgboost
#' @export
predict_missing_data <- function(model, data) {
  require(xgboost)

  # Predict missing values
  predicted_values <- predict(model, newdata = data)

  # Replace NA values with predicted values
  data[is.na(data)] <- predicted_values[is.na(data)]

  return(data)
}



#' Create Gap-Filled Column
#'
#' This function creates a new column in the dataset where missing values are filled with predictions from a specified model.
#'
#' @param data The dataset.
#' @param varname The name of the target variable.
#' @param model The trained model for predicting missing values.
#' @return The dataset with the new gap-filled column.
#' @import xgboost
#' @export
create_gap_filled_column <- function(data, varname, model) {
  require(xgboost)

  # Create new column name
  new_col_name <- paste0(varname, "_f")

  # Copy original values to the new column
  data[[new_col_name]] <- data[[varname]]

  # Find rows with missing values
  missing_rows <- which(is.na(data[[varname]]))

  # If there are missing values, fill them with predictions from the model
  if (length(missing_rows) > 0) {
    # Predict missing values
    predicted_values <- predict(model, newdata = data[missing_rows, ])

    # Replace NA values with predicted values in the new column
    data[missing_rows, new_col_name] <- predicted_values
  }

  return(data)
}




# Assuming 'data' is your dataset and 'my_model' is your trained XGBoost model
# data <- create_gap_filled_column(data, "varname", my_model)




#' Replace Variables in REddyProc Object
#'
#' This function replaces specified variables in a REddyProc object with columns from a dataset.
#'
#' @param processor The REddyProc object.
#' @param data The dataset containing columns to replace variables with.
#' @param varnames Vector of variable names to replace in the REddyProc object.
#' @import data.table
#' @export
replace_variables_in_REddyProc <- function(processor, data, varnames) {
  require(data.table)

  # Convert dataset to data.table for easier column manipulation
  dt_filled <- as.data.table(data)

  # Loop through each variable to replace
  for (varname in varnames) {
    # Get the corresponding column from the dataset
    col_name <- paste0(varname, "_f")

    # Check if the column exists in the dataset
    if (col_name %in% names(dt_filled)) {
      # Replace the variable in the REddyProc object with the corresponding column
      processor$sTEMP[[varname]] <- dt_filled[[col_name]]
    } else {
      warning(paste("Column", col_name, "not found in the dataset. Skipping replacement."))
    }
  }

  return(processor)
}






#' Perform Flux Partitioning
#'
#' This function performs flux partitioning using the specified method and parameters.
#'
#' @param processor The REddyProc object.
#' @param method The partitioning method to use. Options: "sGL", "sMR", "sTK".
#' @param params List of parameters required for the chosen method.
#'   - For "sGL" method:
#'     - useLocaltime: Logical. If TRUE, use local time zone instead of geo-solar time to compute potential radiation.
#'     - debug.l: Deprecated. Use debug instead.
#'     - isWarnReplaceColumns: Logical. Set to FALSE to avoid the warning on replacing output columns.
#'   - For "sMR" method:
#'     - FluxVar: Character. Name of the net ecosystem flux variable.
#'     - QFFluxVar: Character. Name of the quality-flagged net ecosystem flux variable.
#'     - QFFluxValue: Numeric. Quality flag value to indicate bad flux values.
#'     - TempVar: Character. Name of the temperature variable.
#'     - QFTempVar: Character. Name of the quality-flagged temperature variable.
#'     - QFTempValue: Numeric. Quality flag value to indicate bad temperature values.
#'     - RadVar: Character. Name of the radiation variable.
#'     - TRef: Numeric. Reference temperature.
#'     - Suffix.s: Character. Suffix for the partitioning results.
#'   - For "sTK" method:
#'     - controlGLPart: Further default parameters, such as suffix.
#' @export
perform_flux_partitioning <- function(processor, method, params) {
  method <- tolower(method)
  if (method == "sgl") {
    if (!is.null(params$uStarScenKeep)) {
      params <- c(params, isWarnReplaceColumns = FALSE)
    }
    results <- processor$sGLFluxPartition(params)
  } else if (method == "smr") {
    results <- processor$sMRFluxPartition(params)
  } else if (method == "stk") {
    results <- processor$sTKFluxPartition(params)
  } else {
    stop("Invalid partitioning method. Choose from 'sGL', 'sMR', or 'sTK'.")
  }
  return(results)
}





###check

# # Example usage with different processor names
# # Assuming 'my_processor' is the user's REddyProc object
# nighttime_params <- list(FluxVar = "NEE_f", QFFluxVar = "NEE_fqc", QFFluxValue = 0,
#                          TempVar = "Tair_f", QFTempVar = "Tair_fqc", QFTempValue = 0,
#                          RadVar = "Rg", TRef = 273.15 + 15, Suffix.s = "suffix")
# nighttime_results <- perform_flux_partitioning(my_processor, "sMR", nighttime_params)



