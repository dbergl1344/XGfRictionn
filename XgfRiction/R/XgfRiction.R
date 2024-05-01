#' Convert to POSIX time format
#'
#' Converts the specified data columns to POSIXct format using the REddyProc package function `fConvertTimeToPosix`.
#'
#' @param data Data frame containing the data.
#' @param date_format Character string specifying the date format (e.g., "YMDHM").
#' @param year_col Name of the year column.
#' @param month_col Name of the month column.
#' @param day_col Name of the day column.
#' @param hour_col Name of the hour column.
#' @param min_col Name of the minute column.
#' @param posix_time_col Name of the POSIXct column to be created.
#' @return Data frame with the POSIXct column.
#' @importFrom REddyProc fConvertTimeToPosix
#' @export
convert_to_POSIX <- function(data, date_format, year_col, month_col, day_col, hour_col, min_col, posix_time_col) {
  # Convert date and time columns to POSIXct format using fConvertTimeToPosix
  data <- REddyProc::fConvertTimeToPosix(data, date_format, Year = year_col, Month = month_col, Day = day_col, Hour = hour_col, Min = min_col)

  # Convert the specified column to POSIXct format
  data[[posix_time_col]] <- as.POSIXct(data[[posix_time_col]], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  return(data)
}




#' Subtract 15 minutes from DateTime column
#'
#' Subtracts 15 minutes (900 seconds) from the DateTime column of a data frame.
#'
#' @param data Data frame containing the data.
#' @param posix_time_col Name of the POSIXct column (DateTime) to adjust.
#' @return Data frame with adjusted DateTime column.
#' @export
subtract_15_minutes_from_datetime <- function(data, posix_time_col) {
  # Calculate the interval in seconds (15 minutes)
  interval <- 15 * 60  # 15 minutes in seconds

  # Subtract 15 minutes (900 seconds) from the DateTime column
  data[[posix_time_col]] <- data[[posix_time_col]] - interval

  # Return the data frame with the adjusted DateTime column
  return(data)
}


#' Initialize user_class processor
#'
#' Initializes the user_class processor object and sets location information.
#'
#' @param site_ID Site ID from the AmeriFlux website.
#' @param data Data frame containing the eddy data.
#' @param column_names Vector with selected column names.
#' @param posix_time_column Column name with POSIX time stamp.
#' @param LatDeg Latitude of the site.
#' @param LongDeg Longitude of the site.
#' @param TimeZoneHour Time zone offset in hours.
#' @return Initialized user_class object with location information set.
#' @importFrom REddyProc sEddyProc
#' @export
initialize_user_class_processor <- function(site_ID, data, column_names, posix_time_column, LatDeg, LongDeg, TimeZoneHour) {
  # Initialize the user_class processor object
  user_class <- REddyProc::sEddyProc$new(site_ID, data, column_names)

  # Set location information for the user_class processor
  user_class$sSetLocationInfo(LatDeg = LatDeg, LongDeg = LongDeg, TimeZoneHour = TimeZoneHour)

  return(user_class)
}


#' Calculate VPD and Tair
#'
#' Calculates VPD and Tair for the specified data frame.
#'
#' @param data Data frame containing the data.
#' @param rh_col Name of the relative humidity column.
#' @param tair_col Name of the air temperature column.
#' @return Data frame with calculated VPD and Tair columns added.
#' @importFrom REddyProc fCalcVPDfromRHandTair
#' @export
calculate_vpd_and_tair <- function(data, rh_col, tair_col) {
  # Calculate VPD from relative humidity and air temperature using fCalcVPDfromRHandTair function
  vpd <- REddyProc::fCalcVPDfromRHandTair(data[[rh_col]], data[[tair_col]])

  # Add the calculated VPD column to the data frame
  data$VPD <- vpd

  return(data)
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



#' Perform gap-filling of meteorological data for XGFrictionProcessor object
#'
#' This function performs gap-filling of meteorological data using the `sMDSGapFill` method for an XGFrictionProcessor object.
#'
#' @param processor XGFrictionProcessor object to perform gap-filling on.
#' @param variables Vector of variable names to gap-fill.
#' @param fill_all Logical, indicating whether to fill all gaps in the variables or only the missing values between two valid observations.
#' @return Updated XGFrictionProcessor object with gap-filled meteorological data.
#' @importFrom REddyProc sEddyProc_sMDSGapFill
#' @export
gap_fill_met_data <- function(processor, variables, fill_all = FALSE) {
  # Perform gap-filling on specified variables using sMDSGapFill method
  for (variable in variables) {
    processor$sEddyProc_sMDSGapFill(variable, FillAll = fill_all)
  }
  # Return the updated processor object
  return(processor)
}





#' Estimate u* threshold distribution for an XGFrictionProcessor object
#'
#' This function estimates the u* threshold distribution using the specified processor.
#'
#' @param processor XGFrictionProcessor object.
#' @param n_sample Number of samples for estimation (default: 700L).
#' @param probs Probabilities for estimation (default: c(0.05, 0.5, 0.95)).
#' @return Data frame containing the estimated u* threshold distribution.
#' @importFrom REddyProc sEddyProc_sEstUstarThresholdDistribution sEddyProc_sGetEstimatedUstarThresholdDistribution
#' @export
estimate_ustar_threshold_distribution <- function(processor, n_sample = 700L, probs = c(0.05, 0.5, 0.95)) {
  processor$sEstUstarThresholdDistribution(nSample = n_sample, probs = probs)
  uStarTh <- processor$sGetEstimatedUstarThresholdDistribution()
  return(uStarTh)
}


#' Filter u* Threshold Distribution by Aggregation Mode
#'
#' This function filters a data frame containing the estimated u* threshold distribution based on the specified aggregation mode.
#'
#' @param ustar_distribution Data frame containing the estimated u* threshold distribution.
#' @param mode Character string specifying the aggregation mode ("year", "season", or "single"). Default is "year".
#' @return Data frame filtered based on the specified aggregation mode with selected quantiles (5%, 50%, 95%).
#' @importFrom dplyr filter select
#' @export
filter_ustar_by_mode <- function(ustar_distribution, mode = "year") {
  filtered_ustar <- ustar_distribution %>%
    dplyr::filter(aggregationMode == mode) %>%
    dplyr::select(`5%`, `50%`, `95%`)
  return(filtered_ustar)
}



#' Get Annual u* Threshold Map
#'
#' This function obtains the annual or seasonal u* threshold map from a filtered data frame.
#'
#' @param filtered_ustar Data frame filtered based on the specified aggregation mode with selected quantiles (5%, 50%, 95%).
#' @return Data frame containing the annual or seasonal u* threshold map.
#' @export
get_ustar_annual <- function(filtered_ustar) {
  # Return the data frame as is since it contains the quantile columns (5%, 50%, 95%)
  uStarThAnnual <- filtered_ustar
  # Print the data frame to verify its structure
  print(uStarThAnnual)
  # Return the data frame
  return(uStarThAnnual)
}



#' Get u* Suffixes from u* Threshold Map
#'
#' This function retrieves the u* suffixes from the u* threshold map.
#'
#' @param uStarTh Data frame containing the u* threshold map.
#' @return Character vector containing the column names of the u* threshold map.
#' @export
get_ustar_suffixes <- function(uStarTh) {
  # Retrieve column names (u* suffixes) of the u* threshold map
  uStarSuffixes <- colnames(uStarTh)
  return(uStarSuffixes)
}




#' Perform Gap Filling using u* Threshold Distribution
#'
#' This function performs gap filling for a specified variable using the u* threshold distribution.
#'
#' @param processor An instance of the sEddyProc class or similar processor object containing the user's data.
#' @param variable Character string specifying the variable to gap fill (e.g., "NEE").
#' @param ustar_distribution Data frame containing the estimated u* threshold distribution. It should include the quantiles (5%, 50%, 95%) for the specified aggregation mode (e.g., "year").
#' @param fill_all Logical indicating whether to fill all gaps (default: TRUE). If `TRUE`, all gaps in the specified variable will be filled using the u* threshold distribution. If `FALSE`, only gaps matching specific conditions will be filled.
#' @return Updated processor object after gap filling. The processor object will have the specified variable gap-filled using the u* threshold distribution.
#' @importFrom REddyProc sEddyProc_sMDSGapFillAfterUStarDistr
#' @export
perform_gap_filling_with_ustar <- function(processor, variable, ustar_distribution, fill_all = TRUE) {
  # Step 1: Filter the u* threshold distribution for yearly aggregation
  filtered_ustar <- filter_ustar_by_mode(ustar_distribution, mode = "year")

  # Step 2: Get the annual u* threshold map
  uStarThAnnual <- get_ustar_annual(filtered_ustar)

  # Step 3: Get the u* suffixes
  uStarSuffixes <- get_ustar_suffixes(uStarThAnnual)

  # Step 4: Perform gap filling using the u* threshold distribution
  processor$sEddyProc_sMDSGapFillAfterUStarDistr(
    variable,
    uStarTh = uStarThAnnual,
    uStarSuffixes = uStarSuffixes,
    FillAll = fill_all
  )

  # Return the updated processor object
  return(processor)
}




#' Export results from XGFrictionProcessor object
#'
#' This function exports the results from the specified processor after gap filling.
#'
#' @param processor XGFrictionProcessor object.
#' @return Data frame containing the exported results.
#' @importFrom REddyProc sEddyProc_sExportResults
#' @export
export_results_from_processor <- function(processor) {
  results <- processor$sEddyProc_sExportResults()
  return(results)
}






#' Prepare data for model training
#'
#' This function splits the data into training and testing sets.
#'
#' @param data The dataset to be split
#' @return A list containing the training and testing sets
#' @importFrom rsample initial_split training testing
#' @export
prepare_data <- function(data) {
  split <- rsample::initial_split(data)
  train_data <- rsample::training(split)
  test_data <- rsample::testing(split)
  list(train_data = train_data, test_data = test_data)
}

#' Create model specification
#'
#' This function creates a model specification for the XGBoost model.
#'
#' @param n_trees Number of trees for the model
#' @param tree_depth Maximum tree depth
#' @param min_n Minimum observations in terminal nodes
#' @param loss_reduction Minimum loss reduction to make further splits
#' @param sample_size Fraction of the training set to sample
#' @param mtry Number of variables to sample at each split
#' @param learn_rate Learning rate
#' @return A model specification for XGBoost
#' @importFrom parsnip boost_tree set_engine set_mode
#' @export
create_model_spec <- function(n_trees = 2000, tree_depth = tune(), min_n = tune(), loss_reduction = tune(),
                              sample_size = sample_prop(), mtry = tune(), learn_rate = tune()) {
  parsnip::boost_tree(trees = n_trees,
                      tree_depth = tree_depth,
                      min_n = min_n,
                      loss_reduction = loss_reduction,
                      sample_size = sample_size,
                      mtry = mtry,
                      learn_rate = learn_rate) %>%
    parsnip::set_engine("xgboost") %>%
    parsnip::set_mode("regression")
}

#' Create a workflow
#'
#' This function creates a workflow for training an XGBoost model.
#'
#' @param formula The formula for the model
#' @param model_spec The model specification
#' @return A workflow for training the model
#' @importFrom workflows workflow add_formula add_model
#' @export
create_workflow <- function(formula, model_spec) {
  workflows::workflow() %>%
    workflows::add_formula(formula) %>%
    workflows::add_model(model_spec)
}

#' Create cross-validation folds
#'
#' This function creates cross-validation folds for model training.
#'
#' @param data The dataset
#' @param folds Number of folds for cross-validation
#' @param repeats Number of repeats for cross-validation
#' @return Cross-validation folds
#' @importFrom rsample vfold_cv
#' @export
create_cv_folds <- function(data, folds = 5, repeats = 2) {
  rsample::vfold_cv(data, v = folds, repeats = repeats)
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
#' @param n_trees Number of trees for the model
#' @param tree_depth Maximum tree depth
#' @param min_n Minimum observations in terminal nodes
#' @param loss_reduction Minimum loss reduction to make further splits
#' @param sample_size Fraction of the training set to sample
#' @param mtry Number of variables to sample at each split
#' @param learn_rate Learning rate
#' @param seed_count Number of seed values to use
#' @param seed_for_data Seed value for data preparation
#' @return Trained XGBoost model and testing data
#' @importFrom parsnip boost_tree set_engine set_mode
#' @importFrom tune control_grid finalize_model select_best tune_grid
#' @importFrom rsample initial_split training testing vfold_cv
#' @importFrom workflows workflow add_formula add_model
#' @importFrom doParallel registerDoParallel
#' @import xgboost
#' @export
train_xgboost_model <- function(data, formula, grid, folds = 5, repeats = 2, n_trees = 2000,
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

  # Set seed for reproducibility
  set.seed(seed_for_data)

  # Split data into training and testing sets
  split <- rsample::initial_split(data)
  train_data <- rsample::training(split)
  test_data <- rsample::testing(split)

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
  resamples <- lapply(seeds, function(seed) {
    set.seed(seed)
    rsample::vfold_cv(train_data, v = folds, repeats = repeats)
  })

  # Tune grid
  doParallel::registerDoParallel()
  res <- tune::tune_grid(xgb_workflow, resamples = resamples, grid = grid, control = control_grid(save_pred = TRUE))

  # Train model with best parameters
  best_params <- tune::select_best(res, "rmse")
  trained_model <- tune::finalize_model(xgb_spec, best_params)

  return(list(trained_model = trained_model, test_data = test_data))
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


