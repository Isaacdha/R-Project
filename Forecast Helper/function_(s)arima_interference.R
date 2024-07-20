################################################################################
# Required Library
################################################################################
library(forecast)
library(tseries)
library(TSA)
library(readxl)
library(dplyr)

################################################################################
# Transformation Function
################################################################################

# Transform & Scale
stationarity_test_transform <- function(ts.data, scale = 1) {
  # Apply scaling
  ts.data <- ts.data * scale
  
  # Variance Stationarity
  lambda_list <- numeric()
  trans.data <- ts.data
  
  repeat {
    lambda <- BoxCox.lambda(trans.data)
    trans.data <- BoxCox(trans.data, lambda)
    lambda_list <- c(lambda_list, lambda)
    
    if (abs(lambda - 1) <= 0.1) break
  }
  
  # Mean Stationarity
  n_diff <- 0
  diff.data <- trans.data 
  adf_pval <- adf.test(diff.data)$p.value
  
  while (tail(adf_pval, n = 1) > 0.05){
    diff.data <- diff(diff.data)
    n_diff = n_diff + 1
    adf_pval <- append(adf_pval, adf.test(diff.data)$p.value)
  }
  
  # Create a list to store lambdas and scale
  lambda_scale_list <- list(
    lambdas = lambda_list,
    scale = scale
  )
  
  listed <- list(
    trans.data = trans.data,
    lambda_scale_list = lambda_scale_list, 
    mean.trans.data = diff.data, 
    adf.pvalue = adf_pval, 
    d = n_diff
  )
  
  comment(listed) <- "Basic Variance and Mean Stationarity Transformation. 
  Please use Var Stationarity Data in next step. 
  Mean Stationarity only for search (d) order"
  
  return(listed)
}

#------------------------------------------------------------------------------#
# Transform Back the Transformed Data & Plot
transform_back <- function(forecast_object, lambda_scale_list) {
  
  # Extract the scale factor and lambda list
  scale <- lambda_scale_list$scale
  lambda_list <- lambda_scale_list$lambdas
  
  if (length(lambda_list) == 1) {
    lambda <- lambda_list[1]
    forecast_object$upper <- InvBoxCox(forecast_object$upper, lambda)
    forecast_object$lower <- InvBoxCox(forecast_object$lower, lambda)
    forecast_object$mean <- InvBoxCox(forecast_object$mean, lambda)
    forecast_object$x <- InvBoxCox(forecast_object$x, lambda)
  } else {
    for (i in seq_along(lambda_list)) {
      forecast_object$upper <- InvBoxCox(forecast_object$upper, lambda_list[i])
      forecast_object$lower <- InvBoxCox(forecast_object$lower, lambda_list[i])
      forecast_object$mean <- InvBoxCox(forecast_object$mean, lambda_list[i])
      forecast_object$x <- InvBoxCox(forecast_object$x, lambda_list[i])
    }
  }
  
  # Rescale the forecast object
  forecast_object$upper <- forecast_object$upper / scale
  forecast_object$lower <- forecast_object$lower / scale
  forecast_object$mean <- forecast_object$mean / scale
  forecast_object$x <- forecast_object$x / scale
  
  return(forecast_object)
}

################################################################################
# Manual ARIMA FUnction 
################################################################################

# Combination Generation Function
generate_combinations <- function(p, d, q, P = NULL, D = NULL, Q = NULL) {
  if (is.null(P) || is.null(D) || is.null(Q)) {
    combinations <- expand.grid(p, d, q)
    combinations <- combinations[!(combinations[, 1] == 0 & combinations[, 3] == 0), ]
    row.names(combinations) <- NULL
    colnames(combinations) <- c("Order p", "Order d", "Order q")
  } else {
    combinations <- expand.grid(p, d, q, P, D, Q)
    combinations <- combinations[!(combinations[, 1] == 0 & combinations[, 3] == 0), ]
    row.names(combinations) <- NULL
    colnames(combinations) <- c("Order p", "Order d", "Order q", "Order P", "Order D", "Order Q")
  }
  return(as.matrix(combinations))
}

#------------------------------------------------------------------------------#
# Intervention Helper
generate_intervention <- function(ts_data, intervention_point, 
                                  types = c("step", "ramp", "pulse"), 
                                  pulse_start = NULL, pulse_end = NULL, 
                                  mode = c("regular", "forecast"), 
                                  previous_data_last_number = NULL) {
  # Determine the mode
  mode <- match.arg(mode)
  
  # Check if ts_data is a time series or a number (length of the array)
  if (is.numeric(ts_data) && length(ts_data) == 1) {
    time_index <- seq_len(ts_data)
  } else {
    time_index <- time(ts_data)
  }
  
  # Create intervention variables based on the specified types
  interventions <- list()
  
  if (mode == "regular") {
    if ("step" %in% types) {
      step_intervention <- ifelse(time_index >= intervention_point, 1, 0)
      interventions$step_intervention <- step_intervention
    }
    
    if ("ramp" %in% types) {
      ramp_intervention <- pmax(0, time_index - intervention_point + 1)
      interventions$ramp_intervention <- ramp_intervention
    }
    
    if ("pulse" %in% types) {
      if (is.null(pulse_start) | is.null(pulse_end)) {
        stop("For pulse intervention, both pulse_start and pulse_end need to be specified.")
      }
      pulse_intervention <- ifelse(time_index >= pulse_start & time_index <= pulse_end, 1, 0)
      interventions$pulse_intervention <- pulse_intervention
    }
  } else if (mode == "forecast") {
    if ("step" %in% types) {
      step_intervention <- rep(1, length(time_index))
      interventions$step_intervention <- step_intervention
    }
    
    if ("ramp" %in% types) {
      if (is.null(previous_data_last_number)) {
        stop("For ramp intervention in forecast mode, previous_data_last_number must be specified.")
      }
      last_ramp_value <- previous_data_last_number
      ramp_intervention <- seq(last_ramp_value + 1, last_ramp_value + length(time_index))
      interventions$ramp_intervention <- ramp_intervention
    }
    
    if ("pulse" %in% types) {
      pulse_intervention <- rep(0, length(time_index))
      interventions$pulse_intervention <- pulse_intervention
    }
  }
  
  # Combine interventions into a data frame
  intervention_df <- as.data.frame(interventions)
  
  return(intervention_df)
}

#------------------------------------------------------------------------------#
# Tentative Model List & Check
manual_arima <- function(data, p, d, q, s = NULL, 
                         P = NULL, D = NULL, Q = NULL, xreg = NULL) {
  # Parameter Combination
  parameter_combinations <- generate_combinations(p, d, q, P, D, Q)
  
  # Check if it's intervention
  if (!is.null(xreg)) {
    parameter_combinations <- rbind(parameter_combinations, parameter_combinations)
    half_length <- nrow(parameter_combinations) / 2
    intv_status <- matrix(c(rep("INTV", half_length), rep(" ", half_length)), ncol = 1)
  } else {
    intv_status <- matrix(" ", nrow = nrow(parameter_combinations))
  }
  
  # Make Empty Array
  status <- matrix("ERROR", nrow = nrow(parameter_combinations))
  coef_test <- matrix("-", nrow = nrow(parameter_combinations))
  ljung_test <- matrix("-", nrow = nrow(parameter_combinations))
  normal_test <- matrix("-", nrow = nrow(parameter_combinations))
  aic_model <- matrix("-", nrow = nrow(parameter_combinations))
  aicc_model <- matrix("-", nrow = nrow(parameter_combinations))
  bic_model <- matrix("-", nrow = nrow(parameter_combinations))
  
  for (i in 1:nrow(parameter_combinations)) {
    params <- parameter_combinations[i, ]
    cat("Processing combination:", params, "\n")
    
    tryCatch({
      withCallingHandlers({
        if (is.null(s)) {
          if (intv_status[i] == " " | is.null(xreg)) {
            model <- Arima(data, order = params[1:3])
          } else {
            model <- Arima(data, order = params[1:3], xreg = xreg)
          }
        } else {
          if (intv_status[i] == " " | is.null(xreg)) {
            model <- Arima(data, order = params[1:3], seasonal = list(order=params[4:6], period=s))
          } else {
            model <- Arima(data, order = params[1:3], seasonal = list(order=params[4:6], period=s), xreg = xreg)
          }
        }
        
        aicc_model[i] <- round(model$aicc, 3)
        bic_model[i] <- round(model$bic, 3)
        aic_model[i] <- round(model$aic, 3)
        
        coeftesting <- lmtest::coeftest(model)
        if (all(coeftesting[, 4] < 0.05)) {
          coef_test[i] <- "PASS"
        } else {
          coef_test[i] <- "FAIL"
        }
        
        residualtesting <- checkresiduals(model, plot = FALSE)
        if (residualtesting$p.value > 0.05) {
          ljung_test[i] <- paste0("PASS (p-value: ", sprintf("%.3f", residualtesting$p.value), ")")
        } else {
          ljung_test[i] <- paste0("FAIL (p-value: ", sprintf("%.3f", residualtesting$p.value), ")")
        }
        
        residualnormality <- shapiro.test(residuals(model))
        if (residualnormality$p.value > 0.05) {
          normal_test[i] <- paste("PASS (p-value:", sprintf("%.3f", residualnormality$p.value), ")")
        } else {
          normal_test[i] <- paste("FAIL (p-value:", sprintf("%.3f", residualnormality$p.value), ")")
        }
        
        status[i] <- "CLEAR"
      }, warning = function(war) {
        cat("Warning for combination:", params, "\n")
        cat("Warning message:", war$message, "\n")
        invokeRestart("muffleWarning")
      })
    }, error = function(err) {
      cat("Error for combination:", params, "\n")
      cat("Error message:", err$message, "\n")
    })
  }
  
  # Combine Parameter into 1 Column
  parameter_combinations <- cbind(parameter_combinations, intv_status)
  
  if (is.null(s)) {
    parameter_combinations <- apply(parameter_combinations, 1, function(x) {
      paste0("(", x[1], ", ", x[2], ", ", x[3], ") ", x[4])
    })
  } else {
    parameter_combinations <- apply(parameter_combinations, 1, function(x) {
      paste0("(", x[1], ", ", x[2], ", ", x[3], ")  (", x[4], ", ", x[5], ", ", x[6], ")  ", s, "  ", x[7])
    })
  }
  
  # Make the Dataframe
  df <- data.frame(Model = parameter_combinations)
  df$"Status" <- as.vector(status)
  df$"Coefficient Significance Test" <- as.vector(coef_test)
  df$"Ljung-Box Residual Test" <- as.vector(ljung_test)
  df$"Shapiro-Wilk Residual Normality Test" <- as.vector(normal_test)
  df$"AIC Test" <- as.vector(aic_model)
  df$"AICC Test" <- as.vector(aicc_model)
  df$"BIC Test" <- as.vector(bic_model)
  
  arima_list <- list(df)
  names(arima_list) <- c("model_list")
  comment(arima_list) <- "Basic Coefficient Test, Ljung-Box Residual Test, and Shapiro-Wilk Residual Normality Test for each possible tentative model"
  
  return(arima_list)
}

#------------------------------------------------------------------------------#
# Model Select Based on Assumption Test & AIC BIC AICC
model_selector <- function(arima_list){
  df <- arima_list$model_list
  df <- df[df["Status"] == "CLEAR", ]
  df <- df[df["Coefficient Significance Test"] == "PASS", ]
  df %>% dplyr::filter(substr("Ljung-Box Residual Test",1,4) == "PASS")
  df <- df[order(df[,6], df[,7],df[,8]), ]
  best <- df[1,]
  
  if (is.null(best)){
    df <- arima_list$model_list
    df <- df[df["Status"] == "CLEAR", ]
    df <- df[df["Coefficient Significance Test"] == "PASS", ]
    df <- df[order(df[,6], df[,7],df[,8]), ]
    best <- df[1,]
  }
  
  if (is.null(best)){
    df <- arima_list$model_list
    df <- df[df["Status"] == "CLEAR", ]
    df <- df[order(df[,6], df[,7],df[,8]), ]
    best <- df[1,]
  }
  
  return(best)
}

################################################################################
# Unfinished Features
################################################################################
# Lag Sequence Detector
containsMultiples <- function(arr, num) {
  target_numbers <- num * c(1, 2, 3)
  
  # Check if all target numbers are present in the array
  result <- all(purrr::map_lgl(target_numbers, ~ .x %in% arr))
  return(result)
}

# Seasonal Lag Detector
divide_and_filter <- function(arr, num) {
  divided <- arr / num
  rounded <- round(divided)
  is_integer <- rounded == divided
  integer_divisions <- arr[is_integer]
  integer_divisions <- integer_divisions/num
  return(integer_divisions)
}

# Add Zero in Array Function
addZero <- function(inputArray) {
  if (length(inputArray) == 0) {
    return(0)
  } else {
    return(c(0, inputArray))
  }
}

# Lag Detector, Auto Transform if Non-Stationary
order_detector <- function(data, lag, s = NULL){
  data.trans <- data
  acf.positions <- c()
  pacf.positions <- c()
  sacf.positions <- c()
  spacf.positions <- c()
  
  ## Non-Seasonal
  acf_val <- acf(data, lag.max = lag, plot = FALSE)$acf
  cv <- 1.96/sqrt(length(data))
  acf.positions <- which(abs(acf_val) > cv)
  d = 0
  
  while (containsMultiples(acf.positions, 1)){
    data.trans <- diff(data.trans)
    d = d + 1
    acf_val <- acf(data.trans, lag.max = lag, plot = FALSE)$acf
    acf.positions <- which(abs(acf_val) > cv)
  }
  
  pacf_val <- pacf(data.trans, lag.max = lag, plot = FALSE)$acf
  pacf.positions <- which(abs(pacf_val) > cv)
  
  ## Seasonal
  data.trans_s <- data.trans
  
  if (!is.null(s)){
    D = 0
    sacf.positions <- acf.positions
    while (containsMultiples(sacf.positions, s)){
      data.trans_s <- diff(data.trans_s, s)
      D = D + 1
      sacf_val <- acf(data.trans_s, lag.max = lag, plot = FALSE)$acf
      sacf.positions <- which(abs(sacf_val) > cv)
    }
    
    spacf_val <- pacf(data.trans, lag.max = lag, plot = FALSE)$acf
    spacf.positions <- which(abs(spacf_val) > cv)
  }
  
  # Iterate over each array
  sacf.positions <- divide_and_filter(sacf.positions, s)
  spacf.positions <- divide_and_filter(spacf.positions, s)

  acf.positions <- addZero(acf.positions)
  pacf.positions <- addZero(pacf.positions)
  sacf.positions <- addZero(sacf.positions)
  spacf.positions <- addZero(spacf.positions)

  listed <- list("p" = pacf.positions, "d" = d, "q" = acf.positions, 
                 "P" = spacf.positions, "D" = D, "Q" = sacf.positions)
  return(listed)
}

################################################################################



