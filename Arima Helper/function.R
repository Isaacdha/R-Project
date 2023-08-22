# Required Library
library(forecast)
library(tseries)
library(TSA)
library(readxl)
library(dplyr)

#--------------------------------------------------------------------------------------------------------------------------------------

# Transformation Function
stationarity_test_transform <- function(ts.data) {
  # Variance Stationarity
  lambda <- BoxCox.lambda(ts.data)
  lambda_list <- c(lambda)
  trans.data <- ts.data
  
  while (abs(tail(lambda_list, n = 1) - 1) > 0.1) {
    trans.data <- BoxCox(trans.data, tail(lambda_list, n = 1))
    lambda_list <- append(lambda_list, BoxCox.lambda(trans.data))
  }
  
  # Mean Stationarity
  n_diff <- 0
  diff.data <- trans.data 
  adf_pval <- c(adf.test(diff.data)$p.value)
  
  while (tail(adf_pval, n = 1) > 0.05){
    diff.data <- diff(diff.data)
    n_diff = n_diff + 1
    adf_pval <- append(adf_pval, adf.test(diff.data)$p.value)
  }
  
  listed <- list("trans.data" = trans.data, "lambda" = lambda_list, 
                 "mean.trans.data" = diff.data, "adf.pvalue" = adf_pval, 
                 "d" = n_diff)
  comment(listed) <- "Basic Variance and Mean Stationarity Transformation. Please use Var Stationarity Data in next step. Mean Stationarity only for search (d) order"
  return(listed)
}

#--------------------------------------------------------------------------------------------------------------------------------------

# Combination Generation Function
generate_combinations <- function(p, d, q, P = NULL, D = NULL, Q = NULL) {
  if (is.null(P) || is.null(D) || is.null(Q)) {
    combinations <- expand.grid(p, d, q)
    combinations <- combinations[!(combinations[, 1] == 0 
                                   & combinations[, 3] == 0), ]
    row.names(combinations) <- NULL
    colnames(combinations) <- c("Order p", "Order d", "Order q")
  } else {
    combinations <- expand.grid(p, d, q, P, D, Q)
    combinations <- combinations[!(combinations[, 1] == 0 
                                   & combinations[, 3] == 0), ]
    row.names(combinations) <- NULL
    colnames(combinations) <- c("Order p", "Order d", "Order q", "Order P", 
                                "Order D", "Order Q")
  }
  return(as.matrix(combinations))
}

#--------------------------------------------------------------------------------------------------------------------------------------
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

# Tentative Model List & Check
manual_arima <- function(data, p, d, q, s = NULL, 
                         P = NULL, D = NULL, Q = NULL, 
                         xreg = NULL) {
  
  # Parameter Combination
  parameter_combinations <- generate_combinations(p, d, q, P, D, Q)
  
  # Check if its intervention
  if (!is.null(xreg)){
    parameter_combinations <- rbind(parameter_combinations, 
                                    parameter_combinations)
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
    tryCatch({
      params <- parameter_combinations[i, ]
      
      if (is.null(s)){
        if (intv_status[i] == " " | is.null(xreg)) {
          model <- Arima(data, order = params[1:3])
        } else {
          model <- Arima(data, order = params[1:3], xreg = xreg)
        }
      } else
        if (intv_status[i] == " " | is.null(xreg)) {
          model <- Arima(data, order = params[1:3], 
                         seasonal = list(order=params[4:6],period=s))
        } else {
          model <- Arima(data, order = params[1:3], 
                         seasonal = list(order=params[4:6],period=s), 
                         xreg = xreg)
        }
      
      aicc_model[i] <- round(model$aicc, 3)
      bic_model[i] <- round(model$bic, 3)
      aic_model[i] <- round(model$aic, 3)
      
      coeftesting <- lmtest::coeftest(model)
      if (all(coeftesting[,4] < 0.05)) {
        coef_test[i] <- "PASS"
      } else {
        coef_test[i] <- "FAIL"
      }
      
      invisible(capture.output(residualtesting <- checkresiduals(model, 
                                                                 plot = FALSE)))
      if (residualtesting$p.value > 0.05) {
        ljung_test[i] <- paste0("PASS (p-value: ", 
                                sprintf("%.3f", residualtesting$p.value), ")")
      } else {
        ljung_test[i] <- paste0("FAIL (p-value: ", 
                                sprintf("%.3f", residualtesting$p.value), ")")
      }
      
      residualnormality <- shapiro.test(residuals(model))
      if (residualnormality$p.value > 0.05) {
        normal_test[i] <- paste("PASS (p-value:", 
                                sprintf("%.3f", residualnormality$p.value), ")")
      } else {
        normal_test[i] <- paste("FAIL (p-value:", 
                                sprintf("%.3f", residualnormality$p.value), ")")
      }
      
      status[i] <- "CLEAR"
    }, error = function(err) {
    })
  }
  
  # Combine Parameter into 1 Column
  parameter_combinations <- cbind(parameter_combinations, intv_status)
  
  if (is.null(s)){
    parameter_combinations <- apply(parameter_combinations, 1, function(x) {
      paste0("(", x[1], ", ", x[2], ", ", x[3], ") ", x[4])
    })
  } else {
    parameter_combinations <- apply(parameter_combinations, 1, function(x) {
      paste0("(", x[1], ", ", x[2], ", ", x[3], ")  (", x[4]
             , ", ", x[5], ", ", x[6], ")  ", s, "  ", x[7])
    })
  }
  
  # Make the Dataframe
  df <- data.frame(Model = parameter_combinations)
  df$"Status" = as.vector(status)
  df$"Coefficient Significance Test" = as.vector(coef_test)
  df$"Ljung-Box Residual Test" = as.vector(ljung_test)
  df$"Saphiro Wilk Residual Normality Test" = as.vector(normal_test)
  df$"AIC Test" = as.vector(aic_model)
  df$"AICC Test" = as.vector(aicc_model)
  df$"BIC Test" = as.vector(bic_model)
  
  arima_list <- list(df)
  names(arima_list) <- c("model_list")
  comment(arima_list) <- "Basic Coefficient Test, Ljung-Box Residual Test, and Saphiro-Wilk Residual Normality Test for each possible tentative model"
  
  return(arima_list)
}

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

# Transform Back the Transformed Data & Plot
transform_back <- function(forecast_object, lambda){
  n <- length(lambda)
  
  forecast_object$upper <- InvBoxCox(forecast_object$upper, lambda[1])
  forecast_object$lower <- InvBoxCox(forecast_object$lower, lambda[1])
  forecast_object$mean <- InvBoxCox(forecast_object$mean, lambda[1])
  forecast_object$x <- InvBoxCox(forecast_object$x, lambda[1])
  
  return(forecast_object)
}
