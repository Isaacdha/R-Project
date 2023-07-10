library(forecast)
library(tseries)
library(TSA)
library(readxl)

generate_combinations <- function(p, d, q, P = NULL, D = NULL, Q = NULL) {
  # Filter combinations where p and q are both 0
  if (is.null(P) || is.null(D) || is.null(Q)) {
    combinations <- expand.grid(p, d, q)
    combinations <- combinations[!(combinations[, 1] == 0 & combinations[, 3] == 0), ]
  } else {
    combinations <- expand.grid(p, d, q, P, D, Q)
    combinations <- combinations[!(combinations[, 1] == 0 & combinations[, 3] == 0), ]
  }
  
  row.names(combinations) <- NULL
  colnames(combinations) <- c("Order p", "Order d", "Order q", "Order P", "Order D", "Order Q")
  return(as.matrix(combinations))
}

manual_arima <- function(data, p, d, q, s, P = NULL, D = NULL, Q = NULL) {
  parameter_combinations <- generate_combinations(p, d, q, P, D, Q)
  status <- matrix(NA, nrow = nrow(parameter_combinations))
  coef_test <- matrix(NA, nrow = nrow(parameter_combinations))
  ljung_test <- matrix(NA, nrow = nrow(parameter_combinations))
  normal_test <- matrix(NA, nrow = nrow(parameter_combinations))
  aicc_model <- matrix(NA, nrow = nrow(parameter_combinations))
  bic_model <- matrix(NA, nrow = nrow(parameter_combinations))
  
  for (i in 1:nrow(parameter_combinations)) {
    tryCatch({
      params <- parameter_combinations[i, ]
      model <- Arima(data, order = params[1:3], seasonal = list(order=params[4:6],period=s))
      aicc_model[i] <- model$aicc
      bic_model[i] <- model$bic
      
      coeftesting <- lmtest::coeftest(model)
      if (all(coeftesting[,4] < 0.05)) {
        coef_test[i] <- "PASS"
      } else {
        coef_test[i] <- "FAIL"
      }
      
      residualtesting <- checkresiduals(model, plot = FALSE)
      if (residualtesting$p.value > 0.05) {
        ljung_test[i] <- "PASS"
      } else {
        ljung_test[i] <- "FAIL"
      }
      
      residualnormality <- shapiro.test(residuals(model))
      if (residualnormality$p.value > 0.05) {
        normal_test[i] <- "PASS"
      } else {
        normal_test[i] <- "FAIL"
      }
      status[i] <- "CLEAR"
    }, error = function(e) {
      status[i] <- "ERROR"
    })
  }
  
  # Combine the parameters into a single column
  parameter_combinations <- apply(parameter_combinations, 1, function(x) {
    paste0("(", x[1], ", ", x[2], ", ", x[3], ") (", x[4], ", ", x[5], ", ", x[6], ") ")
  })
  
  df <- data.frame(Parameter = parameter_combinations)
  df$"Status" = as.vector(status)
  df$"Coefficient Test" = as.vector(coef_test)
  df$"Ljung-Box Test" = as.vector(ljung_test)
  df$"Residual Normality Test" = as.vector(normal_test)
  df$"AICC Test" = as.vector(aicc_model)
  df$"BIC Test" = as.vector(bic_model)
 
  
  return(df)
}
