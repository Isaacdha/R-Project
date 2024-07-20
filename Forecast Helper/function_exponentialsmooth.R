################################################################################
# Required Library & Initial Settings
################################################################################
library(readxl)
library(ggplot2)
library(forecast)
library(ggeasy)
library(dplyr)
library(openxlsx)
library(zoo)

theme_update(plot.title = element_text(hjust = 0.5))
theme_update(text = element_text(size = 13))
options(warn = -1)

################################################################################
# Plot Function
################################################################################
# Function to plot fitted vs. actual and forecast
plot_forecast <- function(model, forecastnum, fitted_title, forecast_title, actual_title, ts_data, 
                          x_axis_name, y_axis_name, include_fitted_in_forecast, 
                          xrange, yrange, CI, legendtitle, actuallegend,
                          fittedlegend, forecastlegend) {
  fitted_data <- fitted(model)
  forecast_data <- forecast(model, h = forecastnum)
  
  # Plot Actual Data
  p0 <- autoplot(ts_data, series = actuallegend) +
    labs(title = actual_title, x = x_axis_name, y = y_axis_name, colour = legendtitle) +
    theme_minimal() + ggeasy::easy_center_title() +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  
  if (!is.null(xrange)) {
    p0 <- p0 + xlim(xrange)
  }
  if (!is.null(yrange)) {
    p0 <- p0 + ylim(yrange)
  }
  
  # Plot Fitted vs Actual
  p1 <- autoplot(ts_data, series = actuallegend) +
    autolayer(fitted_data, series = fittedlegend, PI = FALSE) +
    labs(title = fitted_title, x = x_axis_name, y = y_axis_name, colour = legendtitle) +
    theme_minimal() + ggeasy::easy_center_title() +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  
  if (!is.null(xrange)) {
    p1 <- p1 + xlim(xrange)
  }
  if (!is.null(yrange)) {
    p1 <- p1 + ylim(yrange)
  }
  
  # Plot Forecast
  p2 <- autoplot(ts_data, series = actuallegend) +
    autolayer(forecast_data, series = forecastlegend, PI = CI) +
    labs(title = forecast_title, x = x_axis_name, y = y_axis_name, colour = legendtitle) +
    theme_minimal() + ggeasy::easy_center_title() +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  
  if (include_fitted_in_forecast) {
    p2 <- p2 + autolayer(fitted_data, series = fittedlegend, PI = FALSE)
  }
  
  if (!is.null(xrange)) {
    p2 <- p2 + coord_cartesian(xlim = xrange)
  }
  if (!is.null(yrange)) {
    p2 <- p2 + coord_cartesian(ylim = yrange)
  }
  
  return(list(actual_plot = p0, fitted_plot = p1, forecast_plot = p2))
}

################################################################################
# Combined Exponential Smoothing Function
################################################################################
# Exponential Smoothing Function
exponential_smoothing <- function(ts_data, method, h = 60, 
                                  alpha = NULL, beta = NULL, gamma = NULL, 
                                  seasonal = "additive", 
                                  fitted_title = NULL, forecast_title = NULL, 
                                  actual_title = "Actual Data Plot",
                                  x_axis_name = NULL, y_axis_name = NULL, 
                                  include_fitted_in_forecast = FALSE, 
                                  xrange = NULL, yrange = NULL, CI = TRUE, 
                                  legendtitle = "Legend", 
                                  actuallegend = "Actual",
                                  fittedlegend = "Fitted", 
                                  forecastlegend = "Forecast") {
  if (method == "single") {
    model <- ses(ts_data, h = h, alpha = alpha)
    default_title <- "Single Exponential Smoothing"
  } else if (method == "double") {
    model <- holt(ts_data, h = h, alpha = alpha, beta = beta)
    default_title <- "Holt's Linear Trend Model"
  } else if (method == "triple") {
    model <- hw(ts_data, seasonal = seasonal, h = h, alpha = alpha, beta = beta, gamma = gamma)
    default_title <- ifelse(seasonal == "additive", "Holt-Winters Method (Additive)", "Holt-Winters Method (Multiplicative)")
  } else {
    stop("Invalid method. Choose 'single', 'double', or 'triple'.")
  }
  
  x_axis_name <- ifelse(is.null(x_axis_name), "Time", x_axis_name)
  y_axis_name <- ifelse(is.null(y_axis_name), deparse(substitute(ts_data)), y_axis_name)
  
  fitted_title <- ifelse(is.null(fitted_title), paste(default_title, "- Fitted vs Actual"), fitted_title)
  forecast_title <- ifelse(is.null(forecast_title), paste(default_title, "- Forecast"), forecast_title)
  
  plots <- plot_forecast(model, h, fitted_title, forecast_title, actual_title, ts_data, 
                         x_axis_name, y_axis_name, include_fitted_in_forecast, 
                         xrange, yrange, CI, legendtitle, actuallegend,
                         fittedlegend, forecastlegend)
  
  fitted_values <- fitted(model)
  forecast_values <- forecast(model, h = h)
  
  # Calculate Error Metrics
  errors <- fitted_values - ts_data
  RMSE <- sqrt(mean(errors^2, na.rm = TRUE))
  MSE <- mean(errors^2, na.rm = TRUE)
  MAE <- mean(abs(errors), na.rm = TRUE)
  MAPE <- mean(abs(errors / ts_data) * 100, na.rm = TRUE)
  
  result <- list(
    alpha = model$model$par["alpha"],
    beta = model$model$par["beta"],
    gamma = model$model$par["gamma"],
    fitted_values = fitted_values,
    real_values = ts_data,
    forecast_values = forecast_values$mean,
    RMSE = RMSE,
    MSE = MSE,
    MAE = MAE,
    MAPE = MAPE,
    actual_plot = plots$actual_plot,
    fitted_plot = plots$fitted_plot,
    forecast_plot = plots$forecast_plot
  )
  
  return(result)
}

################################################################################
# Export Function
################################################################################
# Export to Excel
export_to_excel <- function(result, path = "D:/Result.xlsx") {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Sheet 1: Parameters and Errors
  addWorksheet(wb, "Parameters and Errors")
  param_error_data <- data.frame(
    Parameter = c("Alpha", "Beta", "Gamma"),
    Value = c(result$alpha, result$beta, result$gamma)
  )
  error_data <- data.frame(
    Error_Metric = c("RMSE", "MSE", "MAE", "MAPE"),
    Value = c(result$RMSE, result$MSE, result$MAE, result$MAPE)
  )
  writeData(wb, sheet = "Parameters and Errors", x = param_error_data, startCol = 1, startRow = 1, colNames = TRUE)
  writeData(wb, sheet = "Parameters and Errors", x = error_data, startCol = 4, startRow = 1, colNames = TRUE)
  
  # Extract time information and convert to Date format
  real_time_index <- as.Date(time(result$real_values))
  fitted_time_index <- as.Date(time(result$fitted_values))
  forecast_time_index <- as.Date(time(result$forecast_values))
  
  # Sheet 2: Fitted vs Real Values
  addWorksheet(wb, "Fitted vs Real Values")
  fitted_vs_real <- data.frame(
    Index = seq_along(result$real_values),
    Time = real_time_index,
    Real_Value = as.numeric(result$real_values),
    Fitted_Value = as.numeric(result$fitted_values)
  )
  writeData(wb, sheet = "Fitted vs Real Values", x = fitted_vs_real, startCol = 1, startRow = 1, colNames = TRUE)
  
  # Save the fitted vs real values plot
  fitted_plot_path <- "fitted_plot.png"
  ggsave(fitted_plot_path, plot = result$fitted_plot, width = 8, height = 6)
  insertImage(wb, "Fitted vs Real Values", fitted_plot_path, startRow = nrow(fitted_vs_real) + 3, startCol = 1)
  
  # Sheet 3: Forecasted Values
  addWorksheet(wb, "Forecasted Values")
  forecast_values <- data.frame(
    Index = seq_along(result$forecast_values),
    Time = forecast_time_index,
    Forecasted_Value = as.numeric(result$forecast_values)
  )
  writeData(wb, sheet = "Forecasted Values", x = forecast_values, startCol = 1, startRow = 1, colNames = TRUE)
  
  # Save the forecast plot
  forecast_plot_path <- "forecast_plot.png"
  ggsave(forecast_plot_path, plot = result$forecast_plot, width = 8, height = 6)
  insertImage(wb, "Forecasted Values", forecast_plot_path, startRow = nrow(forecast_values) + 3, startCol = 1)
  
  # Sheet 4: Actual Data Plot
  addWorksheet(wb, "Actual Data Plot")
  actual_plot_path <- "actual_plot.png"
  ggsave(actual_plot_path, plot = result$actual_plot, width = 8, height = 6)
  insertImage(wb, "Actual Data Plot", actual_plot_path, startRow = 1, startCol = 1)
  
  # Save the workbook
  saveWorkbook(wb, file = path, overwrite = TRUE)
}

################################################################################
