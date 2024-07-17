library(readxl)
library(ggplot2)
library(forecast)
library(dplyr)
options(warn = -1)

# Function to plot fitted vs. actual and forecast
plot_forecast <- function(model, fitted_title, forecast_title, ts_data, 
                          x_axis_name, y_axis_name, include_fitted_in_forecast, 
                          xrange, yrange, CI) {
  fitted_data <- fitted(model)
  forecast_data <- forecast(model, h = 60)
  
  # Plot Fitted vs Actual
  p1 <- autoplot(ts_data, series = "Actual") +
    autolayer(fitted_data, series = "Fitted", PI = FALSE) +
    labs(title = fitted_title, x = x_axis_name, y = y_axis_name, colour = "Legend") +
    theme_minimal()
  
  # Apply x and y axis limits if specified
  if (!is.null(xrange)) {
    p1 <- p1 + xlim(xrange)
  }
  if (!is.null(yrange)) {
    p1 <- p1 + ylim(yrange)
  }
  
  # Plot Forecast
  p2 <- autoplot(ts_data, series = "Actual") +
    autolayer(forecast_data, series = "Forecast", PI = CI)
  
  if (include_fitted_in_forecast) {
    p2 <- p2 + autolayer(fitted_data, series = "Fitted", PI = FALSE)
  }
  
  # Apply x and y axis limits if specified
  if (!is.null(xrange)) {
    p2 <- p2 + xlim(xrange)
  }
  if (!is.null(yrange)) {
    p2 <- p2 + ylim(yrange)
  }
  
  p2 <- p2 + labs(title = forecast_title, x = x_axis_name, y = y_axis_name, colour = "Legend") +
    theme_minimal()
  
  return(list(fitted_plot = p1, forecast_plot = p2))
}

# Combined Exponential Smoothing Function
exponential_smoothing <- function(ts_data, method, h = 60, 
                                  alpha = NULL, beta = NULL, gamma = NULL, 
                                  seasonal = "additive", 
                                  fitted_title = NULL, forecast_title = NULL, 
                                  x_axis_name = NULL, y_axis_name = NULL, 
                                  include_fitted_in_forecast = FALSE, 
                                  xrange = NULL, yrange = NULL, CI = TRUE) {
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
  
  plots <- plot_forecast(model, fitted_title, forecast_title, ts_data, x_axis_name, y_axis_name, include_fitted_in_forecast, xrange, yrange, CI)
  
  fitted_values <- fitted(model)
  forecast_values <- forecast(model, h = h)
  
  result <- list(
    alpha = model$model$par["alpha"],
    beta = model$model$par["beta"],
    gamma = model$model$par["gamma"],
    fitted_values = fitted_values,
    real_values = ts_data,
    forecast_values = forecast_values$mean,
    fitted_plot = plots$fitted_plot,
    forecast_plot = plots$forecast_plot
  )
  
  return(result)
}