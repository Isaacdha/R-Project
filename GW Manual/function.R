## Install Packages
install.packages("xlsx")
install.packages("readxl")
install.packages("sp")
install.packages("spdep")
install.packages("ggplot2")
install.packages("MASS")
install.packages("dplyr")
install.packages("magrittr")
  
## Load Package
library(xlsx)
library(readxl)
library(sp)
library(spdep)
library(ggplot2)
library(MASS)
library(dplyr)
library(magrittr)

## Buat Fungsi Otomatis
pemetaan <- function(map, data, tahun_peta){
  map$KABKOTNO <- as.numeric(map$KABKOTNO)
  data_map <- data[,c("no_kabkot","tahun","produksi")]
  colnames(data_map) <- c("KABKOTNO", "tahun", "produksi")
  
  if (is.numeric(tahun_peta)){
    data_map <- data_map[data_map$tahun == tahun_peta, ]
    judul <- paste("Heatmap Produksi Jawa Timur tahun", tahun_peta)
  }
  else if (tahun_peta == "agregat"){
    data_map <- aggregate(produksi ~ KABKOTNO, data = data_map, FUN = mean)
    judul <- "Heatmap Produksi Jawa Timur agregat"
  }
  
  data_map$KABKOTNO[30:length(data_map$KABKOTNO)] <- 
    seq(71, length.out = length(data_map$KABKOTNO)-29)
  data_map <- merge(map, data_map, by = "KABKOTNO")
  
  heatmap <- ggplot() +
    geom_sf(data = data_map, aes(fill = produksi)) +
    scale_fill_gradient(low = "white", high = "green") +
    labs(title = judul, fill = "Produksi")
  
  return(heatmap)
}

################################################################################
## Fungsi Permodelan GWPR MEST
################################################################################

## Fungsi Ambil Variabel & Sortir Data
preprocess_gw <- function(data, idx_area, idx_waktu=NULL, VarY, VarX, Longitude, 
                          Latitude) {
  # Sort Data disesuaikan dengan jenis data masuk
  if (is.null(idx_waktu)) {
    sorted_data <- data[order(data[[idx_area]]), ]
  } else {
    sorted_data <- data[order(data[[idx_area]], data[[idx_waktu]]), ]
  }
  
  sorted_data[[idx_area]] <- as.integer(sorted_data[[idx_area]])
  
  if ("no" %in% names(sorted_data)) {
    sorted_data[["no"]] <- as.integer(sorted_data[["no"]])
  }
  
  
  if (!is.null(idx_waktu)) {
    sorted_data[[idx_waktu]] <- as.integer(sorted_data[[idx_waktu]])
    n_waktu <- sorted_data[[idx_waktu]] %>% n_distinct()
  } else {
    n_waktu <- 1
  }
  
  ## Jumlah Kota & Variabel
  n_area <- sorted_data[[idx_area]] %>% n_distinct()
  n <- nrow(sorted_data)
  
  ## Variabel X, Y, dan banyaknya Prediktor (p)
  y <- sorted_data[[VarY]] %>% as.matrix() 
  x <- sorted_data[, VarX] %>% cbind(1, .) %>% as.matrix()   
  p <- ncol(x)
  
  ## Atur Koordinat
  U <- sorted_data[, Longitude] %>% as.matrix() # Baca & Convert Longitude
  V <- sorted_data[, Latitude] %>% as.matrix()  # Baca & Convert Latitude
  
  ## Menghitung Jarak Euclidean
  d <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      d[i, j] <- sqrt(((U[i] - U[j])^2) + 
                        ((V[i] - V[j])^2))
    }
  }
  
  result <- list(data = sorted_data, n_area = n_area, n_waktu = n_waktu, n = n, 
                 y = y, x = x, p = p, dist_mat = d)
  return(result)
}

## Fungsi M-Estimator Tukey Biweight
biweight_tukey_mest <- function(x, y, delta_inp, max_iter = 1000, eps = 1e-3) {
  n <- nrow(x)
  p <- ncol(x)
  w <- rep(1, n)
  delta <- delta_inp
  for (iter in 1:max_iter) {
    delta_prev <- delta
    residuals <- y - x %*% delta
    scaled_residuals <- ifelse(residuals == 0, 0, 
                               residuals / (4.685 * mad(residuals)))
    weights <- (1 - (scaled_residuals / 4.685)^2)^2
    weights[abs(scaled_residuals) >= 4.685] <- 0
    
    W_tukey <- matrix(0, nrow = n, ncol = n)
    diag(W_tukey) <- weights
    
    tryCatch({
      delta <- solve(t(x) %*% W_tukey %*% x) %*% t(x) %*% W_tukey %*% y
    }, error = function(e) {
      delta <- ginv(t(x) %*% W_tukey %*% x) %*% t(x) %*% W_tukey %*% y
    })
    
    if (max(abs(delta - delta_prev)) < eps) {
      break
    }
  }
  return(delta)
}

## Fungsi Bandwith Adaptive M-Estimator
cari_bw_adapt <- function(data, dist_mat, n_area, n_data, n_waktu=1, 
                          idx_area, x, y, length_params = 25){
  # Inisialisasi Vektor Matrix
  CVmin <- c() # Vektor untuk menyimpan nilai minimum CV
  l_ab <- c() # Vektor untuk menyimpan nilai bandwith setiap area
  
  # Iterasi untuk setiap area
  for (i in 1:n_area){
    # Inisialisasi variabel-variabel awal
    A <- 0.0001
    B <- max(max(dist_mat))
    iter_ab <- 0
    minCV <- 0
    selisih <- 1000
    
    # Melakukan iterasi sampai selisih kecil (konvergen) atau terlalu lama
    while (selisih > 0.0001 && iter_ab <= 1000) {
      # Inisialisasi awalan bandwidth (l_awal) dengan length_params titik pada rentang A B
      l_awal <- seq(A, B, length.out = length_params)
      nl <- length(l_awal) # Banyaknya awalan bandwidth
      CV <- numeric(nl)    # Vektor untuk menyimpan nilai CV
      
      # Perhitungan CV untuk setiap l_awal
      for (k in 1:nl) {
        tryCatch({
          l <- l_awal[k]
          sigma <- l # Setel bandwidth sebagai (sigma) untuk kernel Gaussian
          Wb <- exp(-(dist_mat^2) / (2 * sigma^2)) # Kernel Gaussian
          
          # Temukan observasi yang terkait dengan area i dan ditiadakan
          kota <- which(data[, idx_area] == i)
          
          W <- diag(Wb[i, ])
          W <- as.matrix(W)
          W <- W[-kota, ]
          W <- W[, -kota]
          
          x_cv <- x
          x_cv <- as.matrix(x_cv)
          x_cv <- x_cv[-kota, ]
          
          y_cv <- y
          y_cv <- as.matrix(y_cv)
          y_cv <- y_cv[-kota]
          
          # Estimasi koefisien beta_cv menggunakan MKT
          beta_cv <- ginv(t(x_cv) %*% W %*% x_cv) %*% t(x_cv) %*% W %*% y_cv
          
          # Lanjutkan menggunakan iterasi tukey biweight jika m-estimator
          # Update beta_cv dengan menggunakan bobot W_tukey
          beta_cv <- biweight_tukey_mest(x_cv, y_cv, beta_cv)
          
          # Hitung yhat_cv
          yhat_cv <- x[kota, ] %*% beta_cv
          
          # Hitung CV dan simpan nilai di CV(k)
          CV[k] <- sum((y[kota] - yhat_cv)^2)
        }, error = function(e) {
          return()
        })
      }
      
      hasilCV <- cbind(l_awal, CV)
      A0 <- A
      B0 <- B
      minCV <- min(CV, na.rm = TRUE)    
      l_min <- which(CV == minCV)[1]
      
      if (l_min == 1) {
        A <- l_awal[l_min]
        B <- l_awal[l_min + 1]
      } else if (l_min == nl) {
        A <- l_awal[l_min - 1]
        B <- l_awal[l_min]
      } else {
        A <- l_awal[l_min - 1]
        B <- l_awal[l_min + 1]
      }
      
      selisih <- (B0 - A0) - (B - A)  # Selisih antara A0-B0 dan A-B
      iter_ab <- iter_ab + 1
      
      cat("Iterasi BW (Index Area, Iterasi, Selisih):\n")
      cat(i, iter_ab, selisih, "\n")  # Tampilkan informasi iterasi
    }
    
    hasilCV <- hasilCV[order(hasilCV[, 2]), ]  
    l_ab[i] <- hasilCV[1, 1]                  
    CVmin[i] <- hasilCV[1, 2]   
  }
  return(l_ab)
}

## Fungsi Bandwith Fixed M-Estimator
cari_bw_fixed <- function(data, dist_mat, n_area, n_data, n_waktu=1, 
                          idx_area, x, y, length_params = 25){
  # Inisialisasi Vektor Matrix
  CVmin <- c() # Vektor untuk menyimpan nilai minimum CV
  l_ab <- c() # Vektor untuk menyimpan nilai bandwith setiap area
  
  # Inisialisasi variabel-variabel awal
  A <- 0.0001
  B <- max(max(dist_mat))
  iter_ab <- 0
  minCV <- 0
  selisih <- 1000
    
  # Melakukan iterasi sampai selisih kecil (konvergen) atau terlalu lama
  while (selisih > 0.0001 && iter_ab <= 1000) {
    # Inisialisasi awalan bandwidth (l_awal) dengan length_params titik pada rentang A B
    l_awal <- seq(A, B, length.out = length_params)
    nl <- length(l_awal) # Banyaknya awalan bandwidth
    CV <- numeric(nl)    # Vektor untuk menyimpan nilai CV
      
    # Perhitungan CV untuk setiap l_awal
    for (k in 1:nl) {
      tryCatch({
        l <- l_awal[k]
        sigma <- l # Setel bandwidth sebagai (sigma) untuk kernel Gaussian
        Wb <- exp(-(dist_mat^2) / (2 * sigma^2)) # Kernel Gaussian
        
        for (i in 1:n_area){
          # Temukan observasi yang terkait dengan area i dan ditiadakan
          kota <- which(data[, idx_area] == i)
          
          W <- diag(Wb[i, ])
          W <- as.matrix(W)
          W <- W[-kota, ]
          W <- W[, -kota]
          
          x_cv <- x
          x_cv <- as.matrix(x_cv)
          x_cv <- x_cv[-kota, ]
          
          y_cv <- y
          y_cv <- as.matrix(y_cv)
          y_cv <- y_cv[-kota]
          
          # Estimasi koefisien beta_cv menggunakan MKT
          beta_cv <- ginv(t(x_cv) %*% W %*% x_cv) %*% t(x_cv) %*% W %*% y_cv
          
          # Lanjutkan menggunakan iterasi tukey biweight jika m-estimator
          # Update beta_cv dengan menggunakan bobot W_tukey
          beta_cv <- biweight_tukey_mest(x_cv, y_cv, beta_cv)
          
          # Hitung yhat_cv
          yhat_cv <- x[kota, ] %*% beta_cv
          
          # Hitung CV dan simpan nilai di CV(k)
          CV[k] <- CV[k] + sum((y[kota] - yhat_cv)^2)
        }
        }, error = function(e) {
          return()
        })
      }
      
      hasilCV <- cbind(l_awal, CV)
      A0 <- A
      B0 <- B
      minCV <- min(CV, na.rm = TRUE)    
      l_min <- which(CV == minCV)[1]
      
      if (l_min == 1) {
        A <- l_awal[l_min]
        B <- l_awal[l_min + 1]
      } else if (l_min == nl) {
        A <- l_awal[l_min - 1]
        B <- l_awal[l_min]
      } else {
        A <- l_awal[l_min - 1]
        B <- l_awal[l_min + 1]
      }
      
      selisih <- (B0 - A0) - (B - A)  # Selisih antara A0-B0 dan A-B
      iter_ab <- iter_ab + 1
      
      cat("Iterasi BW (Iterasi, Selisih):\n")
      cat(iter_ab, selisih, "\n")  # Tampilkan informasi iterasi
    }
    
    hasilCV <- hasilCV[order(hasilCV[, 2]), ] 
    l_ab <- rep(unname(hasilCV[1, 1]), times = n_area)[1:n_area]                
    CVmin <- rep(hasilCV[1, 2], times = n_area)[1:n_area]  
  
  return(l_ab)
}


## Fungsi Weight
cari_weight <- function(dist_mat, bw, n, n_waktu=1){
  bw <- rep(bw, each = n_waktu)
  # Inisialisasi matriks weight
  W_all <- matrix(0, nrow = n, ncol = n)  
  for (i in 1:n) {
    for (j in 1:n) {
      # Menghitung matriks weight berdasarkan kernel Gaussian
      W_all[i, j] <- exp(-(dist_mat[i, j]^2) / (2 * (bw[i]^2)))  
    }
  }
  return(W_all)
}


## Membuat Fungsi Koefisien Lokal
proses_lokal <- function(x, y, w, partest="t"){
  # Inisialisasi Koefisien Matrix & SE
  koef <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  var <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  se <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  # Inisialisasi Matrix Test
  t_values <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  wald_values <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  p_values <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  
  # Inisialisasi Matrix Residual
  residual <- matrix(0, nrow = nrow(x), ncol = 1)
  residual_lokal <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  residual_precentage <- matrix(0, nrow = nrow(x), ncol = 1)
  jkri <- matrix(0, nrow = nrow(x), ncol = 1)
  
  # Proses Perhitungan Koefisien & Test Setiap Observasi
  for (i in 1:nrow(x)) {
    # Ambil Weight Observasi
    weight <- diag(w[i,])
    
    # Perhitungan Koefisien Menggunakan Weight, Gagal digunakan invers general
    tryCatch({
      koef[i,] <- solve(t(x) %*% weight %*% x) %*% t(x) %*% weight %*% y
    }, error = function(e) {
      koef[i,] <- ginv(t(x) %*% weight %*% x) %*% t(x) %*% weight %*% y
    })
    
    # Perhitungan Komponen invers t(x)*w*x
    xtwx_inv <- tryCatch({
      solve(t(x) %*% weight %*% x)
    }, error = function(e) {
      ginv(t(x) %*% weight %*% x)
    })
    
    # Koefisien diiterasi kembali dgn M-EST. 
    xw <- t(t(x) %*% weight)
    yw <- t(t(y) %*% weight)
    koef[i,] <- biweight_tukey_mest(xw, yw, koef[i,]) 
    
    # Perhitungan Var & SE setiap model lokal
    residual_lokal[i, ] <- y - x %*% koef[i,]
    sse <- sum(residual_lokal[i,]^2 %*% weight)
    
    var[i, ] <- diag(xtwx_inv) * sse / (nrow(x) - ncol(x))
    se[i, ] <- sqrt(var[i, ])
    
    # Switch sesuai dengan jenis metode pemeriksaan parsial yang digunakan
    if (partest == "t"){
      t_values[i,] <- ifelse(se[i, ] == 0, 1e06, koef[i,] 
                             / se[i, ])
      p_values[i,] <- 2 * (1 - pt(abs(t_values[i,]), df=nrow(x)-ncol(x)))
    } else if (partest == 'wald'){
      wald_values[i,] <- ifelse(se[i, ] == 0, 1e06, (koef[i,]^2) 
                                / se[i, ]^2)
      p_values[i,] <- 1 - pchisq(wald_values[i,], df=1)
    }
    
    # Perhitungan Residual GWPR, Residual Precentage (MAPE), dan JKR tiap model
    residual[i, ] <- y[i, ] - x[i, ] %*% koef[i, ]
    residual_precentage[i, ] <- abs(residual[i, ]/y[i, ])
    jkri[i, ] <- ((x[i, ] %*% koef[i, ]) - mean(y))^2
  }
  
  # Perhitungan AIC
  jkr <- sum(jkri)
  jkt <- sum((y-mean(y))^2)
  jkg <- jkt-jkr
  aic <- (nrow(x) * log(jkg/nrow(x))) + 2*ncol(x)
  
  # Perhitungan Metrik Evaluasi Lain
  rsqr <- jkr/jkt
  mae <- sum(abs(residual))/nrow(x)
  mape <- sum(residual_precentage)*100/nrow(x)
  mse <- sum((residual)^2)/nrow(x)
  
  # Output yang berbeda sesuai dengan fungsi yang digunakan
  if (partest == 't'){
    listed <- list("coefficients" = koef, "t_values" = t_values, 
                   "p_values" = p_values, "residual" = residual, 
                   "var" = var, "se" = se, "AIC" = aic,
                   "Rsquare" = rsqr, "MAE" = mae, "MAPE" = mape, 
                   "MSE" = mse)
  } else {
    listed <- list("coefficients" = koef, "wald_values" = wald_values,
                   "p_values" = p_values, "residual" = residual, 
                   "var" = var, "se" = se, "AIC" = aic,
                   "Rsquare" = rsqr, "MAE" = mae, "MAPE" = mape, 
                   "MSE" = mse)
  }
  
  # Cek NaN & Throw Warning Jika Terdeteksi
  nan_inf_vars <- sapply(listed, function(x) {
    if (is.numeric(x)) {
      return(any(!is.finite(x)))
    } else {
      return(FALSE)
    }
  })
  
  nan_inf_vars_names <- names(nan_inf_vars)[nan_inf_vars]
  nan_inf_vars_names <- names(nan_inf_vars)[nan_inf_vars]
  
  if (length(nan_inf_vars_names) > 0) {
    warning(paste("Harap mengecek input kembali, Variabel berikut terdapat NaN atau Infinity:"
                  , paste(nan_inf_vars_names, collapse = ", ")))
  }
  
  # Return Variabel
  return(listed)
}

proses_global <- function(x, y, partest="t", simtest="F", B=1000){
  # Inisialisasi Matrix dan Vektor
  koef <- numeric(ncol(x))
  var <- 0
  se <- numeric(ncol(x))
  
  # Inisialisasi Test Matrix
  par_values <- numeric(ncol(x))
  p_values_par <- numeric(ncol(x))
  
  # Inisialisasi Vector Residual 
  residual <- numeric(nrow(x))
  
  # Hitung Koefisien
  koef <- solve(t(x) %*% x) %*% t(x) %*% y
  koef <- biweight_tukey_mest(x, y, koef) 
  
  # Hitung Residual & Komponen Lain
  residual <- y - x %*% koef
  var <- sum(residual^2) / (nrow(x) - ncol(x))
  se <- sqrt(diag(var * solve(t(x) %*% x)))
  
  # Switch untuk test Parsial
  if (partest == "t"){
    par = "Partial Test : T-Test"
    par_values <- ifelse(se == 0, 1e06, koef / se)
    p_values_par <- 2 * (1 - pt(abs(par_values), df=nrow(x)-ncol(x)))
  } else if (partest == 'wald'){
    par = "Partial Test : Wald-Test"
    par_values <- ifelse(se == 0, 1e06, (koef^2) / se^2)
    p_values_par <- 1 - pchisq(par_values, df=1)
  }
  
  # Switch untuk test Simultan
  if (simtest == "F"){
    # F-test
    sim = "Simultaneously Test : F-Test"
    f_value <- (sum((x %*% koef - mean(y))^2) / (ncol(x) - 1)) / var
    p_value_sim <- 1 - pf(f_value, df1=ncol(x)-1, df2=nrow(x)-ncol(x))
  } else if (simtest == "Bootstrap"){
    # Bootstrap test
    sim = "Simultaneously Test : Bootstrap Significance Test"
    bootstrap_values <- replicate(B, {
      index <- sample(1:nrow(x), nrow(x), replace = TRUE)
      x_b <- x[index, ]
      y_b <- y[index, ]
      koef_b <- solve(t(x_b) %*% x_b) %*% t(x_b) %*% y_b
      sum((x_b %*% koef_b - mean(y_b))^2) / (ncol(x_b) - 1)
    })
    p_value_sim <- mean(bootstrap_values >= var)
  }
  
  # AIC
  aic <- (nrow(x) * log(var)) + 2*ncol(x)
  
  # MAPE, MAE, MSE
  residual_percentage <- abs(residual/y)
  mape <- sum(residual_percentage)*100/nrow(x)
  mae <- sum(abs(residual))/nrow(x)
  mse <- sum((residual)^2)/nrow(x)
  
  # R-square
  ssr <- sum((x %*% koef - mean(y))^2)
  sst <- sum((y - mean(y))^2)
  rsqr <- ssr/sst
  
  # Return Variables
  listed <- list("coefficients" = koef, "partial_type" = par, 
                 "par_values" = par_values, "p_values_par" = p_values_par, 
                 "simultaneously_type" = sim, "p_value_sim" = p_value_sim,
                 "residual" = residual, "var" = var, "se" = se,
                 "AIC" = aic, "MAPE" = mape, "MAE" = mae, 
                 "MSE" = mse, "Rsquare" = rsqr
  )
  
  # Cek NaN & Throw Warning Jika Terdeteksi
  nan_inf_vars <- sapply(listed, function(x) {
    if (is.numeric(x)) {
      return(any(!is.finite(x)))
    } else {
      return(FALSE)
    }
  })
  
  nan_inf_vars_names <- names(nan_inf_vars)[nan_inf_vars]
  
  if (length(nan_inf_vars_names) > 0) {
    warning(paste("Harap mengecek input kembali, Variabel berikut terdapat NaN atau Infinity:"
                  , paste(nan_inf_vars_names, collapse = ", ")))
  }
  return(listed)
}

