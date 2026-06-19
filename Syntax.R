
# =========================
# 0) packages
# =========================
library(readxl)
library(dplyr)
library(tidyr)
library(cmdstanr)
library(posterior)
library(Matrix)
library(mvtnorm)
library(ggplot2)
library(bayesplot)
library(plotly)
library(knitr)
library(kableExtra)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

setwd("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/Syntax")

# =========================
# 1) fungsi
# =========================

# Triangle wide -> long
long_fun <- function(data){
  data %>%
    pivot_longer(
      cols = -c(`Treaty Year`, `Ultimate Premium`),
      names_to = "dev_year",
      values_to = "paid_cum"
    ) %>%
    mutate(dev_year = as.integer(dev_year)) %>%
    filter(!is.na(paid_cum))
}

# AY index (based on Treaty Year)
add_ay_index <- function(df){
  df %>%
    arrange(`Treaty Year`) %>%
    mutate(ay_index = as.integer(factor(`Treaty Year`, levels = sort(unique(`Treaty Year`)))))
}

# Statistical descriptives by dev_year
stat_des <- function(data){
  data %>%
    group_by(dev_year) %>%
    summarise(
      n     = n(),
      mean  = mean(paid_cum),
      sd    = sd(paid_cum),
      min   = min(paid_cum),
      max   = max(paid_cum),
      cv    = sd/mean,
      .groups = "drop"
    )
}

# z-score response
add_zscore <- function(df, y_col = "paid_cum"){
  mu <- mean(df[[y_col]])
  sdv <- sd(df[[y_col]])
  df <- df %>% mutate(z = ( .data[[y_col]] - mu ) / sdv)
  list(df = df, mu = mu, sd = sdv)
}

# min-max normalize input
add_minmax_xy <- function(df){
  min_ay <- min(df$ay_index); max_ay <- max(df$ay_index)
  min_dev <- min(df$dev_year); max_dev <- max(df$dev_year)
  
  df %>%
    mutate(
      ay_norm  = (ay_index - min_ay) / (max_ay - min_ay),
      dev_norm = (dev_year - min_dev) / (max_dev - min_dev)
    )
}

# Empirical semivariogram 2D (using paid_cum on original scale)
semivariogram_2d <- function(df, bin_width = 1, min_pairs = 10){
  stopifnot(all(c("ay_index", "dev_year", "paid_cum") %in% names(df)))
  
  coords <- as.matrix(df[, c("dev_year","ay_index")])
  y <- df$paid_cum
  
  D <- as.matrix(dist(coords))
  Ydiff2 <- (outer(y,y,"-"))^2
  
  idx <- which(upper.tri(D), arr.ind=TRUE)
  h <- D[idx]
  gamma <- 0.5 * Ydiff2[idx]
  
  breaks <- seq(0, max(h, na.rm = TRUE) + bin_width, by = bin_width)
  bin <- cut(h, breaks = breaks, include.lowest = TRUE, right = FALSE)
  
  out <- data.frame(h = h, gamma = gamma, bin = bin) %>%
    group_by(bin) %>%
    summarise(
      h_mid = mean(h),
      gamma_hat = mean(gamma),
      npairs = n(),
      .groups = "drop"
    ) %>%
    filter(npairs >= min_pairs)  
  
  out
}

# 3D surface helper from long -> matrix
to_surface_matrix <- function(df_long, value_col, ay_col = "ay_index", dev_col = "dev_year"){
  wide <- df_long %>%
    select(all_of(c(ay_col, dev_col, value_col))) %>%
    pivot_wider(names_from = all_of(dev_col), values_from = all_of(value_col)) %>%
    arrange(.data[[ay_col]])
  
  Z <- as.matrix(wide[, -1, drop = FALSE])
  list(wide = wide, Z = Z)
}


#3d vis
graph_3d<- function(data, name){
  dev_cols <- setdiff(names(data),c("Treaty Year","Ultimate Premium"))
  Z <- as.matrix(data[,dev_cols])
  storage.mode(Z) <- "double"
  n <- nrow(Z)
  x <- 1:(nrow(Z)+1)
  y <- 1:(ncol(Z)+1)
  p_upper <- plot_ly(
    x = x, y = y, z = Z,
    type = "surface"
  ) %>%
    layout(
      title = paste("Upper Triangle", name), 
      scene = list(
        xaxis = list(title = "DY"),
        yaxis = list(title = "AY"),
        zaxis = list(title = "Loss")
      )
    )
  return(p_upper)
}

# =================================================
# 2) Kernel Functions (untuk prediksi di R)
# =================================================
#SE kernel without warping
cov_SE <- function(X1, X2, psi, eta2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  K <- matrix(0, n1, n2)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (X1[i,1] - X2[j,1])^2 +
        psi[2] * (X1[i,2] - X2[j,2])^2
      K[i,j] <- eta2 * exp(-dsq)
    }
  }
  
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# SE kernel with warping
cov_SE_warp <- function(X1, X2, psi, eta2, a1, b1, a2, b2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  
  W1 <- cbind(pbeta(X1[,1], a1, b1), pbeta(X1[,2], a2, b2))
  W2 <- cbind(pbeta(X2[,1], a1, b1), pbeta(X2[,2], a2, b2))
  
  K <- matrix(0, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (W1[i,1] - W2[j,1])^2 +
        psi[2] * (W1[i,2] - W2[j,2])^2
      K[i,j] <- eta2 * exp(-dsq)
    }
  }
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# Matern 3/2 without warping
cov_Matern32 <- function(X1, X2, psi, eta2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  K <- matrix(0, n1, n2)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (X1[i,1] - X2[j,1])^2 +
        psi[2] * (X1[i,2] - X2[j,2])^2
      r <- sqrt(dsq)
      K[i,j] <- eta2 * (1 + sqrt(3) * r) * exp(-sqrt(3) * r)
    }
  }
  
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# Matern 3/2 with warping
cov_Matern32_warp <- function(X1, X2, psi, eta2, a1, b1, a2, b2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  
  W1 <- cbind(pbeta(X1[,1], a1, b1), pbeta(X1[,2], a2, b2))
  W2 <- cbind(pbeta(X2[,1], a1, b1), pbeta(X2[,2], a2, b2))
  
  K <- matrix(0, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (W1[i,1] - W2[j,1])^2 +
        psi[2] * (W1[i,2] - W2[j,2])^2
      r <- sqrt(dsq)
      K[i,j] <- eta2 * (1 + sqrt(3) * r) * exp(-sqrt(3) * r)
    }
  }
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# Matern 3/2 without warping
cov_Matern52 <- function(X1, X2, psi, eta2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  K <- matrix(0, n1, n2)
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (X1[i,1] - X2[j,1])^2 +
        psi[2] * (X1[i,2] - X2[j,2])^2
      r <- sqrt(dsq)
      K[i,j] <- eta2 *
        (1 + sqrt(5) * r + (5/3) * r^2) *
        exp(-sqrt(5) * r)
    }
  }
  
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# Matern 5/2 with warping
cov_Matern52_warp <- function(X1, X2, psi, eta2, a1, b1, a2, b2, sigma2 = 0) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  
  W1 <- cbind(pbeta(X1[,1], a1, b1), pbeta(X1[,2], a2, b2))
  W2 <- cbind(pbeta(X2[,1], a1, b1), pbeta(X2[,2], a2, b2))
  
  K <- matrix(0, n1, n2)
  for (i in 1:n1) {
    for (j in 1:n2) {
      dsq <- psi[1] * (W1[i,1] - W2[j,1])^2 +
        psi[2] * (W1[i,2] - W2[j,2])^2
      r <- sqrt(dsq)
      K[i,j] <- eta2 * (1 + sqrt(5) * r + (5/3) * r^2) * exp(-sqrt(5) * r)
    }
  }
  if (n1 == n2 && sigma2 > 0) K <- K + diag(sigma2, n1)
  K
}

# =========================
# 3) Baca Data
# =========================
df_prop <- read_excel(
  "D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/Data/data_fix.xlsx",
  sheet = "properti"
)

df_cred <- read_excel(
  "D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/Data/data_fix.xlsx",
  sheet = "cred"
)

df_marine <- read_excel(
  "D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/Data/data_fix.xlsx",
  sheet = "marine"
)




# =========================
# 4) Fungsi Preprocessing per LOB
# =========================
prep_lob <- function(df_raw){
  
  df_long <- long_fun(df_raw) %>%
    add_ay_index()
  
  zinfo <- add_zscore(df_long, "paid_cum")
  df_long <- zinfo$df
  
  # simpan min/max sebelum add_minmax_xy (pakai index asli)
  min_ay  <- min(df_long$ay_index);  max_ay  <- max(df_long$ay_index)
  min_dev <- min(df_long$dev_year);  max_dev <- max(df_long$dev_year)
  
  df_long <- df_long %>%
    mutate(
      ay_norm  = (ay_index - min_ay) / (max_ay - min_ay),
      dev_norm = (dev_year - min_dev) / (max_dev - min_dev)
    )
  
  list(
    data = df_long,
    mu = zinfo$mu,
    sd = zinfo$sd,
    min_ay = min_ay, max_ay = max_ay,
    min_dev = min_dev, max_dev = max_dev
  )
}

# Preprocess semua LOB
lob_prop <- prep_lob(df_prop)
lob_cred  <- prep_lob(df_cred)
lob_marine  <- prep_lob(df_marine)

lob_list <- list(prop = lob_prop, cred = lob_cred, marine = lob_marine)

lob_list



# =========================
# 3) 3D + Semivariogram
# =========================

View(lob_list$prop$data)

#=========================
# Properti 
#=========================
prop_3d <- to_surface_matrix(lob_list$prop$data, value_col = "paid_cum")
x_dev_p <- 1:ncol(prop_3d$Z)
y_ay_p  <- prop_3d$wide$ay_index

# copas trus isi NA sama nol
Z_clean_p <- prop_3d$Z
Z_clean_p[is.na(Z_clean_p)] <- 0

# mask biar si 12 nya muncul
mask_p <- matrix(1, nrow = nrow(Z_clean_p), ncol = ncol(Z_clean_p))
mask_p[is.na(prop_3d$Z)] <- 0  # 0 = transparent

viridis_rev <- list(
  c(0, "#FDE725"),  # kecil = kuning terang
  c(0.25, "#5DC863"),
  c(0.5, "#21918C"),
  c(0.75, "#3B528B"),
  c(1, "#440154")  # besar = ungu gelap
)


plot_ly(
  x = x_dev_p,
  y = y_ay_p,
  z = Z_clean_p,
  type = "surface",
  opacity = mask_p,  # Use mask to control transparency
  colorscale = viridis_rev
) %>%
  layout(
    title = "Permukaan 3D Properti (Upper Triangle Only)",
    scene = list(
      xaxis = list(
        title = "Development (dev_year)",
        tickmode = "array",
        tickvals = x_dev_p,
        ticktext = as.character(x_dev_p)
      ),
      yaxis = list(
        title = "Accident Year Index (ay_index)",
        tickmode = "array",
        tickvals = y_ay_p,
        ticktext = as.character(y_ay_p)
      ),
      zaxis = list(title = "paid_cum")
    )
  )

#Heatmap
plot_ly(
  x = x_dev_p,
  y = y_ay_p,
  z = prop_3d$Z,
  type = "heatmap",
  colorscale = viridis_rev
) %>%
  layout(
    title = "Heatmap Cumulative Paid Claim - Properti",
    xaxis = list(
      title = "Development Year (DY)",
      tickmode = "array",
      tickvals = 1:12,
      ticktext = as.character(1:12)
    ),
    yaxis = list(
      title = "Accident Year Index (AY)",
      tickmode = "linear",
      dtick = 1
    )
  )
prop_3d


ggplot(lob_list$prop$data,
       aes(x = dev_year,
           y = z,
           group = ay_index,
           color = factor(ay_index))) +
  geom_line(size = 1.1) +
  geom_point() +
  labs(
    title = "Pola Perkembangan Z-score Properti",
    x = "Development Year",
    y = "Klaim yang Distandarisasi (z)",
    color = "AY"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(2, 14, by = 2)
  )
#=========================
# Credit
#=========================
cred_3d <- to_surface_matrix(lob_list$cred$data, value_col = "paid_cum")
x_dev_c <- 1:ncol(cred_3d$Z)
y_ay_c <- cred_3d$wide$ay_index

# copas trus isi NA sama nol
Z_clean_c <- cred_3d$Z
Z_clean_c[is.na(Z_clean_c)] <- 0

# mask biar si 12 nya muncul
mask_c <- matrix(1, nrow = nrow(Z_clean_c), ncol = ncol(Z_clean_c))
mask_c[is.na(cred_3d$Z)] <- 0  # 0 = transparent

plot_ly(
  x = x_dev_c,
  y = y_ay_c,
  z = Z_clean_c,
  type = "surface",
  opacity = mask_c,  # Use mask to control transparency
  colorscale = viridis_rev
) %>%
  layout(
    title = "Permukaan 3D Credit (Upper Triangle Only)",
    scene = list(
      xaxis = list(
        title = "Development (dev_year)",
        tickmode = "array",
        tickvals = x_dev_c,
        ticktext = as.character(x_dev_c)
      ),
      yaxis = list(
        title = "Accident Year Index (ay_index)",
        tickmode = "array",
        tickvals = y_ay_c,
        ticktext = as.character(y_ay_c)
      ),
      zaxis = list(title = "paid_cum")
    )
  )

#Heatmap
plot_ly(
  x = x_dev_c,
  y = y_ay_c,
  z = cred_3d$Z,
  type = "heatmap",
  colorscale = viridis_rev
) %>%
  layout(
    title = "Heatmap Cumulative Paid Claim - Credit",
    xaxis = list(
      title = "Development Year (DY)",
      tickmode = "array",
      tickvals = 1:12,
      ticktext = as.character(1:12)
    ),
    yaxis = list(
      title = "Accident Year Index (AY)",
      tickmode = "linear",
      dtick = 1
    )
  )
cred_3d

ggplot(lob_list$cred$data,
       aes(x = dev_year,
           y = z,
           group = ay_index,
           color = factor(ay_index))) +
  geom_line(size = 1.1) +
  geom_point() +
  labs(
    title = "Pola Perkembangan Z-score Kredit",
    x = "Development Year",
    y = "Klaim yang Distandarisasi (z)",
    color = "AY"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(2, 14, by = 2)
  )
#=========================
# Marine
#=========================
mar_3d <- to_surface_matrix(lob_list$marine$data, value_col = "paid_cum")
x_dev_m <- 1:ncol(mar_3d$Z)
y_ay_m  <- mar_3d$wide$ay_index

# copas trus isi NA sama nol
Z_clean_m <- mar_3d$Z
Z_clean_m[is.na(Z_clean_m)] <- 0

# mask biar si 12 nya muncul
mask_m <- matrix(1, nrow = nrow(Z_clean_m), ncol = ncol(Z_clean_m))
mask_m[is.na(mar_3d$Z)] <- 0  # 0 = transparent

plot_ly(
  x = x_dev_m,
  y = y_ay_m,
  z = Z_clean_m,
  type = "surface",
  opacity = mask_m,  # Use mask to control transparency
  colorscale = viridis_rev
) %>%
  layout(
    title = "Permukaan 3D Marine (Upper Triangle Only)",
    scene = list(
      xaxis = list(
        title = "Development (dev_year)",
        tickmode = "array",
        tickvals = x_dev_m,
        ticktext = as.character(x_dev_m)
      ),
      yaxis = list(
        title = "Accident Year Index (ay_index)",
        tickmode = "array",
        tickvals = y_ay_m,
        ticktext = as.character(y_ay_m)
      ),
      zaxis = list(title = "paid_cum")
    )
  )

#Heatmap
plot_ly(
  x = x_dev_m,
  y = y_ay_m,
  z = mar_3d$Z,
  type = "heatmap",
  colorscale = viridis_rev
) %>%
  layout(
    title = "Heatmap Cumulative Paid Claim - Marine",
    xaxis = list(
      title = "Development Year (DY)",
      tickmode = "array",
      tickvals = 1:12,
      ticktext = as.character(1:12)
    ),
    yaxis = list(
      title = "Accident Year Index (AY)",
      tickmode = "linear",
      dtick = 1
    )
  )
mar_3d

ggplot(lob_list$marine$data,
       aes(x = dev_year,
           y = z,
           group = ay_index,
           color = factor(ay_index))) +
  geom_line(size = 1.1) +
  geom_point() +
  labs(
    title = "Pola Perkembangan Z-score Maritim",
    x = "Development Year",
    y = "Klaim yang Distandarisasi (z)",
    color = "AY"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(2, 14, by = 2)
  )
# Semivariogram (paid_cum)

#=========================
# Properti 
#=========================
prop_svg <- semivariogram_2d(lob_list$prop$data, bin_width = 1)
plot(prop_svg$h_mid, prop_svg$gamma_hat,
     type = "o",
     lwd  = 2.5,        
     pch  = 16,         
     cex  = 1.2,        
     xlab = "Jarak h",
     ylab = expression(hat(gamma)(h)),
     main = "Empirical Semivariogram Properti"
)

#=========================
# Credit
#=========================
cred_svg <- semivariogram_2d(lob_list$cred$data, bin_width = 1)
plot(cred_svg$h_mid, cred_svg$gamma_hat,
     type = "o",
     lwd  = 2.5,        
     pch  = 16,         
     cex  = 1.2,        
     xlab = "Jarak h",
     ylab = expression(hat(gamma)(h)),
     main = "Empirical Semivariogram Credit"
)

#=========================
# Marine
#=========================
mar_svg <- semivariogram_2d(lob_list$marine$data, bin_width = 1)
plot(mar_svg$h_mid, mar_svg$gamma_hat,
     type = "o",
     lwd  = 2.5,        
     pch  = 16,         
     cex  = 1.2,        
     xlab = "Jarak h",
     ylab = expression(hat(gamma)(h)),
     main = "Empirical Semivariogram Marine"
)

# library(writexl)
# 
# write_xlsx(
#   list(
#     properti_cln = lob_list$prop$data,
#     cred_cln     = lob_list$cred$data,
#     marine_cln   = lob_list$marine$data
#   ),
#   path = "D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/Data/data_fix.xlsx"
# )

View(lob_prop$data)


# =========================
# 5) Stack Data untuk Stan
# =========================
X_stack <- do.call(rbind, lapply(lob_list, \(x) as.matrix(x$data[,c("dev_norm","ay_norm")])))
z_stack <- unlist(lapply(lob_list, \(x) x$data$z))

N_l   <- sapply(lob_list, \(x) nrow(x$data))

start <- cumsum(c(1, head(N_l, -1)))
names(start) <- names(N_l)
start
start + N_l - 1

# =========================
# MODEL I tanpa Warp Estimasi Parameter
# =========================

# =========================
# Squared Exponential
# =========================
set_cmdstan_path("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/cmdstan-2.37.0")

mod1_se <- cmdstan_model("model_1_SE.stan", force_recompile = TRUE)

fit_lob_se <- function(lob_obj, lob_name,
                       iter_warmup = 1500,
                       iter_sampling = 1500,
                       seed = 123,
                       adapt_delta = 0.95,
                       max_treedepth = 15){
  
  stan_data <- list(
    N = nrow(lob_obj$data),
    X = as.matrix(lob_obj$data[, c("dev_norm", "ay_norm")]),
    z = as.vector(lob_obj$data$z)
  )
  
  fit <- mod1_se$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 1 SE -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = c("eta_sq", "psi", "sigma_sq"))
  
  summ <- fit$summary(
    variables = c("eta_sq", "psi", "sigma_sq")
  )[, c("mean","sd","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m1_prop_se   <- fit_lob_se(lob_list$prop,   "Properti")
m1_cred_se   <- fit_lob_se(lob_list$cred,   "Credit")
m1_marine_se <- fit_lob_se(lob_list$marine, "Marine")


# =========================
# Trace Plot
# =========================
library(posterior)
library(bayesplot)
color_scheme_set("blue")

plot_trace_se <- function(res_obj){
  draws <- res_obj$fit$draws(variables = c("eta_sq", "psi", "sigma_sq"))
  mcmc_trace(draws) + ggtitle(paste("Traceplot Model 1 SE -", res_obj$name))
}

plot_trace_se(m1_prop_se)
plot_trace_se(m1_cred_se)
plot_trace_se(m1_marine_se)

# =========================
# PPC
# =========================

ppc_se <- function(res_obj){
  z_obs <- res_obj$stan_data$z
  z_rep <- res_obj$fit$draws("z_rep", format = "matrix")
  
  print(
    ppc_dens_overlay(
      y = z_obs,
      yrep = z_rep[1:200, ]
    ) + ggtitle(paste("PPC Density Overlay Model 1 -", res_obj$name))
  )
  
  D_obs <- sum(z_obs^2)
  D_rep <- apply(z_rep, 1, function(x) sum(x^2))
  p_B <- mean(D_rep >= D_obs)
  
  cat("\nBayesian p-value -", res_obj$name, ":", p_B, "\n")
  
  print(ppc_stat(z_obs, z_rep, stat = "sd") + ggtitle(paste("PPC sd Model 1 -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "max") + ggtitle(paste("PPC max Model 1 -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "median") + ggtitle(paste("PPC median Model 1 -", res_obj$name)))
}

ppc_se(m1_prop_se)
ppc_se(m1_cred_se)
ppc_se(m1_marine_se)


m1_prop_se
# =========================
# Matern 3/2
# =========================
set_cmdstan_path("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/cmdstan-2.37.0")

mod1_m32 <- cmdstan_model("model_1_Matern32.stan", force_recompile = TRUE)

fit_lob_m32 <- function(lob_obj, lob_name,
                        iter_warmup = 1500,
                        iter_sampling = 1500,
                        seed = 123,
                        adapt_delta = 0.95,
                        max_treedepth = 15){
  
  stan_data <- list(
    N = nrow(lob_obj$data),
    X = as.matrix(lob_obj$data[, c("dev_norm", "ay_norm")]),
    z = as.vector(lob_obj$data$z)
  )
  
  fit <- mod1_m32$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 1 Matern 3/2 -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = c("eta_sq", "psi", "sigma_sq"))
  
  summ <- fit$summary(
    variables = c("eta_sq", "psi", "sigma_sq")
  )[, c("mean","sd","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m1_prop_m32   <- fit_lob_m32(lob_list$prop,   "Properti")
m1_cred_m32   <- fit_lob_m32(lob_list$cred,   "Credit")
m1_marine_m32 <- fit_lob_m32(lob_list$marine, "Marine")


# =========================
# Trace Plot
# =========================
library(posterior)
library(bayesplot)
color_scheme_set("blue")

plot_trace_m32 <- function(res_obj){
  draws <- res_obj$fit$draws(variables = c("eta_sq", "psi", "sigma_sq"))
  mcmc_trace(draws) + ggtitle(paste("Traceplot Model 1 Matern 3/2 -", res_obj$name))
}

plot_trace_m32(m1_prop_m32)
plot_trace_m32(m1_cred_m32)
plot_trace_m32(m1_marine_m32)

# =========================
# PPC
# =========================
ppc_m32 <- function(res_obj){
  z_obs <- res_obj$stan_data$z
  z_rep <- res_obj$fit$draws("z_rep", format = "matrix")
  
  print(
    ppc_dens_overlay(
      y = z_obs,
      yrep = z_rep[1:200, ]
    ) + ggtitle(paste("PPC Density Overlay Model 1 -", res_obj$name))
  )
  
  D_obs <- sum(z_obs^2)
  D_rep <- apply(z_rep, 1, function(x) sum(x^2))
  p_B <- mean(D_rep >= D_obs)
  
  cat("\nBayesian p-value -", res_obj$name, ":", p_B, "\n")
  
  print(ppc_stat(z_obs, z_rep, stat = "sd") + ggtitle(paste("PPC sd Model 1 -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "max") + ggtitle(paste("PPC max Model 1  -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "median") + ggtitle(paste("PPC median Model 1  -", res_obj$name)))
}

ppc_m32(m1_prop_m32)
ppc_m32(m1_cred_m32)
ppc_m32(m1_marine_m32)

# =========================
# Matern 5/2
# =========================
set_cmdstan_path("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/cmdstan-2.37.0")

mod1_m52 <- cmdstan_model("model_1_Matern52.stan", force_recompile = TRUE)

fit_lob_m52 <- function(lob_obj, lob_name,
                        iter_warmup = 1500,
                        iter_sampling = 1500,
                        seed = 123,
                        adapt_delta = 0.95,
                        max_treedepth = 15){
  
  stan_data <- list(
    N = nrow(lob_obj$data),
    X = as.matrix(lob_obj$data[, c("dev_norm", "ay_norm")]),
    z = as.vector(lob_obj$data$z)
  )
  
  fit <- mod1_m52$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 1 Matern 5/2 -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = c("eta_sq", "psi", "sigma_sq"))
  
  summ <- fit$summary(
    variables = c("eta_sq", "psi", "sigma_sq")
  )[, c("mean","sd","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m1_prop_m52   <- fit_lob_m52(lob_list$prop,   "Properti")
m1_cred_m52   <- fit_lob_m52(lob_list$cred,   "Credit")
m1_marine_m52 <- fit_lob_m52(lob_list$marine, "Marine")

# =========================
# Trace Plot
# =========================

library(posterior)
library(bayesplot)
color_scheme_set("blue")

plot_trace_m52 <- function(res_obj){
  draws <- res_obj$fit$draws(variables = c("eta_sq", "psi", "sigma_sq"))
  mcmc_trace(draws) + ggtitle(paste("Traceplot Model 1 Matern 5/2 -", res_obj$name))
}

plot_trace_m52(m1_prop_m52)
plot_trace_m52(m1_cred_m52)
plot_trace_m52(m1_marine_m52)

# =========================
# PPC
# =========================
ppc_m52 <- function(res_obj){
  z_obs <- res_obj$stan_data$z
  z_rep <- res_obj$fit$draws("z_rep", format = "matrix")
  
  print(
    ppc_dens_overlay(
      y = z_obs,
      yrep = z_rep[1:200, ]
    ) + ggtitle(paste("PPC Density Overlay Model 1 -", res_obj$name))
  )
  
  D_obs <- sum(z_obs^2)
  D_rep <- apply(z_rep, 1, function(x) sum(x^2))
  p_B <- mean(D_rep >= D_obs)
  
  cat("\nBayesian p-value -", res_obj$name, ":", p_B, "\n")
  
  print(ppc_stat(z_obs, z_rep, stat = "sd") + ggtitle(paste("PPC sd Model 1 -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "max") + ggtitle(paste("PPC max Model 1 -", res_obj$name)))
  print(ppc_stat(z_obs, z_rep, stat = "median") + ggtitle(paste("PPC median Model 1 -", res_obj$name)))
}

ppc_m52(m1_prop_m52)
ppc_m52(m1_cred_m52)
ppc_m52(m1_marine_m52)


# =========================
# MODEL 2 dengan Warp Estimasi Parameter
# =========================

# =========================================================
# MODEL 2 = GP + Input Warping (Beta CDF)
# Kernel:
# 1. Squared Exponential
# 2. Matern 3/2
# 3. Matern 5/2
# =========================================================

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

set_cmdstan_path("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/cmdstan-2.37.0")
color_scheme_set("blue")

# =========================================================
# Helper umum Model 2
# =========================================================

# summary parameter model 2
summary_vars_m2 <- c("eta_sq", "psi", "sigma_sq", "a1", "b1", "a2", "b2")

# buat stan_data dari object LoB
make_stan_data_m2 <- function(lob_obj){
  list(
    N = nrow(lob_obj$data),
    X = as.matrix(lob_obj$data[, c("dev_norm", "ay_norm")]),
    z = as.vector(lob_obj$data$z)
  )
}

# helper traceplot
plot_trace_m2 <- function(res_obj, model_name){
  draws <- res_obj$fit$draws(
    variables = c("eta_sq", "psi", "sigma_sq", "a1", "b1", "a2", "b2")
  )
  
  mcmc_trace(draws) +
    ggtitle(paste("Traceplot", model_name, "-", res_obj$name))
}

# helper PPC
ppc_m2 <- function(res_obj, model_name){
  z_obs <- res_obj$stan_data$z
  z_rep <- res_obj$fit$draws("z_rep", format = "matrix")
  
  print(
    ppc_dens_overlay(
      y = z_obs,
      yrep = z_rep[1:200, ]
    ) + ggtitle(paste("PPC Density Overlay", model_name, "-", res_obj$name))
  )
  
  D_obs <- sum(z_obs^2)
  D_rep <- apply(z_rep, 1, function(x) sum(x^2))
  p_B <- mean(D_rep >= D_obs)
  
  cat("\nBayesian p-value", model_name, "-", res_obj$name, ":", p_B, "\n")
  
  print(
    ppc_stat(z_obs, z_rep, stat = "sd") +
      ggtitle(paste("PPC sd", model_name, "-", res_obj$name))
  )
  
  print(
    ppc_stat(z_obs, z_rep, stat = "max") +
      ggtitle(paste("PPC max", model_name, "-", res_obj$name))
  )
  
  print(
    ppc_stat(z_obs, z_rep, stat = "median") +
      ggtitle(paste("PPC median", model_name, "-", res_obj$name))
  )
}

# =========================================================
# 1) MODEL 2 - Squared Exponential Warp
# =========================================================

mod2_se <- cmdstan_model("model_2_SE_warp.stan", force_recompile = TRUE)

fit_lob_m2_se <- function(lob_obj, lob_name,
                          iter_warmup = 1500,
                          iter_sampling = 1500,
                          seed = 123,
                          adapt_delta = 0.95,
                          max_treedepth = 15){
  
  stan_data <- make_stan_data_m2(lob_obj)
  
  fit <- mod2_se$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 2 SE Warp -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = summary_vars_m2)
  
  summ <- fit$summary(
    variables = summary_vars_m2
  )[, c("mean","median","sd","mad","q5","q95","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m2_prop_se   <- fit_lob_m2_se(lob_list$prop,   "Properti")
m2_cred_se   <- fit_lob_m2_se(lob_list$cred,   "Credit")
m2_marine_se <- fit_lob_m2_se(lob_list$marine, "Marine")

# Trace Plot
plot_trace_m2(m2_prop_se,   "Model 2 SE Warp")
plot_trace_m2(m2_cred_se,   "Model 2 SE Warp")
plot_trace_m2(m2_marine_se, "Model 2 SE Warp")

# PPC
ppc_m2(m2_prop_se,   "Model 2 SE Warp")
ppc_m2(m2_cred_se,   "Model 2 SE Warp")
ppc_m2(m2_marine_se, "Model 2 SE Warp")


# =========================================================
# 2) MODEL 2 - Matern 3/2 Warp
# =========================================================

mod2_m32 <- cmdstan_model("model_2_Matern32_warp.stan", force_recompile = TRUE)

fit_lob_m2_m32 <- function(lob_obj, lob_name,
                           iter_warmup = 1500,
                           iter_sampling = 1500,
                           seed = 123,
                           adapt_delta = 0.95,
                           max_treedepth = 15){
  
  stan_data <- make_stan_data_m2(lob_obj)
  
  fit <- mod2_m32$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 2 Matern 3/2 Warp -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = summary_vars_m2)
  
  summ <- fit$summary(
    variables = summary_vars_m2
  )[, c("mean","median","sd","mad","q5","q95","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m2_prop_m32   <- fit_lob_m2_m32(lob_list$prop,   "Properti")
m2_cred_m32   <- fit_lob_m2_m32(lob_list$cred,   "Credit")
m2_marine_m32 <- fit_lob_m2_m32(lob_list$marine, "Marine")

# Trace Plot
plot_trace_m2(m2_prop_m32,   "Model 2 Matern 3/2 Warp")
plot_trace_m2(m2_cred_m32,   "Model 2 Matern 3/2 Warp")
plot_trace_m2(m2_marine_m32, "Model 2 Matern 3/2 Warp")

# PPC
ppc_m2(m2_prop_m32,   "Model 2 Matern 3/2 Warp")
ppc_m2(m2_cred_m32,   "Model 2 Matern 3/2 Warp")
ppc_m2(m2_marine_m32, "Model 2 Matern 3/2 Warp")


# =========================================================
# 3) MODEL 2 - Matern 5/2 Warp
# =========================================================

mod2_m52 <- cmdstan_model("model_2_Matern52_warp.stan", force_recompile = TRUE)

fit_lob_m2_m52 <- function(lob_obj, lob_name,
                           iter_warmup = 1500,
                           iter_sampling = 1500,
                           seed = 123,
                           adapt_delta = 0.95,
                           max_treedepth = 15){
  
  stan_data <- make_stan_data_m2(lob_obj)
  
  fit <- mod2_m52$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  cat("\n=============================\n")
  cat("Hasil Model 2 Matern 5/2 Warp -", lob_name, "\n")
  cat("=============================\n")
  
  print(fit, variables = summary_vars_m2)
  
  summ <- fit$summary(
    variables = summary_vars_m2
  )[, c("mean","median","sd","mad","q5","q95","rhat","ess_bulk","ess_tail")]
  
  print(summ)
  
  list(
    name = lob_name,
    stan_data = stan_data,
    fit = fit,
    summary = summ
  )
}

# Fit per LoB
m2_prop_m52   <- fit_lob_m2_m52(lob_list$prop,   "Properti")
m2_cred_m52   <- fit_lob_m2_m52(lob_list$cred,   "Credit")
m2_marine_m52 <- fit_lob_m2_m52(lob_list$marine, "Marine")

# Trace Plot
plot_trace_m2(m2_prop_m52,   "Model 2 Matern 5/2 Warp")
plot_trace_m2(m2_cred_m52,   "Model 2 Matern 5/2 Warp")
plot_trace_m2(m2_marine_m52, "Model 2 Matern 5/2 Warp")

# PPC
ppc_m2(m2_prop_m52,   "Model 2 Matern 5/2 Warp")
ppc_m2(m2_cred_m52,   "Model 2 Matern 5/2 Warp")
ppc_m2(m2_marine_m52, "Model 2 Matern 5/2 Warp")



haSIL
m2_prop_se$summary

# =========================
# MODEL 3 Hierarchical GP + Input Warping
# Squared Exponential Kernel
# =========================

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

set_cmdstan_path("D:/Kuliah/TA/IFRS 17 Bayesian hierarhcial/cmdstan-2.37.0")
color_scheme_set("blue")

# =========================================================
# MODEL 3 HELPER (PAKAI SEKALI AJA)
# =========================================================

summary_vars_m3_global <- c(
  "mu_psi1","mu_psi2","mu_eta","mu_sigma",
  "tau_psi1","tau_psi2","tau_eta","tau_sigma"
)

summary_vars_m3_local <- c(
  "log_psi1[1]","log_psi1[2]","log_psi1[3]",
  "log_psi2[1]","log_psi2[2]","log_psi2[3]",
  "log_eta_sq[1]","log_eta_sq[2]","log_eta_sq[3]",
  "log_sigma_sq[1]","log_sigma_sq[2]","log_sigma_sq[3]",
  "a1[1]","a1[2]","a1[3]",
  "b1[1]","b1[2]","b1[3]",
  "a2[1]","a2[2]","a2[3]",
  "b2[1]","b2[2]","b2[3]"
)

summary_vars_m3_nc_latent <- c(
  "z_psi1[1]","z_psi1[2]","z_psi1[3]",
  "z_psi2[1]","z_psi2[2]","z_psi2[3]",
  "z_eta[1]","z_eta[2]","z_eta[3]",
  "z_sigma[1]","z_sigma[2]","z_sigma[3]"
)

make_stan_data_m3 <- function(lob_list){
  
  X_stack <- do.call(
    rbind,
    lapply(
      lob_list,
      function(x)
        as.matrix(
          x$data[,c("dev_norm","ay_norm")]
        )
    )
  )
  
  z_stack <- unlist(
    lapply(
      lob_list,
      function(x)x$data$z
    )
  )
  
  N_l <- sapply(
    lob_list,
    function(x)nrow(x$data)
  )
  
  start <- cumsum(
    c(1,head(N_l,-1))
  )
  
  list(
    L=length(lob_list),
    N_total=nrow(X_stack),
    N_l=as.integer(N_l),
    start=as.integer(start),
    X=X_stack,
    z=as.vector(z_stack)
  )
  
}

# =========================================================
# HELPER PPC MODEL 3
# =========================================================

split_indices_m3 <- function(stan_data){
  
  lapply(
    seq_along(stan_data$N_l),
    function(l){
      
      s <- stan_data$start[l]
      n <- stan_data$N_l[l]
      
      s:(s+n-1)
      
    }
  )
  
}


ppc_m3 <- function(
    res_obj,
    lob_names=c(
      "Properti",
      "Credit",
      "Marine"
    ),
    model_name=NULL){
  
  if(is.null(model_name)){
    model_name <- res_obj$name
  }
  
  z_obs <- res_obj$stan_data$z
  
  z_rep <- res_obj$fit$draws(
    "z_rep",
    format="matrix"
  )
  
  # =================================================
  # GLOBAL PPC
  # =================================================
  
  print(
    
    ppc_dens_overlay(
      y=z_obs,
      yrep=z_rep[1:200,]
    )+
      
      ggtitle(
        paste(
          "PPC Density -",
          model_name,
          "Global"
        )
      )
    
  )
  
  D_obs <- sum(
    z_obs^2
  )
  
  D_rep <- apply(
    z_rep,
    1,
    function(x)
      sum(x^2)
  )
  
  p_B <- mean(
    D_rep>=D_obs
  )
  
  cat(
    "\nBayesian p-value:",
    model_name,
    "(Global) =",
    round(p_B,4),
    "\n"
  )
  
  print(
    
    ppc_stat(
      z_obs,
      z_rep,
      stat="sd"
    )+
      
      ggtitle(
        paste(
          "PPC SD -",
          model_name
        )
      )
    
  )
  
  print(
    
    ppc_stat(
      z_obs,
      z_rep,
      stat="max"
    )+
      
      ggtitle(
        paste(
          "PPC MAX -",
          model_name
        )
      )
    
  )
  
  print(
    
    ppc_stat(
      z_obs,
      z_rep,
      stat="median"
    )+
      
      ggtitle(
        paste(
          "PPC Median -",
          model_name
        )
      )
    
  )
  
  # =================================================
  # PER LOB
  # =================================================
  
  idx_list <- split_indices_m3(
    res_obj$stan_data
  )
  
  for(i in seq_along(idx_list)){
    
    idx <- idx_list[[i]]
    
    y_i <- z_obs[idx]
    
    yrep_i <- z_rep[
      1:200,
      idx
    ]
    
    print(
      
      ppc_dens_overlay(
        y=y_i,
        yrep=yrep_i
      )+
        
        ggtitle(
          paste(
            "PPC Density -",
            model_name,
            lob_names[i]
          )
        )
      
    )
    
    D_obs_i <- sum(
      y_i^2
    )
    
    D_rep_i <- apply(
      yrep_i,
      1,
      function(x)
        sum(x^2)
    )
    
    p_B_i <- mean(
      D_rep_i>=D_obs_i
    )
    
    cat(
      "\nBayesian p-value:",
      model_name,
      lob_names[i],
      "=",
      round(
        p_B_i,
        4
      ),
      "\n"
    )
    
    print(
      
      ppc_stat(
        y_i,
        yrep_i,
        stat="sd"
      )+
        
        ggtitle(
          paste(
            "PPC SD -",
            lob_names[i]
          )
        )
      
    )
    
    print(
      
      ppc_stat(
        y_i,
        yrep_i,
        stat="max"
      )+
        
        ggtitle(
          paste(
            "PPC MAX -",
            lob_names[i]
          )
        )
      
    )
    
    print(
      
      ppc_stat(
        y_i,
        yrep_i,
        stat="median"
      )+
        
        ggtitle(
          paste(
            "PPC Median -",
            lob_names[i]
          )
        )
      
    )
    
  }
  
}

# =========================================================
# GENERIC FIT FUNCTION
# =========================================================

fit_model3 <- function(
    model_object,
    model_name,
    lob_list,
    iter_warmup=1500,
    iter_sampling=1500,
    seed=123,
    adapt_delta=0.95,
    max_treedepth=15){
  
  stan_data <- make_stan_data_m3(
    lob_list
  )
  
  fit <- model_object$sample(
    data=stan_data,
    chains=4,
    parallel_chains=4,
    iter_warmup=iter_warmup,
    iter_sampling=iter_sampling,
    seed=seed,
    adapt_delta=adapt_delta,
    max_treedepth=max_treedepth
  )
  
  cat("\n===================\n")
  cat(model_name,"\n")
  cat("===================\n")
  
  print(
    fit,
    variables=summary_vars_m3_global
  )
  
  summ_global <- fit$summary(
    variables=summary_vars_m3_global
  )[,c(
    "mean","median",
    "sd","mad",
    "q5","q95",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )]
  
  summ_local <- fit$summary(
    variables=summary_vars_m3_local
  )[,c(
    "mean","median",
    "sd","mad",
    "q5","q95",
    "rhat",
    "ess_bulk",
    "ess_tail"
  )]
  
  diag_sum <- fit$diagnostic_summary()
  
  list(
    name=model_name,
    stan_data=stan_data,
    fit=fit,
    summary_global=summ_global,
    summary_local=summ_local,
    diagnostic=diag_sum
  )
  
}

# =========================================================
# COMPILE MODEL
# =========================================================

# SE
mod3_se <- cmdstan_model(
  "model_3_SE_hier_warp2.stan",
  force_recompile=TRUE
)

mod3_se_nc <- cmdstan_model(
  "model_3_SE_hier_warp2_nc.stan",
  force_recompile=TRUE
)

# M32
mod3_m32 <- cmdstan_model(
  "model_3_Matern32_hier_warp.stan",
  force_recompile=TRUE
)

mod3_m32_nc <- cmdstan_model(
  "model_3_Matern32_hier_warp_nc.stan",
  force_recompile=TRUE
)

# M52
mod3_m52 <- cmdstan_model(
  "model_3_Matern52_hier_warp.stan",
  force_recompile=TRUE
)

mod3_m52_nc <- cmdstan_model(
  "model_3_Matern52_hier_warp_nc.stan",
  force_recompile=TRUE
)

# =========================================================
# RUN CENTERED
# =========================================================

m3_se <- fit_model3(
  mod3_se,
  "Model 3 SE Centered",
  lob_list
)

m3_m32 <- fit_model3(
  mod3_m32,
  "Model 3 Matern32 Centered",
  lob_list
)

m3_m52 <- fit_model3(
  mod3_m52,
  "Model 3 Matern52 Centered",
  lob_list
)

# =========================================================
# RUN NON CENTERED
# =========================================================

m3_se_nc <- fit_model3(
  mod3_se_nc,
  "Model 3 SE NonCentered",
  lob_list
)

m3_m32_nc <- fit_model3(
  mod3_m32_nc,
  "Model 3 Matern32 NonCentered",
  lob_list
)

m3_m52_nc <- fit_model3(
  mod3_m52_nc,
  "Model 3 Matern52 NonCentered",
  lob_list
)

# =========================================================
# QUICK CHECK
# =========================================================

model_list <- list(
  
  m3_se,
  m3_se_nc,
  
  m3_m32,
  m3_m32_nc,
  
  m3_m52,
  m3_m52_nc
  
)

for(obj in model_list){
  
  cat("\n====================\n")
  cat(obj$name)
  cat("\n====================\n")
  
  print(
    obj$fit,
    variables=summary_vars_m3_global
  )
  
  print(
    obj$fit,
    variables=summary_vars_m3_local,
    max_rows=1000
  )
  
  print(
    obj$diagnostic
  )
  
}

# NC latent only
print(
  m3_se_nc$fit,
  variables=summary_vars_m3_nc_latent
)

print(
  m3_m32_nc$fit,
  variables=summary_vars_m3_nc_latent
)

print(
  m3_m52_nc$fit,
  variables=summary_vars_m3_nc_latent
)

# =========================================================
# PPC CENTERED
# =========================================================

ppc_m3(m3_se)
ppc_m3(m3_m32)
ppc_m3(m3_m52)
m3_se$diagnostic
m3_m32$diagnostic
m3_m52$diagnostic

# =========================================================
# PPC NON CENTERED
# =========================================================

ppc_m3(m3_se_nc)
ppc_m3(m3_m32_nc)
ppc_m3(m3_m52_nc)

m3_se_nc$diagnostic
m3_m32_nc$diagnostic
m3_m52_nc$diagnostic

# =========================================================
# Helper: extract posterior draws for Model 3 SE final
# lob_index: 1=prop, 2=cred, 3=marine
# =========================================================
extract_draws_m3_se <- function(res_obj, lob_index, M = 1000, seed = 123){
  set.seed(seed)
  
  vars_need <- c(
    paste0("log_psi1[", lob_index, "]"),
    paste0("log_psi2[", lob_index, "]"),
    paste0("log_eta_sq[", lob_index, "]"),
    paste0("log_sigma_sq[", lob_index, "]"),
    paste0("a1[", lob_index, "]"),
    paste0("b1[", lob_index, "]"),
    paste0("a2[", lob_index, "]"),
    paste0("b2[", lob_index, "]")
  )
  
  ddf <- as.data.frame(res_obj$fit$draws(variables = vars_need, format = "draws_df"))
  
  keep <- sample(seq_len(nrow(ddf)), size = min(M, nrow(ddf)), replace = FALSE)
  ddf <- ddf[keep, ]
  
  out <- data.frame(
    psi1    = exp(ddf[[paste0("log_psi1[", lob_index, "]")]]),
    psi2    = exp(ddf[[paste0("log_psi2[", lob_index, "]")]]),
    eta2    = exp(ddf[[paste0("log_eta_sq[", lob_index, "]")]]),
    sigma2  = exp(ddf[[paste0("log_sigma_sq[", lob_index, "]")]]),
    a1      = ddf[[paste0("a1[", lob_index, "]")]],
    b1      = ddf[[paste0("b1[", lob_index, "]")]],
    a2      = ddf[[paste0("a2[", lob_index, "]")]],
    b2      = ddf[[paste0("b2[", lob_index, "]")]]
  )
  
  out
}



# =========================================================
# Build full grid and identify future cells
# =========================================================
build_pred_map <- function(lob_obj){
  df <- lob_obj$data
  
  ay_all  <- sort(unique(df$ay_index))
  dev_all <- sort(unique(df$dev_year))
  
  full_grid <- tidyr::expand_grid(
    ay_index = ay_all,
    dev_year = dev_all
  )
  
  obs_keys <- df %>%
    distinct(ay_index, dev_year) %>%
    mutate(observed = TRUE)
  
  grid_tagged <- full_grid %>%
    left_join(obs_keys, by = c("ay_index", "dev_year")) %>%
    mutate(observed = ifelse(is.na(observed), FALSE, TRUE))
  
  min_ay  <- min(ay_all);  max_ay  <- max(ay_all)
  min_dev <- min(dev_all); max_dev <- max(dev_all)
  
  grid_tagged <- grid_tagged %>%
    mutate(
      ay_norm  = (ay_index - min_ay) / (max_ay - min_ay),
      dev_norm = (dev_year - min_dev) / (max_dev - min_dev)
    )
  
  pred_map <- grid_tagged %>%
    filter(!observed) %>%
    arrange(ay_index, dev_year) %>%
    mutate(cell = row_number())
  
  list(
    full_grid = grid_tagged,
    pred_map = pred_map
  )
}


# =========================================================
# Posterior predictive lower triangle for Model 3 SE
# Uses cov_SE_warp() from your existing script
# =========================================================
simulate_lower_triangle_m3_se <- function(res_obj,
                                          lob_obj,
                                          lob_index,
                                          M = 1000,
                                          jitter = 1e-8,
                                          seed = 123,
                                          show_progress = TRUE){
  set.seed(seed)
  
  # posterior draws parameter
  post_draws <- extract_draws_m3_se(res_obj, lob_index = lob_index, M = M, seed = seed)
  M_use <- nrow(post_draws)
  
  # observed data
  obs_df <- lob_obj$data %>%
    arrange(ay_index, dev_year)
  
  X_obs <- as.matrix(obs_df[, c("dev_norm", "ay_norm")])
  z_obs <- obs_df$z
  
  # prediction map
  pm <- build_pred_map(lob_obj)
  pred_map <- pm$pred_map
  X_pred <- as.matrix(pred_map[, c("dev_norm", "ay_norm")])
  
  Npred <- nrow(X_pred)
  z_pred <- matrix(NA_real_, nrow = M_use, ncol = Npred)
  
  if(show_progress){
    pb <- txtProgressBar(min = 0, max = M_use, style = 3)
  }
  
  for(m in 1:M_use){
    psi   <- c(post_draws$psi1[m], post_draws$psi2[m])
    eta2  <- post_draws$eta2[m]
    sig2  <- post_draws$sigma2[m]
    a1    <- post_draws$a1[m]
    b1    <- post_draws$b1[m]
    a2    <- post_draws$a2[m]
    b2    <- post_draws$b2[m]
    
    K_oo <- cov_SE_warp(X_obs,  X_obs,  psi, eta2, a1, b1, a2, b2, sigma2 = sig2)
    K_op <- cov_SE_warp(X_obs,  X_pred, psi, eta2, a1, b1, a2, b2, sigma2 = 0)
    K_pp <- cov_SE_warp(X_pred, X_pred, psi, eta2, a1, b1, a2, b2, sigma2 = 0)
    
    K_oo <- K_oo + diag(jitter, nrow(K_oo))
    
    R <- chol(K_oo)
    
    # alpha = K_oo^{-1} z_obs
    alpha <- backsolve(R, forwardsolve(t(R), z_obs))
    
    # predictive mean
    mu_pred <- t(K_op) %*% alpha
    
    # predictive covariance
    v <- forwardsolve(t(R), K_op)
    Sigma_pred <- K_pp - crossprod(v)
    
    # tambahkan observation/process noise
    Sigma_obs <- Sigma_pred + diag(sig2, nrow(Sigma_pred))
    
    # stabilize
    Sigma_obs <- (Sigma_obs + t(Sigma_obs)) / 2
    Sigma_obs <- Sigma_obs + diag(jitter, nrow(Sigma_obs))
    
    # posterior predictive sampling
    z_pred[m, ] <- mvtnorm::rmvnorm(
      1,
      mean  = as.numeric(mu_pred),
      sigma = Sigma_obs
    )
    
    if(show_progress) setTxtProgressBar(pb, m)
  }
  
  if(show_progress) close(pb)
  
  # back-transform ke skala asli
  y_pred <- z_pred * lob_obj$sd + lob_obj$mu
  
  list(
    lob_index = lob_index,
    post_draws = post_draws,
    obs_df = obs_df,
    pred_map = pred_map,
    z_pred = z_pred,
    y_pred = y_pred,
    mu_y = lob_obj$mu,
    sd_y = lob_obj$sd
  )
}




# =========================================================
# Build full triangle per draw and enforce monotone cumulative
# =========================================================
build_full_triangle_draws <- function(sim_obj, lob_name = "LOB"){
  
  obs_long <- sim_obj$obs_df %>%
    select(ay_index, dev_year, paid_cum) %>%
    arrange(ay_index, dev_year)
  
  M_use <- nrow(sim_obj$y_pred)
  
  # observed replicated for each draw
  obs_rep <- tidyr::crossing(
    draw = 1:M_use,
    obs_long
  ) %>%
    mutate(
      paid_cum_raw = paid_cum,
      is_future = FALSE
    ) %>%
    select(draw, ay_index, dev_year, paid_cum_raw, is_future)
  
  # predicted future cells
  pred_long <- as.data.frame(sim_obj$y_pred)
  pred_long$draw <- 1:nrow(pred_long)
  
  pred_long <- tidyr::pivot_longer(
    pred_long,
    cols = -draw,
    names_to = "cell",
    values_to = "paid_cum_raw"
  ) %>%
    mutate(cell = as.integer(sub("^V", "", cell))) %>%
    left_join(sim_obj$pred_map %>% select(cell, ay_index, dev_year), by = "cell") %>%
    mutate(is_future = TRUE) %>%
    select(draw, ay_index, dev_year, paid_cum_raw, is_future)
  
  full_triangle <- bind_rows(obs_rep, pred_long) %>%
    arrange(draw, ay_index, dev_year) %>%
    group_by(draw, ay_index) %>%
    mutate(
      # enforce cumulative monotone
      paid_cum = cummax(paid_cum_raw),
      # derive incremental from full cumulative path
      incr = paid_cum - lag(paid_cum, default = 0)
    ) %>%
    ungroup() %>%
    mutate(lob = lob_name)
  
  full_triangle
}



# =========================================================
# Prepare discount curve
# =========================================================
make_discount_curve <- function(max_horizon, discount_rates = 0){
  if(length(discount_rates) == 1){
    rate_vec <- rep(discount_rates, max_horizon)
  } else {
    if(length(discount_rates) < max_horizon){
      stop("Length discount_rates kurang dari max_horizon")
    }
    rate_vec <- discount_rates[1:max_horizon]
  }
  
  data.frame(
    time_to_payment = 1:max_horizon,
    rate = rate_vec
  )
}

# =========================================================
# Compute BEL / RA / LIC from full triangle draws
# conf = confidence level, e.g. 0.75
# =========================================================
compute_bel_ra_lic <- function(full_triangle_draws,
                               conf = 0.75,
                               discount_rates = 0){
  
  max_ay <- max(full_triangle_draws$ay_index)
  
  # valuation diagonal for square triangle:
  # observed cells satisfy ay_index + dev_year <= max_ay + 1
  valuation_cutoff <- max_ay + 1
  
  tmp <- full_triangle_draws %>%
    mutate(
      time_to_payment = pmax(ay_index + dev_year - valuation_cutoff, 0)
    )
  
  max_horizon <- max(tmp$time_to_payment)
  
  if(max_horizon == 0){
    stop("Tidak ada future cell yang teridentifikasi.")
  }
  
  disc_tbl <- make_discount_curve(max_horizon = max_horizon,
                                  discount_rates = discount_rates)
  
  tmp <- tmp %>%
    left_join(disc_tbl, by = "time_to_payment") %>%
    mutate(
      rate = ifelse(is.na(rate), 0, rate),
      disc_factor = ifelse(time_to_payment == 0, 1, (1 + rate)^time_to_payment),
      incr_future = ifelse(is_future, incr, 0),
      pv_future = incr_future / disc_factor
    )
  
  # BEL per AY per draw
  bel_ay_draw <- tmp %>%
    filter(is_future) %>%
    group_by(draw, lob, ay_index) %>%
    summarise(
      BEL_AY = sum(pv_future, na.rm = TRUE),
      .groups = "drop"
    )
  
  # total BEL per draw
  bel_draw <- bel_ay_draw %>%
    group_by(draw, lob) %>%
    summarise(
      BEL = sum(BEL_AY, na.rm = TRUE),
      .groups = "drop"
    )
  
  # point estimate + RA
  q_conf <- as.numeric(quantile(bel_draw$BEL, probs = conf, na.rm = TRUE))
  mean_bel <- mean(bel_draw$BEL, na.rm = TRUE)
  med_bel  <- median(bel_draw$BEL, na.rm = TRUE)
  
  RA_VaR <- q_conf - mean_bel
  RA_CTE <- mean(bel_draw$BEL[bel_draw$BEL >= q_conf], na.rm = TRUE) - mean_bel
  
  summary_tbl <- data.frame(
    lob = unique(bel_draw$lob),
    BEL_mean = mean_bel,
    BEL_median = med_bel,
    BEL_sd = sd(bel_draw$BEL, na.rm = TRUE),
    BEL_q_conf = q_conf,
    RA_VaR = RA_VaR,
    RA_CTE = RA_CTE,
    LIC_VaR = mean_bel + RA_VaR,
    LIC_CTE = mean_bel + RA_CTE
  )
  
  list(
    full_data = tmp,
    bel_ay_draw = bel_ay_draw,
    bel_draw = bel_draw,
    summary = summary_tbl
  )
}


# =========================================================
# Predictive mean lower triangle
# =========================================================
predictive_mean_lower_triangle <- function(sim_obj){
  pred_mean <- colMeans(sim_obj$y_pred)
  
  sim_obj$pred_map %>%
    mutate(paid_cum_pred_mean = pred_mean)
}

# =========================================================
# Full triangle mean (observed + predictive mean)
# =========================================================
full_triangle_mean <- function(sim_obj){
  obs_long <- sim_obj$obs_df %>%
    select(ay_index, dev_year, paid_cum)
  
  pred_cell_mean <- predictive_mean_lower_triangle(sim_obj) %>%
    select(ay_index, dev_year, paid_cum_pred_mean)
  
  full_mean <- full_join(
    obs_long,
    pred_cell_mean,
    by = c("ay_index", "dev_year")
  ) %>%
    mutate(
      paid_cum_mean = coalesce(paid_cum, paid_cum_pred_mean)
    ) %>%
    arrange(ay_index, dev_year)
  
  # enforce monotone on mean triangle as well
  full_mean <- full_mean %>%
    group_by(ay_index) %>%
    mutate(paid_cum_mean = cummax(paid_cum_mean)) %>%
    ungroup()
  
  full_mean
}


# misal final fitted object lo = m3_se
# 1 = Properti, 2 = Credit, 3 = Marine

sim_prop <- simulate_lower_triangle_m3_se(
  res_obj = m3_se_nc,
  lob_obj = lob_list$prop,
  lob_index = 1,
  M = 1000,
  seed = 123
)

sim_cred <- simulate_lower_triangle_m3_se(
  res_obj = m3_se_nc,
  lob_obj = lob_list$cred,
  lob_index = 2,
  M = 1000,
  seed = 123
)

sim_marine <- simulate_lower_triangle_m3_se(
  res_obj = m3_se_nc,
  lob_obj = lob_list$marine,
  lob_index = 3,
  M = 1000,
  seed = 123
)

tri_prop   <- build_full_triangle_draws(sim_prop,   lob_name = "Properti")
tri_cred   <- build_full_triangle_draws(sim_cred,   lob_name = "Credit")
tri_marine <- build_full_triangle_draws(sim_marine, lob_name = "Marine")

tri_prop
tri_cred
tri_marine

make_mean_cum_triangle <- function(tri_obj){
  
  tri_obj %>%
    group_by(ay_index, dev_year) %>%
    summarise(
      mean_paid_cum = mean(paid_cum, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = dev_year,
      values_from = mean_paid_cum
    ) %>%
    arrange(ay_index)
  
}
make_mean_incr_triangle <- function(tri_obj){
  
  tri_obj %>%
    group_by(ay_index, dev_year) %>%
    summarise(
      mean_incr = mean(incr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = dev_year,
      values_from = mean_incr
    ) %>%
    arrange(ay_index)
  
}

cum_tri_prop <- make_mean_cum_triangle(tri_prop)
incr_tri_prop <- make_mean_incr_triangle(tri_prop)
cum_tri_cred <- make_mean_cum_triangle(tri_cred)
incr_tri_cred <- make_mean_incr_triangle(tri_cred)
cum_tri_marine <- make_mean_cum_triangle(tri_marine)
incr_tri_marine <- make_mean_incr_triangle(tri_marine)

cum_tri_prop
incr_tri_prop

cum_tri_cred
incr_tri_cred

cum_tri_marine
incr_tri_marine
#============Misal
disc_curve <- c(0.0252, 0.0252, 0.0252, 0.0252, 0.0252, 0.0252,
                0.0252, 0.0252, 0.0252, 0.0252, 0.0252)

#=====================

res_prop <- compute_bel_ra_lic(
  full_triangle_draws = tri_prop,
  conf = 0.75,
  discount_rates = disc_curve
)

res_cred <- compute_bel_ra_lic(
  full_triangle_draws = tri_cred,
  conf = 0.75,
  discount_rates = disc_curve
)

res_marine <- compute_bel_ra_lic(
  full_triangle_draws = tri_marine,
  conf = 0.75,
  discount_rates = disc_curve
)

res_prop$summary
res_cred$summary
res_marine$summary

#===========================
#88th Percentile
#===========================


res_prop88 <- compute_bel_ra_lic(
  full_triangle_draws = tri_prop,
  conf = 0.88,
  discount_rates = disc_curve
)

res_cred88 <- compute_bel_ra_lic(
  full_triangle_draws = tri_cred,
  conf = 0.88,
  discount_rates = disc_curve
)

res_marine88 <- compute_bel_ra_lic(
  full_triangle_draws = tri_marine,
  conf = 0.88,
  discount_rates = disc_curve
)

res_prop88$summary
res_cred88$summary
res_marine88$summary



#Plot distribusi BEL
# =========================
# Properti
# =========================

hist(res_prop$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Properti",
     xlab = "BEL",
     col = "gray",
     border = "black")

# Mean BEL
abline(v = res_prop$summary$BEL_mean,
       col = "blue",
       lwd = 2)

# VaR
abline(v = res_prop$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

# CTE
abline(v = res_prop$summary$BEL_mean + res_prop$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL", "VaR 75%", "CTE 75%"),
       col = c("blue", "red", "darkgreen"),
       lwd = 2,
       lty = c(1,2,3),
       bty = "n")


# =========================
# Kredit
# =========================

hist(res_cred$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Kredit",
     xlab = "BEL",
     col = "gray",
     border = "black")

abline(v = res_cred$summary$BEL_mean,
       col = "blue",
       lwd = 2)

abline(v = res_cred$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

abline(v = res_cred$summary$BEL_mean + res_cred$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL", "VaR 75%", "CTE 75%"),
       col = c("blue", "red", "darkgreen"),
       lwd = 2,
       lty = c(1,2,3),
       bty = "n")


# =========================
# Maritim
# =========================

hist(res_marine$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Maritim",
     xlab = "BEL",
     col = "gray",
     border = "black")

abline(v = res_marine$summary$BEL_mean,
       col = "blue",
       lwd = 2)

abline(v = res_marine$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

abline(v = res_marine$summary$BEL_mean + res_marine$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL", "VaR 75%", "CTE 75%"),
       col = c("blue", "red", "darkgreen"),
       lwd = 2,
       lty = c(1,2,3),
       bty = "n")

#=========
# 88 75 
#Properti
hist(res_prop$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Properti",
     xlab = "BEL",
     col = "gray",
     border = "black")

# Mean BEL
abline(v = res_prop$summary$BEL_mean,
       col = "blue",
       lwd = 2)

# VaR 75%
abline(v = res_prop$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

# CTE 75%
abline(v = res_prop$summary$BEL_mean + res_prop$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

# VaR 88%
abline(v = res_prop88$summary$BEL_q_conf,
       col = "purple",
       lwd = 2,
       lty = 2)

# CTE 88%
abline(v = res_prop88$summary$BEL_mean + res_prop88$summary$RA_CTE,
       col = "orange",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL",
                  "VaR 75%",
                  "CTE 75%",
                  "VaR 88%",
                  "CTE 88%"),
       col = c("blue",
               "red",
               "darkgreen",
               "purple",
               "orange"),
       lwd = 2,
       lty = c(1,2,3,2,3),
       bty = "n")

#Kredit
hist(res_cred$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Kredit",
     xlab = "BEL",
     col = "gray",
     border = "black")

# Mean BEL
abline(v = res_cred$summary$BEL_mean,
       col = "blue",
       lwd = 2)

# VaR 75%
abline(v = res_cred$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

# CTE 75%
abline(v = res_cred$summary$BEL_mean + res_cred$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

# VaR 88%
abline(v = res_cred88$summary$BEL_q_conf,
       col = "purple",
       lwd = 2,
       lty = 2)

# CTE 88%
abline(v = res_cred88$summary$BEL_mean + res_cred88$summary$RA_CTE,
       col = "orange",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL",
                  "VaR 75%",
                  "CTE 75%",
                  "VaR 88%",
                  "CTE 88%"),
       col = c("blue",
               "red",
               "darkgreen",
               "purple",
               "orange"),
       lwd = 2,
       lty = c(1,2,3,2,3),
       bty = "n")


#Marine
hist(res_marine$bel_draw$BEL,
     breaks = 40,
     main = "Posterior Distribution of BEL - Maritim",
     xlab = "BEL",
     col = "gray",
     border = "black")

# Mean BEL
abline(v = res_marine$summary$BEL_mean,
       col = "blue",
       lwd = 2)

# VaR 75%
abline(v = res_marine$summary$BEL_q_conf,
       col = "red",
       lwd = 2,
       lty = 2)

# CTE 75%
abline(v = res_marine$summary$BEL_mean + res_marine$summary$RA_CTE,
       col = "darkgreen",
       lwd = 2,
       lty = 3)

# VaR 88%
abline(v = res_marine88$summary$BEL_q_conf,
       col = "purple",
       lwd = 2,
       lty = 2)

# CTE 88%
abline(v = res_marine88$summary$BEL_mean + res_marine88$summary$RA_CTE,
       col = "orange",
       lwd = 2,
       lty = 3)

legend("topright",
       legend = c("Mean BEL",
                  "VaR 75%",
                  "CTE 75%",
                  "VaR 88%",
                  "CTE 88%"),
       col = c("blue",
               "red",
               "darkgreen",
               "purple",
               "orange"),
       lwd = 2,
       lty = c(1,2,3,2,3),
       bty = "n")
