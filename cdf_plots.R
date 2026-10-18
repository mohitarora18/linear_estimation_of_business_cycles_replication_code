rm(list=ls())
options(scipen = 999)

# =====================================================================
# --- PLOT SETTINGS (shared by all CDF figures) ---
# Figures are saved at 6.5 x 7.8 inches so text stays legible when each
# panel is shrunk to one third of the page width in LaTeX.
# =====================================================================
CEX_AXIS  <- 1.7   # tick labels
CEX_LAB   <- 1.7   # axis titles
CEX_MARK  <- 1.6   # value printed at x = 1
open_png <- function(file) {
  png(file, width = 6.5, height = 7.8, units = "in", res = 400)
  par(mar = c(8.5, 5, 1, 1.5), mgp = c(3, 1.1, 0))
}
# Multi-line x-axis title: "<name> Ratio" on top, "(nonlinear/linear)" below.
# Long names with a bracketed part, e.g. "Downturn Severity (Max Crash)",
# are split over two lines so the title fits the panel width.
xlab2 <- function(name) {
  if (grepl("\\(", name) && nchar(name) > 22) name <- sub(" \\(", "\n(", name)
  paste0(name, " Ratio\n(nonlinear/linear)")
}
# Write the x-axis title one line at a time below the tick labels, so
# two- and three-line titles start at the same height and never overlap
# the tick labels
title_x <- function(lab) {
  lines <- strsplit(lab, "\n")[[1]]
  for (i in seq_along(lines))
    mtext(lines[i], side = 1, line = 3.3 + (i - 1) * 1.5, cex = CEX_LAB)
}


# =====================================================================
# --- 0. BBQ-STYLE TURNING POINT / AMPLITUDE FUNCTIONS ---
# (Harding & Pagan, 2002, "Dissecting the Cycle: A Methodological
#  Investigation," Journal of Monetary Economics 49(2), 365-381,
#  building on Bry & Boschan, 1971, NBER Technical Paper No. 20.)
#
# Parameters chosen for these simulations (300-period, non-seasonal
# AR(3)/VAR(1) draws -- NOT quarterly macro data, so the literature's
# quarterly defaults k=2/min_phase=2/min_cycle=5 are relaxed):
#   k         = 1   (local-extremum window)
#   min_phase = 1   (minimum periods between consecutive turning points)
#   min_cycle = 2   (minimum periods between two turning points of the
#                    same type, i.e. peak-to-peak or trough-to-trough)
#
# bbq_turning_points(): implements the three core BBQ steps --
#   (1) locate candidate local maxima/minima in a +/-k window,
#   (2) enforce strict alternation between peaks and troughs,
#   (3) apply minimum-phase / minimum-cycle censoring rules.
#
# bbq_metrics(): given a series, returns
#   amplitude = mean |A_i| = mean absolute peak-to-trough AND
#               trough-to-peak movement across all dated phases
#               (both phases averaged)
# =====================================================================

bbq_turning_points <- function(y, k = 1, min_phase = 1, min_cycle = 2) {
  n <- length(y)
  if (n < (2 * k + 3)) {
    return(list(idx = integer(0), type = character(0)))
  }
  
  # --- Step 1: candidate local extrema in a +/-k window ---
  cand_type <- rep(NA_character_, n)
  for (t in (k + 1):(n - k)) {
    win <- y[(t - k):(t + k)]
    if (y[t] == max(win)) cand_type[t] <- "P"
    if (y[t] == min(win)) cand_type[t] <- "T"
  }
  idx  <- which(!is.na(cand_type))
  type <- cand_type[idx]
  if (length(idx) < 2) return(list(idx = idx, type = type))
  
  enforce_alternation <- function(idx, type) {
    repeat {
      changed <- FALSE
      i <- 1
      while (i < length(type)) {
        if (type[i] == type[i + 1]) {
          if (type[i] == "P") {
            drop <- if (y[idx[i]] >= y[idx[i + 1]]) i + 1 else i
          } else {
            drop <- if (y[idx[i]] <= y[idx[i + 1]]) i + 1 else i
          }
          idx  <- idx[-drop]
          type <- type[-drop]
          changed <- TRUE
          break
        }
        i <- i + 1
      }
      if (!changed) break
    }
    list(idx = idx, type = type)
  }
  
  # --- Step 2: enforce strict peak/trough alternation ---
  out  <- enforce_alternation(idx, type)
  idx  <- out$idx; type <- out$type
  
  # --- Step 3a: minimum phase length censoring ---
  repeat {
    if (length(idx) < 2) break
    d <- diff(idx)
    bad <- which(d < min_phase)
    if (length(bad) == 0) break
    i <- bad[1]
    idx  <- idx[-(i + 1)]
    type <- type[-(i + 1)]
    out  <- enforce_alternation(idx, type)
    idx  <- out$idx; type <- out$type
  }
  
  # --- Step 3b: minimum full-cycle length censoring (same-type spacing) ---
  repeat {
    if (length(idx) < 3) break
    cyc_d <- idx[3:length(idx)] - idx[1:(length(idx) - 2)]
    bad <- which(cyc_d < min_cycle)
    if (length(bad) == 0) break
    i <- bad[1]
    idx  <- idx[-(i + 1)]
    type <- type[-(i + 1)]
  }
  
  list(idx = idx, type = type)
}

bbq_metrics <- function(y, k = 1, min_phase = 1, min_cycle = 2) {
  tp <- bbq_turning_points(y, k = k, min_phase = min_phase, min_cycle = min_cycle)
  idx <- tp$idx
  
  if (length(idx) < 2) {
    return(c(amplitude = NA_real_))
  }
  
  amp <- abs(y[idx[-1]] - y[idx[-length(idx)]])  # |A_i| for every phase, both directions
  
  c(amplitude = mean(amp))
}

# =====================================================================
# --- 0b. AUTOCORRELATION FUNCTION ---
# Lag-k sample autocorrelation (standard estimator used by stats::acf).
# Returns NA if the series is constant (autocorrelation undefined).
# =====================================================================
# Lag reported in the table: 1

acf_k <- function(y, k) {
  y <- as.numeric(y)
  if (sd(y) == 0) return(NA_real_)
  as.numeric(acf(y, lag.max = k, plot = FALSE)$acf[k + 1])
}

# =====================================================================
# --- 1. Define Fixed Parameters (Beaudry 2017) ---
# =====================================================================
sigma_eps <- 0.25
p <- 3
Burn_In <- 600     # Number of periods to discard (transient phase)
Num_Obs <- 300     # Number of periods to actually plot/analyze
Total_Sim <- Burn_In + Num_Obs
# =====================================================================
# --- 2. MONTE CARLO SIMULATION ---
# =====================================================================
M <- 10000 # number of simulations
# Initialize storage vectors for the metrics
mc_sd_nl_beaudry <- numeric(M)
mc_sd_lin_beaudry <- numeric(M)
mc_eigen_lin_beaudry <- numeric(M)
true_eigen_dgp_beaudry <- numeric(M)
mc_min_nl_beaudry <- numeric(M)
mc_min_lin_beaudry <- numeric(M)
mc_Q1_nl_beaudry <- numeric(M);
mc_Q1_lin_beaudry <- numeric(M);
# peak-trough movement storage
mc_amp_nl_beaudry   <- numeric(M); mc_amp_lin_beaudry   <- numeric(M)

# autocorrelation storage
mc_ac1_nl_beaudry <- numeric(M); mc_ac1_lin_beaudry <- numeric(M)  # lag 1
# unconditional mean storage
mc_mean_nl_beaudry <- numeric(M); mc_mean_lin_beaudry <- numeric(M)

#CDF data storage: nonlinear to linear ratios
peak_trough_movement_ratio_beaudry <- numeric(M) 
volatility_ratio_beaudry <- numeric(M) 
downturn_first_Q_ratio_beaudry <- numeric(M) 
max_crash_ratio_beaudry <- numeric(M) 
ac1_ratio_beaudry <- numeric(M) 
max_eigenvalue_ratio_beaudry <- numeric(M) 
mean_ratio_beaudry <- numeric(M) 

set.seed(123) # Set master seed for reproducibility of the MC loop
for (m in 1:M) {
  # A. Generate model parameters   
  a1    <- runif(1, 0.60, 1.50)    #1.5
  a2    <- runif(1, -0.625, -0.575)  #-0.6
  a3    <- runif(1, -0.325, -0.275)  #-0.3
  a_cub <- runif(1, -0.010, -0.005) #-0.01
  
  Comp_Mat_DGP <- rbind(c(a1, a2, a3),
                        c(1, 0, 0),
                        c(0, 1, 0))
  true_eigen_dgp_beaudry[m] <- max(Mod(eigen(Comp_Mat_DGP)$values))
  
  # B. Generate fresh noise for this iteration
  eps_m <- rnorm(Total_Sim + p)
  x_m <- numeric(Total_Sim + p)
  
  # C. Simulate true Nonlinear DGP
  for(i in (p+1):(Total_Sim+p)){
    x_m[i] <- a1*x_m[i-1] + a2*x_m[i-2] + a3*x_m[i-3] + a_cub*(x_m[i-1]^3) + sigma_eps*eps_m[i]
  }
  
  # Isolate post-burn-in
  series_m <- x_m[(Burn_In + p + 1):(Total_Sim + p)]
  ts_m <- ts(series_m)
  
  # D. Estimate Linear AR(3) on the isolated series
  d_m <- ts.intersect(ts_m, ts.L1=stats::lag(ts_m,-1), ts.L2=stats::lag(ts_m,-2), ts.L3=stats::lag(ts_m,-3))
  lin_mod_m <- lm(ts_m ~ ts.L1 + ts.L2 + ts.L3, data=d_m)
  
  c_m  <- coef(lin_mod_m)[1]
  a1_m <- coef(lin_mod_m)[2]
  a2_m <- coef(lin_mod_m)[3]
  a3_m <- coef(lin_mod_m)[4]
  
  # E. Calculate Max Eigenvalue using a Companion Matrix
  Comp_Mat <- rbind(c(a1_m, a2_m, a3_m),
                    c(1, 0, 0),
                    c(0, 1, 0))
  mc_eigen_lin_beaudry[m] <- max(Mod(eigen(Comp_Mat)$values))
  
  # F. Simulate Linear Forward with concurrent noise
  x_lin_m <- numeric(Num_Obs)
  x_lin_m[1:3] <- series_m[1:3] # Start from same condition
  for (t in 4:Num_Obs) {
    x_lin_m[t] <- c_m + a1_m*x_lin_m[t-1] + a2_m*x_lin_m[t-2] + a3_m*x_lin_m[t-3] + sigma_eps*eps_m[Burn_In + p + t]
  }
  
  # G. Calculate & Store the metrics
  # 1. Volatility (Standard Deviation)
  mc_sd_nl_beaudry[m] <- sd(series_m)
  mc_sd_lin_beaudry[m] <- sd(x_lin_m)
  
  # 2. Maximum Crash (Minimum Value)
  mc_min_nl_beaudry[m] <- min(series_m)
  mc_min_lin_beaudry[m] <- min(x_lin_m)
  
  
  # 3. First quatrile
  mc_Q1_nl_beaudry[m] <- as.numeric(quantile(series_m,0.25))
  mc_Q1_lin_beaudry[m] <- as.numeric(quantile(x_lin_m,0.25))
  
  # 4. Peak-trough movement (Harding & Pagan, 2002)
  bbq_nl  <- bbq_metrics(series_m, k = 1, min_phase = 1, min_cycle = 2)
  bbq_lin <- bbq_metrics(x_lin_m,  k = 1, min_phase = 1, min_cycle = 2)
  mc_amp_nl_beaudry[m]    <- bbq_nl["amplitude"]
  mc_amp_lin_beaudry[m]   <- bbq_lin["amplitude"]
  
  
  # 5. Autocorrelation (lag 1)
  mc_ac1_nl_beaudry[m]  <- acf_k(series_m, 1)
  mc_ac1_lin_beaudry[m] <- acf_k(x_lin_m, 1)
  
  # 6. Unconditional mean
  mc_mean_nl_beaudry[m]  <- mean(series_m)
  mc_mean_lin_beaudry[m] <- mean(x_lin_m)
  
  
  #CDF data:
  peak_trough_movement_ratio_beaudry[m] <-  mc_amp_nl_beaudry[m] / mc_amp_lin_beaudry[m]
  volatility_ratio_beaudry[m] <-  mc_sd_nl_beaudry[m] / mc_sd_lin_beaudry[m]
  downturn_first_Q_ratio_beaudry[m] <-    mc_Q1_nl_beaudry[m]/mc_Q1_lin_beaudry[m]
  max_crash_ratio_beaudry[m] <-  mc_min_nl_beaudry[m] /   mc_min_lin_beaudry[m]
  ac1_ratio_beaudry[m] <- mc_ac1_nl_beaudry[m]/mc_ac1_lin_beaudry[m] 
  max_eigenvalue_ratio_beaudry[m] <-  true_eigen_dgp_beaudry[m] /  mc_eigen_lin_beaudry[m] 
  mean_ratio_beaudry[m] <- mc_mean_nl_beaudry[m] / mc_mean_lin_beaudry[m]
  
  
}


# =====================================================================
# --- 4. PLOT CDFs (Beaudry) ---
# =====================================================================

# Automatic x-range: middle 99% of runs (0.5th to 99.5th percentile),
# always including 1, with a small margin on each side
auto_xlim <- function(x) {
  r <- quantile(x, c(0.005, 0.995), na.rm = TRUE)
  r <- range(r, 1)
  r + c(-1, 1) * 0.05 * diff(r)
}

# Draw the x-axis ourselves: R's usual tick positions plus a tick at 1 that
# is always labelled. Ticks whose labels would overlap the "1" label are
# dropped, and the "1" tick is drawn in its own call so R never omits it.
# Labels are written as plain numbers (no scientific notation).
draw_x_axis <- function(log = FALSE) {
  fmt <- function(v) sapply(v, function(z) format(z, scientific = FALSE, trim = TRUE,
                                                  drop0trailing = TRUE, big.mark = ","))
  ticks  <- axTicks(1)
  ticks  <- ticks[abs(ticks - 1) > 1e-9]
  pos_in <- function(v) grconvertX(v, "user", "inches")
  w_in   <- function(lab) strwidth(lab, units = "inches", cex = CEX_AXIS)
  keep   <- abs(pos_in(ticks) - pos_in(1)) > 0.6 * (w_in(fmt(ticks)) + w_in("1")) + 0.05
  ticks  <- ticks[keep]
  if (length(ticks)) axis(1, at = ticks, labels = fmt(ticks), cex.axis = CEX_AXIS)
  axis(1, at = 1, labels = "1", cex.axis = CEX_AXIS)
}

# Mark the CDF at x = 1: vertical segment from the x-axis up to the curve,
# horizontal segment from the curve to the y-axis, and the y-value printed
mark_at_one <- function(x) {
  F1 <- ecdf(x)(1)
  segments(x0 = 1, y0 = par("usr")[3], x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  segments(x0 = par("usr")[1], y0 = F1, x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  text(x = par("usr")[1], y = F1, labels = sprintf("%.2f", F1),
       adj = c(-0.15, -0.5), col = "blue", cex = CEX_MARK)
}


# peak-trough movement ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_peak_trough_movement_ratio_CDF.png")
plot(ecdf(peak_trough_movement_ratio_beaudry), xlim=auto_xlim(peak_trough_movement_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Peak-Trough Movement"))
mark_at_one(peak_trough_movement_ratio_beaudry)
dev.off()  

# volatility ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_volatility_ratio_CDF.png")
plot(ecdf(volatility_ratio_beaudry), xlim=auto_xlim(volatility_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Volatility"))
mark_at_one(volatility_ratio_beaudry)
dev.off()  

# downturn first  quartile ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_downturn_first_Q_ratio_CDF.png")
plot(ecdf(downturn_first_Q_ratio_beaudry), xlim=auto_xlim(downturn_first_Q_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Downturn Severity (First Quartile)"))
mark_at_one(downturn_first_Q_ratio_beaudry)
dev.off()    

# max crash ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_max_crash_ratio_CDF.png")
plot(ecdf(max_crash_ratio_beaudry), xlim=auto_xlim(max_crash_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Downturn Severity (Max Crash)"))
mark_at_one(max_crash_ratio_beaudry)
dev.off() 

#ac1 ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_ac1_ratio_CDF.png")
plot(ecdf(ac1_ratio_beaudry), xlim=auto_xlim(ac1_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Autocorrelation (lag 1)"))
mark_at_one(ac1_ratio_beaudry)
dev.off() 

# max eig val ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_max_eig_val_ratio_CDF.png")
plot(ecdf(max_eigenvalue_ratio_beaudry), xlim=auto_xlim(max_eigenvalue_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Max Eigenvalue (Modulus)"))
mark_at_one(max_eigenvalue_ratio_beaudry)
dev.off()

# unconditional mean ratio CDF 
open_png("~/Desktop/Applied_EL_Revisions/Beaudry/Beaudry_mean_ratio_CDF.png")
plot(ecdf(mean_ratio_beaudry), xlim=auto_xlim(mean_ratio_beaudry), xaxt = "n", xlab = "", ylab = "", main = "", lwd=2.5, cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
draw_x_axis()
title_x(xlab2("Unconditional Mean"))
mark_at_one(mean_ratio_beaudry)
dev.off()




#SM39
library(mvtnorm)

# =====================================================================
# --- 1. Define Fixed Parameters (SM39 Bivariate System) ---
# =====================================================================

p <- 1 # 1 lag in the system

Burn_In <- 600     # Discard transient phase
Num_Obs <- 300     # Periods to plot/analyze
Total_Sim <- Burn_In + Num_Obs
# Multivariate Shock Covariance Matrix
eps_var <- matrix(c(0.0001, 0.00,
                    0.00, 0.0001),
                  nrow = 2, byrow = TRUE)
# =====================================================================
# --- 2. MONTE CARLO SIMULATION ---
# =====================================================================
M <- 10000 # number of simulations
# Initialize system storage
mc_eigen_lin_2D <- numeric(M)
true_eigen_dgp_2D <- numeric(M)
# Initialize variables storage (only Output 'y' is retained for metrics)
mc_sd_nl_y <- numeric(M);
mc_sd_lin_y <- numeric(M);
mc_min_nl_y <- numeric(M);
mc_min_lin_y <- numeric(M);
mc_Q1_nl_y <- numeric(M);
mc_Q1_lin_y <- numeric(M);
# peak-trough movement storage
mc_amp_nl_y   <- numeric(M); mc_amp_lin_y   <- numeric(M)
# autocorrelation storage
mc_ac1_nl_y <- numeric(M); mc_ac1_lin_y <- numeric(M)  # lag 1
# unconditional mean storage
mc_mean_nl_y <- numeric(M); mc_mean_lin_y <- numeric(M)

#CDF data storage: nonlinear to linear ratios
peak_trough_movement_ratio_y <- numeric(M) #nonlinear to linear
volatility_ratio_y <- numeric(M) 
downturn_first_Q_ratio_y <- numeric(M) 
max_crash_ratio_y <- numeric(M) 
ac1_ratio_y <- numeric(M) 
max_eigenvalue_ratio_y <- numeric(M) 
mean_ratio_y <- numeric(M) 

set.seed(132)
for (m in 1:M) {
  # A. Generate model parameters
  y_max <- runif(1,0,150)  #100
  I_min <- runif(1,-75,0) #-50
  c <- runif(1,0,1)  #0.9
  beta <- runif(1,0,5) #2
  
  B_true <- matrix(c((c * beta) + c, -beta,
                     c,              0),
                   nrow=2, ncol=2, byrow=TRUE)
  true_eigen_dgp_2D[m] <- max(Mod(eigen(B_true)$values))
  
  
  
  # B. Generate multivariate noise
  shocks_m <- t(rmvnorm(Total_Sim + p, mean = c(0, 0), sigma = eps_var))
  
  ym <- numeric(Total_Sim + p)
  cm <- numeric(Total_Sim + p)
  ym[1] <- 0; cm[1] <- 0
  
  # C. Simulate DGP
  for (t in 2:(Total_Sim + p)) {
    
    cm[t] <- c*ym[t-1] + shocks_m[2, t]
    I_tilde <- beta*(cm[t]-cm[t-1])
    ym[t] <- min(max(I_tilde,I_min) + cm[t] + shocks_m[1, t], y_max)
    
  }
  
  # D. Extract & Estimate
  ts_y <- ts(ym[(Burn_In + p + 1):(Total_Sim + p)])
  ts_c <- ts(cm[(Burn_In + p + 1):(Total_Sim + p)])
  
  d_m <- ts.intersect(ts_y, ts_c, L1.y=stats::lag(ts_y,-1), L1.c=stats::lag(ts_c,-1))
  r1_m <- lm(ts_y ~ L1.y + L1.c, data = d_m)
  r2_m <- lm(ts_c ~ L1.y + L1.c, data = d_m)
  
  # E. VAR(1) Jacobian & Eigenvalue
  B_mat <- matrix(c(coef(r1_m)["L1.y"], coef(r1_m)["L1.c"],
                    coef(r2_m)["L1.y"], coef(r2_m)["L1.c"]),
                  nrow=2, ncol=2, byrow=TRUE)
  mc_eigen_lin_2D[m] <- max(Mod(eigen(B_mat)$values))
  
  # F. Simulate Linear Forward
  ylin_m <- numeric(Num_Obs); clin_m <- numeric(Num_Obs)
  ylin_m[1] <- ts_y[1]; clin_m[1] <- ts_c[1]
  
  for (t in 2:Num_Obs) {
    ylin_m[t] <- coef(r1_m)[1] + coef(r1_m)[2]*ylin_m[t-1] + coef(r1_m)[3]*clin_m[t-1] + shocks_m[1, Burn_In + p + t]
    clin_m[t] <- coef(r2_m)[1] + coef(r2_m)[2]*ylin_m[t-1] + coef(r2_m)[3]*clin_m[t-1] + shocks_m[2, Burn_In + p + t]
  }
  
  # G. Metrics Calculation (Output only)
  mc_sd_nl_y[m] <- sd(ts_y)
  mc_sd_lin_y[m] <- sd(ylin_m)
  mc_min_nl_y[m] <- min(ts_y)
  mc_min_lin_y[m] <- min(ylin_m)
  mc_Q1_nl_y[m] <- as.numeric(quantile(ts_y,0.25))
  mc_Q1_lin_y[m] <- as.numeric(quantile(ylin_m,0.25))
  
  # H. Peak-trough movement (Harding & Pagan, 2002)
  bbq_nl  <- bbq_metrics(as.numeric(ts_y), k = 1, min_phase = 1, min_cycle = 2)
  bbq_lin <- bbq_metrics(ylin_m,           k = 1, min_phase = 1, min_cycle = 2)
  mc_amp_nl_y[m]    <- bbq_nl["amplitude"]
  mc_amp_lin_y[m]   <- bbq_lin["amplitude"]
  
  # I. Autocorrelation (lag 1)
  mc_ac1_nl_y[m]  <- acf_k(ts_y, 1)
  mc_ac1_lin_y[m] <- acf_k(ylin_m, 1)
  
  # J. Unconditional mean
  mc_mean_nl_y[m]  <- mean(ts_y)
  mc_mean_lin_y[m] <- mean(ylin_m)
  
  
  
  #CDF data:
  peak_trough_movement_ratio_y[m] <-  mc_amp_nl_y[m] / mc_amp_lin_y[m]
  volatility_ratio_y[m] <-  mc_sd_nl_y[m] / mc_sd_lin_y[m]
  downturn_first_Q_ratio_y[m] <-    mc_Q1_nl_y[m]/mc_Q1_lin_y[m]
  max_crash_ratio_y[m] <-  mc_min_nl_y[m] /   mc_min_lin_y[m]
  ac1_ratio_y[m] <- mc_ac1_nl_y[m]/mc_ac1_lin_y[m] 
  max_eigenvalue_ratio_y[m] <-  true_eigen_dgp_2D[m] /  mc_eigen_lin_2D[m] 
  mean_ratio_y[m] <- mc_mean_nl_y[m] / mc_mean_lin_y[m]
  
}


# =====================================================================
# --- 4. PLOT CDFs (2D System) ---
# =====================================================================

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

# x-range covering the middle 99% of runs (0.5th to 99.5th percentile) by
# default, or the middle 95% (2.5th to 97.5th percentile) if coverage = 0.95,
# always including 1, with a margin of 5% of the range's width on each side
# With log = TRUE, the range and margin are computed
# on the log scale.
cdf_xlim <- function(x, log = FALSE, coverage = 0.99) {
  x <- x[is.finite(x)]
  if (log) x <- x[x > 0]
  r <- quantile(x, c((1 - coverage) / 2, 1 - (1 - coverage) / 2))
  r <- range(r, 1)
  if (log) {
    lr <- log10(r)
    10^(lr + c(-1, 1) * 0.05 * diff(lr))
  } else {
    r + c(-1, 1) * 0.05 * diff(r)
  }
}

# Decide whether a variable needs a log axis. A log axis is used only if
#  (a) all values are positive (a log axis cannot show zero or negatives), and
#  (b) on a regular axis the middle 80% of runs (10th-90th percentile) would
#      fill less than 35% of the axis width, i.e. long tails squash the plot.
needs_log <- function(x, min_share = 0.35) {
  x <- x[is.finite(x)]
  if (any(x <= 0)) return(FALSE)
  r <- range(quantile(x, c(0.005, 0.995)), 1)
  diff(quantile(x, c(0.10, 0.90))) / diff(r) < min_share
}

# Mark the CDF at x = 1: vertical line from the x-axis up to the curve,
# horizontal line from the curve to the y-axis, and the CDF value printed
mark_at_one <- function(x) {
  F1  <- ecdf(x)(1)
  usr <- par("usr")
  x_left <- if (par("xlog")) 10^usr[1] else usr[1]
  segments(x0 = 1, y0 = usr[3], x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  segments(x0 = x_left, y0 = F1, x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  text(x = x_left, y = F1, labels = sprintf("%.2f", F1),
       adj = c(-0.15, -0.5), col = "blue", cex = CEX_MARK)
}

# Draw the x-axis ourselves: R's usual tick positions plus a tick at 1 that
# is always labelled. Ticks whose labels would overlap the "1" label are
# dropped, and the "1" tick is drawn in its own call so R never omits it.
# Labels are written as plain numbers (no scientific notation).
draw_x_axis <- function(log = FALSE) {
  fmt <- function(v) sapply(v, function(z) format(z, scientific = FALSE, trim = TRUE,
                                                  drop0trailing = TRUE, big.mark = ","))
  ticks  <- axTicks(1)
  ticks  <- ticks[abs(ticks - 1) > 1e-9]
  pos_in <- function(v) grconvertX(v, "user", "inches")
  w_in   <- function(lab) strwidth(lab, units = "inches", cex = CEX_AXIS)
  keep   <- abs(pos_in(ticks) - pos_in(1)) > 0.6 * (w_in(fmt(ticks)) + w_in("1")) + 0.05
  ticks  <- ticks[keep]
  if (length(ticks)) axis(1, at = ticks, labels = fmt(ticks), cex.axis = CEX_AXIS)
  axis(1, at = 1, labels = "1", cex.axis = CEX_AXIS)
}

# Draw one CDF plot; the log axis is chosen automatically by needs_log().
# coverage sets the share of runs shown on a regular axis (default 99%)
plot_cdf <- function(x, file, xlab, coverage = 0.99) {
  lg <- needs_log(x)
  cat(sprintf("%-45s %s\n", file, if (lg) "log x-axis" else "regular x-axis"))
  open_png(file)
  plot(ecdf(x), xlim = cdf_xlim(x, log = lg, coverage = coverage), log = if (lg) "x" else "", xaxt = "n",
       xlab = "", ylab = "", main = "", lwd = 2.5,
       cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
  draw_x_axis(log = lg)
  title_x(xlab)
  mark_at_one(x)
  dev.off()
}

# ---------------------------------------------------------------------
# The seven CDF plots
# ---------------------------------------------------------------------

# peak-trough movement ratio CDF
plot_cdf(peak_trough_movement_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_peak_trough_movement_ratio_y_CDF.png",
         xlab2("Peak-Trough Movement"))

# volatility ratio CDF
plot_cdf(volatility_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_volatility_ratio_CDF.png",
         xlab2("Volatility"))

# downturn first quartile ratio CDF
plot_cdf(downturn_first_Q_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_downturn_first_Q_ratio_CDF.png",
         xlab2("Downturn Severity (First Quartile)"), coverage = 0.95)

# max crash ratio CDF
plot_cdf(max_crash_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_max_crash_ratio_CDF.png",
         xlab2("Downturn Severity (Max Crash)"))

# ac1 ratio CDF
plot_cdf(ac1_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_ac1_ratio_CDF.png",
         xlab2("Autocorrelation (lag 1)"), coverage = 0.95)

# max eig val ratio CDF
plot_cdf(max_eigenvalue_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_max_eig_val_ratio_CDF.png",
         xlab2("Max Eigenvalue (Modulus)"))

# unconditional mean ratio CDF
plot_cdf(mean_ratio_y, "~/Desktop/Applied_EL_Revisions/SM39/SM39_mean_ratio_CDF.png",
         xlab2("Unconditional Mean"))




#Three-dimensional BGST model
library(mvtnorm)
# =====================================================================
# --- 1. Fixed Define Parameters (BGST 3D System) ---
# =====================================================================

# Shock Covariance Matrix
eps_var <- matrix(c(0.0001, 0.00, 0.00,
                    0.00, 0.0001, 0.00,
                    0.00, 0.00, 0.0001),
                  nrow = 3, byrow = TRUE)
p <- 1             # Lags in the estimated VAR(1)
Burn_In <- 600     # Discard transient phase
Num_Obs <- 300     # Periods to analyze/plot
Total_Sim <- Burn_In + Num_Obs

# =====================================================================
# --- 3. MONTE CARLO SIMULATION ---
# =====================================================================
M <- 10000 # number of simulations
# Storage
mc_eigen_lin_3D <- numeric(M)
true_eigen_dgp_3D <- numeric(M)
mc_sd_nl_inv <- numeric(M); mc_sd_lin_inv <- numeric(M)
mc_min_nl_inv <- numeric(M); mc_min_lin_inv <- numeric(M)
mc_Q1_nl_inv <- numeric(M); mc_Q1_lin_inv <- numeric(M)
# peak-trough movement storage
mc_amp_nl_inv   <- numeric(M); mc_amp_lin_inv   <- numeric(M)
# autocorrelation storage
mc_ac1_nl_inv <- numeric(M); mc_ac1_lin_inv <- numeric(M)  # lag 1
# unconditional mean storage
mc_mean_nl_inv <- numeric(M); mc_mean_lin_inv <- numeric(M)
#CDF data storage: nonlinear to linear ratios
peak_trough_movement_ratio_BGST <- numeric(M) #nonlinear to linear
volatility_ratio_BGST <- numeric(M) 
downturn_first_Q_ratio_BGST <- numeric(M) 
max_crash_ratio_BGST <- numeric(M) 
ac1_ratio_BGST <- numeric(M) 
max_eigenvalue_ratio_BGST <- numeric(M) 
mean_ratio_BGST <- numeric(M) 
set.seed(456)
for (m in 1:M) {
  
  ###  
  rho      <- runif(1,0.90,0.98) #0.95
  theta    <- runif(1,0.5,1.5) #0.913
  eta      <- runif(1,5,10) #8
  phi_inv  <- runif(1,0.1,0.5) #0.25
  phi_cost <- runif(1,0.5,1) #0.67
  repay    <- runif(1,0.05,0.2) #0.1
  gamma <- runif(1, 0.1, 0.6)  # 0.5
  alpha <- runif(1, 1.25, 5)  # 2.5
  
  # =====================================================================
  # --- 2. ANALYTICAL PROOF: True DGP Stability ---
  # =====================================================================
  ss_equation <- function(b) { (repay * b) + (gamma * alpha * phi_cost * exp(eta * b)) }
  b_star <- uniroot(ss_equation, lower = -2, upper = 0)$root
  E_star <- alpha * exp(eta * b_star)
  I_star <- -phi_cost * E_star
  J_true <- matrix(0, nrow=3, ncol=3)
  rownames(J_true) <- c("TFP_t", "Debt_t", "Inv_t")
  colnames(J_true) <- c("TFP_t-1", "Debt_t-1", "Inv_t-1")
  J_true[1, 1] <- rho
  J_true[3, 1] <- phi_inv*(rho^2) + (eta * rho * phi_cost * E_star)
  J_true[3, 2] <- -eta * phi_cost * E_star
  J_true[2, 1] <- gamma * J_true[3, 1]
  J_true[2, 2] <- (1 - repay) + (gamma * J_true[3, 2])
  true_eigen_dgp_3D[m] <- max(Mod(eigen(J_true)$values))
  
  # =====================================================================
  # --- 3. Simulate Dynamics ---
  # =====================================================================
  
  shocks_m <- t(rmvnorm(Total_Sim + p, mean=c(0,0,0), sigma = eps_var))
  xm_NL <- matrix(0, nrow = 3, ncol = Total_Sim + p)
  
  for(t in 2:(Total_Sim + p)){
    news <- shocks_m[1,t]
    xm_NL[1,t] <- rho * xm_NL[1,t-1] + news
    belief <- rho * xm_NL[1,t] + theta * news
    spread_NL <- alpha*(exp(eta * (xm_NL[2,t-1] - xm_NL[1,t])))
    xm_NL[3,t] <- phi_inv * belief - phi_cost * spread_NL + shocks_m[3,t]
    xm_NL[2,t] <- (1-repay) * xm_NL[2,t-1] + gamma * xm_NL[3,t] + shocks_m[2,t]
  }
  
  ts_full_m <- ts(t(xm_NL))
  colnames(ts_full_m) <- c("TFP", "Debt", "Invest")
  series_m <- window(ts_full_m, start = Burn_In + p + 1)
  
  dm1 <- ts.intersect(series_m, series_m.L1=stats::lag(series_m,-1))
  
  rm_TFP <- lm(series_m.TFP ~ series_m.L1.TFP + series_m.L1.Debt + series_m.L1.Invest, data = dm1)
  rm_Debt <- lm(series_m.Debt ~ series_m.L1.TFP + series_m.L1.Debt + series_m.L1.Invest, data = dm1)
  rm_Inv <- lm(series_m.Invest ~ series_m.L1.TFP + series_m.L1.Debt + series_m.L1.Invest, data = dm1)
  
  J_m <- rbind(coef(rm_TFP)[2:4], coef(rm_Debt)[2:4], coef(rm_Inv)[2:4])
  mc_eigen_lin_3D[m] <- max(Mod(eigen(J_m)$values))
  
  cm_real <- c(coef(rm_TFP)[1], coef(rm_Debt)[1], coef(rm_Inv)[1])
  States_m <- matrix(0, nrow = 3, ncol = Num_Obs)
  States_m[, 1] <- series_m[1, ]
  
  for (t in 2:Num_Obs) {
    States_m[,t] <- cm_real + J_m %*% States_m[,t-1]
    States_m[, t] <- States_m[, t] + shocks_m[, Burn_In + p + t]
  }
  
  mc_sd_nl_inv[m] <- sd(series_m[, "Invest"])
  mc_sd_lin_inv[m] <- sd(States_m[3, ])
  
  mc_min_nl_inv[m] <- min(series_m[, "Invest"])
  mc_min_lin_inv[m] <- min(States_m[3, ])
  
  mc_Q1_nl_inv[m] <- as.numeric(quantile(series_m[, "Invest"],0.25))
  mc_Q1_lin_inv[m] <- as.numeric(quantile(States_m[3, ],0.25))
  
  # Peak-trough movement (Harding & Pagan, 2002)
  bbq_nl  <- bbq_metrics(as.numeric(series_m[, "Invest"]), k = 1, min_phase = 1, min_cycle = 2)
  bbq_lin <- bbq_metrics(States_m[3, ],                     k = 1, min_phase = 1, min_cycle = 2)
  mc_amp_nl_inv[m]    <- bbq_nl["amplitude"]
  mc_amp_lin_inv[m]   <- bbq_lin["amplitude"]
  
  # Autocorrelation (lag 1)
  mc_ac1_nl_inv[m]  <- acf_k(series_m[, "Invest"], 1)
  mc_ac1_lin_inv[m] <- acf_k(States_m[3, ], 1)
  
  # Unconditional mean
  mc_mean_nl_inv[m]  <- mean(series_m[, "Invest"])
  mc_mean_lin_inv[m] <- mean(States_m[3, ])
  
  
  #CDF data:
  peak_trough_movement_ratio_BGST[m] <-  mc_amp_nl_inv[m] / mc_amp_lin_inv[m]
  volatility_ratio_BGST[m] <-  mc_sd_nl_inv[m] / mc_sd_lin_inv[m]
  downturn_first_Q_ratio_BGST[m] <-    mc_Q1_nl_inv[m]/mc_Q1_lin_inv[m]
  max_crash_ratio_BGST[m] <-  mc_min_nl_inv[m] /   mc_min_lin_inv[m]
  ac1_ratio_BGST[m] <-  mc_ac1_nl_inv[m] / mc_ac1_lin_inv[m]
  max_eigenvalue_ratio_BGST[m] <-  true_eigen_dgp_3D[m] /  mc_eigen_lin_3D[m] 
  mean_ratio_BGST[m] <- mc_mean_nl_inv[m] / mc_mean_lin_inv[m]
  
}

# =====================================================================
# --- 4. PLOT CDFs (BGST) ---
# =====================================================================

# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

# x-range covering the middle 99% of runs (0.5th to 99.5th percentile),
# always including 1, with a margin of 5% of the range's width on each side
# With log = TRUE, the range and margin are computed
# on the log scale.
cdf_xlim <- function(x, log = FALSE) {
  x <- x[is.finite(x)]
  if (log) x <- x[x > 0]
  r <- quantile(x, c(0.005, 0.995))
  r <- range(r, 1)
  if (log) {
    lr <- log10(r)
    10^(lr + c(-1, 1) * 0.05 * diff(lr))
  } else {
    r + c(-1, 1) * 0.05 * diff(r)
  }
}

# Decide whether a variable needs a log axis. A log axis is used only if
#  (a) all values are positive (a log axis cannot show zero or negatives), and
#  (b) on a regular axis the middle 80% of runs (10th-90th percentile) would
#      fill less than 35% of the axis width, i.e. long tails squash the plot.
needs_log <- function(x, min_share = 0.35) {
  x <- x[is.finite(x)]
  if (any(x <= 0)) return(FALSE)
  r <- range(quantile(x, c(0.005, 0.995)), 1)
  diff(quantile(x, c(0.10, 0.90))) / diff(r) < min_share
}

# Mark the CDF at x = 1: vertical line from the x-axis up to the curve,
# horizontal line from the curve to the y-axis, and the CDF value printed
mark_at_one <- function(x) {
  F1  <- ecdf(x)(1)
  usr <- par("usr")
  x_left <- if (par("xlog")) 10^usr[1] else usr[1]
  segments(x0 = 1, y0 = usr[3], x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  segments(x0 = x_left, y0 = F1, x1 = 1, y1 = F1, col = "blue", lty = 2, lwd = 1.5)
  text(x = x_left, y = F1, labels = sprintf("%.2f", F1),
       adj = c(-0.15, -0.5), col = "blue", cex = CEX_MARK)
}

# Draw the x-axis ourselves: R's usual tick positions plus a tick at 1 that
# is always labelled. Ticks whose labels would overlap the "1" label are
# dropped, and the "1" tick is drawn in its own call so R never omits it.
# Labels are written as plain numbers (no scientific notation).
draw_x_axis <- function(log = FALSE) {
  fmt <- function(v) sapply(v, function(z) format(z, scientific = FALSE, trim = TRUE,
                                                  drop0trailing = TRUE, big.mark = ","))
  ticks  <- axTicks(1)
  ticks  <- ticks[abs(ticks - 1) > 1e-9]
  pos_in <- function(v) grconvertX(v, "user", "inches")
  w_in   <- function(lab) strwidth(lab, units = "inches", cex = CEX_AXIS)
  keep   <- abs(pos_in(ticks) - pos_in(1)) > 0.6 * (w_in(fmt(ticks)) + w_in("1")) + 0.05
  ticks  <- ticks[keep]
  if (length(ticks)) axis(1, at = ticks, labels = fmt(ticks), cex.axis = CEX_AXIS)
  axis(1, at = 1, labels = "1", cex.axis = CEX_AXIS)
}

# Draw one CDF plot; the log axis is chosen automatically by needs_log()
plot_cdf <- function(x, file, xlab) {
  lg <- needs_log(x)
  cat(sprintf("%-45s %s\n", file, if (lg) "log x-axis" else "regular x-axis"))
  open_png(file)
  plot(ecdf(x), xlim = cdf_xlim(x, log = lg), log = if (lg) "x" else "", xaxt = "n",
       xlab = "", ylab = "", main = "", lwd = 2.5,
       cex.lab = CEX_LAB, cex.axis = CEX_AXIS)
  draw_x_axis(log = lg)
  title_x(xlab)
  mark_at_one(x)
  dev.off()
}

# ---------------------------------------------------------------------
# The seven CDF plots
# ---------------------------------------------------------------------

# peak-trough movement ratio CDF
plot_cdf(peak_trough_movement_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_peak_trough_movement_ratio_CDF.png",
         xlab2("Peak-Trough Movement"))

# volatility ratio CDF
plot_cdf(volatility_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_volatility_ratio_CDF.png",
         xlab2("Volatility"))

# downturn first quartile ratio CDF
plot_cdf(downturn_first_Q_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_downturn_first_Q_ratio_CDF.png",
         xlab2("Downturn Severity (First Quartile)"))

# max crash ratio CDF
plot_cdf(max_crash_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_max_crash_ratio_CDF.png",
         xlab2("Downturn Severity (Max Crash)"))

# ac1 ratio CDF
plot_cdf(ac1_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_ac1_ratio_CDF.png",
         xlab2("Autocorrelation (lag 1)"))

# max eig val ratio CDF
plot_cdf(max_eigenvalue_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_max_eig_val_ratio_CDF.png",
         xlab2("Max Eigenvalue (Modulus)"))

# unconditional mean ratio CDF
plot_cdf(mean_ratio_BGST, "~/Desktop/Applied_EL_Revisions/BGST/BGST_mean_ratio_CDF.png",
         xlab2("Unconditional Mean"))
