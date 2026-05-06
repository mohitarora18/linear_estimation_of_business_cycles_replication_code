rm(list=ls())
options(scipen = 999)
# =====================================================================
# --- 1. Define Parameters (Beaudry 2017) ---
# =====================================================================
a1 <- 1.5      
a2 <- -0.6
a3 <- -0.3
a_cub <- -0.01 
sigma_eps <- 0.25
p <- 3

Burn_In <- 600     # Number of periods to discard (transient phase)
Num_Obs <- 300     # Number of periods to actually plot/analyze
Total_Sim <- Burn_In + Num_Obs

# =====================================================================
# --- 2. MONTE CARLO SIMULATION (1000 Iterations) ---
# =====================================================================
M <- 1000

# Initialize storage vectors for the metrics
mc_sd_nl_beaudry <- numeric(M)
mc_sd_lin_beaudry <- numeric(M)
mc_eigen_lin_beaudry <- numeric(M)
mc_min_nl_beaudry <- numeric(M)
mc_min_lin_beaudry <- numeric(M)
mc_Q1_nl_beaudry <- numeric(M); 
mc_Q1_lin_beaudry <- numeric(M); 

set.seed(123) # Set master seed for reproducibility of the MC loop

for (m in 1:M) {
  
  # A. Generate fresh noise for this iteration
  eps_m <- rnorm(Total_Sim + p)
  x_m <- numeric(Total_Sim + p)
  
  # B. Simulate true Nonlinear DGP
  for(i in (p+1):(Total_Sim+p)){
    x_m[i] <- a1*x_m[i-1] + a2*x_m[i-2] + a3*x_m[i-3] + a_cub*(x_m[i-1]^3) + sigma_eps*eps_m[i]
  }
  
  # Isolate post-burn-in
  series_m <- x_m[(Burn_In + p + 1):(Total_Sim + p)]
  ts_m <- ts(series_m)
  
  # C. Estimate Linear AR(3) on the isolated series
  d_m <- ts.intersect(ts_m, ts.L1=stats::lag(ts_m,-1), ts.L2=stats::lag(ts_m,-2), ts.L3=stats::lag(ts_m,-3))
  lin_mod_m <- lm(ts_m ~ ts.L1 + ts.L2 + ts.L3, data=d_m)
  
  c_m  <- coef(lin_mod_m)[1]
  a1_m <- coef(lin_mod_m)[2]
  a2_m <- coef(lin_mod_m)[3]
  a3_m <- coef(lin_mod_m)[4]
  
  # D. Calculate Max Eigenvalue using a Companion Matrix
  Comp_Mat <- rbind(c(a1_m, a2_m, a3_m),
                    c(1, 0, 0),
                    c(0, 1, 0))
  mc_eigen_lin_beaudry[m] <- max(Mod(eigen(Comp_Mat)$values))
  
  # E. Simulate Linear Forward with concurrent noise
  x_lin_m <- numeric(Num_Obs)
  x_lin_m[1:3] <- series_m[1:3] # Start from same condition
  for (t in 4:Num_Obs) {
    x_lin_m[t] <- c_m + a1_m*x_lin_m[t-1] + a2_m*x_lin_m[t-2] + a3_m*x_lin_m[t-3] + sigma_eps*eps_m[Burn_In + p + t]
  }
  
  # F. Calculate & Store the metrics
  # 1. Volatility (Standard Deviation)
  mc_sd_nl_beaudry[m] <- sd(series_m)
  mc_sd_lin_beaudry[m] <- sd(x_lin_m)
  
  # 2. Maximum Crash (Minimum Value)
  mc_min_nl_beaudry[m] <- min(series_m)
  mc_min_lin_beaudry[m] <- min(x_lin_m)
  
  # 3. First quatrile
  mc_Q1_nl_beaudry[m] <- as.numeric(quantile(series_m,0.25))
  mc_Q1_lin_beaudry[m] <- as.numeric(quantile(x_lin_m,0.25))
}

# =====================================================================
# --- 4. PRINT OUT AVERAGED METRICS (Beaudry) ---
# =====================================================================
Comp_Mat_DGP <- rbind(c(a1, a2, a3),
                      c(1, 0, 0),
                      c(0, 1, 0))
true_eigen_dgp_beaudry <- max(Mod(eigen(Comp_Mat_DGP)$values))

cat("--- BEAUDRY MODEL ---\n")
cat("   Linear AR(3) Volatility: ", round(mean(mc_sd_lin_beaudry), 4), "(", round(sd(mc_sd_lin_beaudry), 4), ")\n")
cat("   Linear AR(3) Q1: ", round(mean(mc_Q1_lin_beaudry), 4), "(", round(sd(mc_Q1_lin_beaudry), 4), ")\n")
cat("   Linear AR(3) Min: ", round(mean(mc_min_lin_beaudry), 4), "(", round(sd(mc_min_lin_beaudry), 4), ")\n\n")


#SM39
library(mvtnorm)

# =====================================================================
# --- 1. Define Parameters (SM39 Bivariate System) ---
# =====================================================================
y_max <- 100
I_min <- -50
c <- 0.9
beta <- 2
p <- 1 # 1 lag in the system

Burn_In <- 600     # Discard transient phase
Num_Obs <- 300     # Periods to plot/analyze
Total_Sim <- Burn_In + Num_Obs

# Multivariate Shock Covariance Matrix
eps_var <- matrix(c(0.0001, 0.00,  
                    0.00, 0.0001), 
                  nrow = 2, byrow = TRUE)

# =====================================================================
# --- 2. MONTE CARLO SIMULATION (1000 Iterations) ---
# =====================================================================
M <- 1000

# Initialize system storage
mc_eigen_lin_2D <- numeric(M)
# Initialize variables storage (only Output 'y' is retained for metrics)
mc_sd_nl_y <- numeric(M); 
mc_sd_lin_y <- numeric(M); 

mc_min_nl_y <- numeric(M); 
mc_min_lin_y <- numeric(M); 

mc_Q1_nl_y <- numeric(M); 
mc_Q1_lin_y <- numeric(M); 


set.seed(132)

for (m in 1:M) {
  
  # A. Generate multivariate noise
  shocks_m <- t(rmvnorm(Total_Sim + p, mean = c(0, 0), sigma = eps_var))
  
  ym <- numeric(Total_Sim + p)
  cm <- numeric(Total_Sim + p)
  ym[1] <- 0; cm[1] <- 0
  
  # B. Simulate DGP
  for (t in 2:(Total_Sim + p)) { 
    
    cm[t] <- c*ym[t-1] + shocks_m[2, t] 
    I_tilde <- beta*(cm[t]-cm[t-1]) 
    ym[t] <- min(max(I_tilde,I_min) + cm[t] + shocks_m[1, t], y_max) 
    
  }

  
  # C. Extract & Estimate
  ts_y <- ts(ym[(Burn_In + p + 1):(Total_Sim + p)])
  ts_c <- ts(cm[(Burn_In + p + 1):(Total_Sim + p)])
  
  d_m <- ts.intersect(ts_y, ts_c, L1.y=stats::lag(ts_y,-1), L1.c=stats::lag(ts_c,-1))
  r1_m <- lm(ts_y ~ L1.y + L1.c, data = d_m)
  r2_m <- lm(ts_c ~ L1.y + L1.c, data = d_m)
  
  # D. VAR(1) Jacobian & Eigenvalue
  B_mat <- matrix(c(coef(r1_m)["L1.y"], coef(r1_m)["L1.c"],
                    coef(r2_m)["L1.y"], coef(r2_m)["L1.c"]), 
                  nrow=2, ncol=2, byrow=TRUE)
  mc_eigen_lin_2D[m] <- max(Mod(eigen(B_mat)$values))
  
  # E. Simulate Linear Forward
  ylin_m <- numeric(Num_Obs); clin_m <- numeric(Num_Obs)
  ylin_m[1] <- ts_y[1]; clin_m[1] <- ts_c[1]
  
  for (t in 2:Num_Obs) {
    ylin_m[t] <- coef(r1_m)[1] + coef(r1_m)[2]*ylin_m[t-1] + coef(r1_m)[3]*clin_m[t-1] + shocks_m[1, Burn_In + p + t]
    clin_m[t] <- coef(r2_m)[1] + coef(r2_m)[2]*ylin_m[t-1] + coef(r2_m)[3]*clin_m[t-1] + shocks_m[2, Burn_In + p + t]
  }
  
  # F. Metrics Calculation (Output only)
  mc_sd_nl_y[m] <- sd(ts_y)
  mc_sd_lin_y[m] <- sd(ylin_m)
  mc_min_nl_y[m] <- min(ts_y)
  mc_min_lin_y[m] <- min(ylin_m)  
  mc_Q1_nl_y[m] <- as.numeric(quantile(ts_y,0.25))
  mc_Q1_lin_y[m] <- as.numeric(quantile(ylin_m,0.25))
}

# =====================================================================
# --- 4. PRINT AVERAGED METRICS (2D System) ---
# =====================================================================

B_true <- matrix(c((c * beta) + c, -beta, 
                   c,              0), 
                 nrow=2, ncol=2, byrow=TRUE)
true_eigen_dgp_2D <- max(Mod(eigen(B_true)$values))

cat("--- 2D SM39 MODEL ---\n")
cat("   [Series Y] Linear VAR(1) Volatility: ", round(mean(mc_sd_lin_y), 4), "(", round(sd(mc_sd_lin_y), 4), ")\n")
cat("   [Series Y] Linear VAR(1) Q1: ", round(mean(mc_Q1_lin_y), 4), "(", round(sd(mc_Q1_lin_y), 4), ")\n")
cat("   [Series Y] Linear VAR(1) Min: ", round(mean(mc_min_lin_y), 4), "(", round(sd(mc_min_lin_y), 4), ")\n\n")


#Three-dimensional BGST model
library(mvtnorm)

# =====================================================================
# --- 1. Define Parameters (BGST 3D System) ---
# =====================================================================
rho      <- 0.95   
theta    <- 0.913    
eta      <- 8    
phi_inv  <- 0.25    
phi_cost <- 0.67    
repay    <- 0.1   
gamma    <- 0.5
alpha <- 2.5
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

true_eigen_dgp_3D <- max(Mod(eigen(J_true)$values))

# =====================================================================
# --- 3. MONTE CARLO SIMULATION (1000 Iterations) ---
# =====================================================================
M <- 1000

# Storage
mc_eigen_lin_3D <- numeric(M)
mc_sd_nl_inv <- numeric(M); mc_sd_lin_inv <- numeric(M)
mc_min_nl_inv <- numeric(M); mc_min_lin_inv <- numeric(M)
mc_Q1_nl_inv <- numeric(M); mc_Q1_lin_inv <- numeric(M)

set.seed(456)

for (m in 1:M) {
  
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
  
}

cat("--- 3D BGST MODEL ---\n")
cat("   Linear VAR(1) Investment Volatility: ", round(mean(mc_sd_lin_inv), 4), "(", round(sd(mc_sd_lin_inv), 4), ")\n")
cat("   Linear VAR(1) Q1: ", round(mean(mc_Q1_lin_inv), 4), "(", round(sd(mc_Q1_lin_inv), 4), ")\n")
cat("   Linear VAR(1) Min: ", round(mean(mc_min_lin_inv), 4), "(", round(sd(mc_min_lin_inv), 4), ")\n\n")


# =====================================================================
# --- 6. COMPILE HIGH-QUALITY LATEX TABLE (Means and SDs) ---
# =====================================================================
# Install if necessary: install.packages("kableExtra")
library(kableExtra)

# 1. BUILD THE DATA FRAME (Explicitly writing out Mean and SD with paste0 and round)
results_df <- data.frame(
  Metric = c("Volatility (Std. Dev.)", 
             "Downturn Severity (First Quartile)",
             "Downturn Severity (Max Crash)", 
             "Max Eigenvalue (Modulus)"),
  
  # Model 1: Beaudry AR(3) System 
  Beaudry_Lin = c(paste0(round(mean(mc_sd_lin_beaudry), 4), " (", round(sd(mc_sd_lin_beaudry), 4), ")"), 
 		  paste0(round(mean(mc_Q1_lin_beaudry), 4), " (", round(sd(mc_Q1_lin_beaudry), 4), ")"), 
                  paste0(round(mean(mc_min_lin_beaudry), 4), " (", round(sd(mc_min_lin_beaudry), 4), ")"), 
                  paste0(round(mean(mc_eigen_lin_beaudry), 4), " (", round(sd(mc_eigen_lin_beaudry), 4), ")")),
  
  Beaudry_NL  = c(paste0(round(mean(mc_sd_nl_beaudry), 4), " (", round(sd(mc_sd_nl_beaudry), 4), ")"), 
  		  paste0(round(mean(mc_Q1_nl_beaudry), 4), " (", round(sd(mc_Q1_nl_beaudry), 4), ")"),
                  paste0(round(mean(mc_min_nl_beaudry), 4), " (", round(sd(mc_min_nl_beaudry), 4), ")"), 
                  round(true_eigen_dgp_beaudry, 4)),
  
  # Model 2: 2D System (Reporting OUTPUT 'y' ONLY)
  TwoD_Out_Lin = c(paste0(round(mean(mc_sd_lin_y), 4), " (", round(sd(mc_sd_lin_y), 4), ")"),
                   paste0(round(mean(mc_Q1_lin_y), 4), " (", round(sd(mc_Q1_lin_y), 4), ")"),  
                   paste0(round(mean(mc_min_lin_y), 4), " (", round(sd(mc_min_lin_y), 4), ")"), 
                   paste0(round(mean(mc_eigen_lin_2D), 4), " (", round(sd(mc_eigen_lin_2D), 4), ")")),
  
  TwoD_Out_NL  = c(paste0(round(mean(mc_sd_nl_y), 4), " (", round(sd(mc_sd_nl_y), 4), ")"), 
                   paste0(round(mean(mc_Q1_nl_y), 4), " (", round(sd(mc_Q1_nl_y), 4), ")"), 
                   paste0(round(mean(mc_min_nl_y), 4), " (", round(sd(mc_min_nl_y), 4), ")"), 
                   round(true_eigen_dgp_2D, 4)),
  
  # Model 3: Three-Variable System (Reporting INVESTMENT ONLY)
  ThreeVar_Lin = c(paste0(round(mean(mc_sd_lin_inv), 4), " (", round(sd(mc_sd_lin_inv), 4), ")"), 
                   paste0(round(mean(mc_Q1_lin_inv), 4), " (", round(sd(mc_Q1_lin_inv), 4), ")"), 
                   paste0(round(mean(mc_min_lin_inv), 4), " (", round(sd(mc_min_lin_inv), 4), ")"), 
                   paste0(round(mean(mc_eigen_lin_3D), 4), " (", round(sd(mc_eigen_lin_3D), 4), ")")),
  
  ThreeVar_NL  = c(paste0(round(mean(mc_sd_nl_inv), 4), " (", round(sd(mc_sd_nl_inv), 4), ")"), 
                   paste0(round(mean(mc_Q1_nl_inv), 4), " (", round(sd(mc_Q1_nl_inv), 4), ")"), 
                   paste0(round(mean(mc_min_nl_inv), 4), " (", round(sd(mc_min_nl_inv), 4), ")"), 
                   round(true_eigen_dgp_3D, 4))
)


# 2. COMPILE INTO LATEX CODE
latex_table <- kbl(results_df, 
                   format = "latex", 
                   booktabs = TRUE, 
                   caption = "Monte Carlo Simulation Results: Linear Estimation vs. Nonlinear DGP (Mean and Std. Deviation)",
                   label = "mc_results", 
                   col.names = c("Metric", 
                                 "Linear", "Nonlinear", 
                                 "Linear", "Nonlinear", 
                                 "Linear", "Nonlinear"),
                   align = c("l", rep("c", 6))) %>%
  
  # BOTTOM HEADER ROW: Specifies the Variable
  add_header_above(c(" " = 1, 
                     "$x_t$" = 2, 
                     "$Y_t$" = 2, 
                     "$Inv_t$" = 2)) %>%
  
  # TOP HEADER ROW: Specifies the Model
  add_header_above(c(" " = 1, 
                     "beaudry2017macroeconomy" = 2, 
                     "sm39/hk50" = 2, 
                     "bordalo2026real" = 2)) %>%
  
  # Make the Metric names bold for readability
  column_spec(1, bold = TRUE) %>%
  
  # Automatically scale down if it exceeds page width in Overleaf
  kable_styling(latex_options = c("hold_position", "scale_down"))

# Print to console so you can copy/paste directly into Overleaf
cat(latex_table)
