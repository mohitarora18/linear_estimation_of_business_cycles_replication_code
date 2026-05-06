#Command line arguments
args<-commandArgs(TRUE)

library(mvtnorm)

#Three-dimensional BGST model (no limit cycle behavior)
# =====================================================================
# --- 1. Define Parameters (BGST 3D System) ---
# =====================================================================
# Model parameters
rho      <- 0.95   
theta    <- 0.913    
eta      <- as.numeric(args[1]) 
phi_inv  <- 0.25  
phi_cost <- 0.67    
repay    <- 0.1   
gamma    <- 0.5
alpha    <- 2.5

# Shock Covariance Matrix
eps_var <- matrix(c(0.0001, 0.00, 0.00,   
                    0.00, 0.0001, 0.00,  
                    0.00, 0.00, 0.0001), 
                  nrow = 3, byrow = TRUE)

# Simulation Controls
p <- 1             # Lags in the estimated VAR(1)
Burn_In <- 600     # Discard transient phase
Num_Obs <- 300     # Periods to analyze/plot
Total_Sim <- Burn_In + Num_Obs

# Helper function to extract Jacobian blocks 
get_lag_matrix <- function(lag_num, m_TFP, m_Debt, m_Inv) {
  vars <- c("TFP", "Debt", "Invest")
  cols <- paste0("series.L", lag_num, ".", vars)
  rbind(coef(m_TFP)[cols], coef(m_Debt)[cols], coef(m_Inv)[cols])
}


# =====================================================================
# --- 3. SINGLE REALIZATION & PLOTTING ---
# =====================================================================
set.seed(456)

# A. Generate Shocks & Simulate Full Sequence
shocks <- t(rmvnorm(Total_Sim + p-1, mean=c(0,0,0), sigma = eps_var))
x_NL <- matrix(0, nrow = 3, ncol = Total_Sim + p)

for(t in (p+1):(Total_Sim + p)){
  news <- shocks[1,t-1] 
  x_NL[1,t] <- rho * x_NL[1,t-1] + news
  belief <- rho * x_NL[1,t] + theta * news
  dist_to_default <- x_NL[2,t-1] - x_NL[1,t]
  spread_NL <- alpha * exp(eta * dist_to_default)
  x_NL[3,t] <- phi_inv * belief - phi_cost * spread_NL + shocks[3,t-1]
  x_NL[2,t] <- (1-repay) * x_NL[2,t-1] + gamma * x_NL[3,t] + shocks[2,t-1]
}

# B. Convert to multivariate ts, then safely subset to remove burn-in period using window()
ts_full <- ts(t(x_NL))
colnames(ts_full) <- c("TFP", "Debt", "Invest")
series <- window(ts_full, start = Burn_In + p + 1) # THIS PRESERVES THE TS CLASS

# C. Estimate Linear VAR(1) Model safely with ts.intersect
d1 <- ts.intersect(series, series.L1=stats::lag(series,-1))

r_TFP <- lm(series.TFP ~ series.L1.TFP + series.L1.Debt + series.L1.Invest, data = d1)
r_Debt <- lm(series.Debt ~ series.L1.TFP + series.L1.Debt + series.L1.Invest, data = d1)
r_Inv <- lm(series.Invest ~ series.L1.TFP + series.L1.Debt + series.L1.Invest, data = d1)

# D. Build the 3x3 Jacobian
Jacobian <- get_lag_matrix(1, r_TFP, r_Debt, r_Inv)

# E. Simulate Linear Forecast (With Concurrent Noise)
c_real <- c(coef(r_TFP)["(Intercept)"], coef(r_Debt)["(Intercept)"], coef(r_Inv)["(Intercept)"])
States_Linear <- matrix(0, nrow = 3, ncol = Num_Obs)
States_Linear[, 1] <- series[1, ]

for (t in (p+1):Num_Obs) {
  States_Linear[,t] <- c_real + Jacobian %*% States_Linear[,t-1]
  States_Linear[, t] <- States_Linear[, t] + shocks[, Burn_In+p + t-1]
}



# F. Plotting (Focused on Investment)
png(paste("BGST_Invest_Plot_HD_eta=",eta,".png",sep=""), width = 2500, height = 1500, res = 300)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex = 1.1, cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1) 

Zoom_Window <- 1:300
Inv_Nonlinear_zoom <- as.numeric(series[Zoom_Window, "Invest"])
Inv_Linear_zoom    <- States_Linear[3, Zoom_Window]

y_lims <- range(c(Inv_Nonlinear_zoom, Inv_Linear_zoom), na.rm = TRUE)
y_lims[1] <- y_lims[1] - diff(y_lims) * 0.30 

plot(Zoom_Window,Inv_Nonlinear_zoom, type="l", col="blue", lwd=2, ylim=y_lims, 
     ylab="Investment Level", xlab="Time")

lines(Inv_Linear_zoom, col="red", lwd=2, lty=2)

legend("bottom", legend = c("Nonlinear Dynamics", "Linear Dynamics"), 
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, horiz = TRUE, bty = "n", cex = 1)
dev.off()
