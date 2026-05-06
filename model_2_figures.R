#Command line arguments
args<-commandArgs(TRUE)

library(mvtnorm)

#Samuelson 1939/Hicks 1950
# =====================================================================
# --- 1. Define Parameters (SM39 Bivariate System) ---
# =====================================================================
y_max <- 100
I_min <- -50
p <- 1 # 1 lag in the system

Burn_In <- 600     # Discard transient phase
Num_Obs <- 300     # Periods to plot/analyze
Total_Sim <- Burn_In + Num_Obs

# Multivariate Shock Covariance Matrix
eps_var <- matrix(c(0.0001, 0.00,  
                    0.00, 0.0001), 
                  nrow = 2, byrow = TRUE)

# =====================================================================
# --- 2. SINGLE REALIZATION & PLOTTING ---
# =====================================================================
set.seed(42) 

# Generate correlated multivariate shocks
shocks <- t(rmvnorm(Total_Sim + p-1 , mean = c(0, 0), sigma = eps_var))

y_full <- numeric(Total_Sim + p)
c_full <- numeric(Total_Sim + p)

# Initial conditions
y_full[1] <- 0
c_full[1] <- 0

beta_param <- 2 
c_param <- as.numeric(args[1])

# Simulate True Bivariate Nonlinear DGP
for (t in (p+1):(Total_Sim + p)) {

  c_full[t] <- c_param*y_full[t-1] + shocks[2, t-1]
  I_tilde <- beta_param*(c_full[t]-c_full[t-1])
  y_full[t] <- min(max(I_tilde,I_min) + c_full[t] + shocks[1, t-1], y_max)  
  
}

series_y <- ts(y_full[(Burn_In + p + 1):(Total_Sim + p)])
series_c <- ts(c_full[(Burn_In + p + 1):(Total_Sim + p)])

# Estimate Linear VAR(1) Model
d1 <- ts.intersect(series_y, series_c, L1.y=stats::lag(series_y,-1), L1.c=stats::lag(series_c,-1))
r1 <- lm(series_y ~ L1.y + L1.c, data = d1)
r2 <- lm(series_c ~ L1.y + L1.c, data = d1)

# Extract Coefficients
c_y  <- coef(r1)["(Intercept)"]
b_yy <- coef(r1)["L1.y"]
b_yc <- coef(r1)["L1.c"]

c_c  <- coef(r2)["(Intercept)"]
b_cy <- coef(r2)["L1.y"]
b_cc <- coef(r2)["L1.c"]

# Simulate Linear Bivariate Forecast after Burn-in Period (with Concurrent Noise)
y_lin <- numeric(Num_Obs)
c_lin <- numeric(Num_Obs)
y_lin[1] <- series_y[1] 
c_lin[1] <- series_c[1]

for (t in (p+1):Num_Obs) {
  y_lin[t] <- c_y + b_yy*y_lin[t-1] + b_yc*c_lin[t-1] + shocks[1, Burn_In+p + t-1]
  c_lin[t] <- c_c + b_cy*y_lin[t-1] + b_cc*c_lin[t-1] + shocks[2, Burn_In+p + t-1]
}

# Plotting
png(paste("SM39_Bivariate_Plot_HD_c=",c_param,".png",sep=""), width = 2500, height = 3000, res = 300)
par(mfrow = c(2, 1), mar = c(5, 5, 4, 2), cex = 1.1, cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1) 

Zoom_Window <- 1:200

series_y_zoom <- as.numeric(series_y)[Zoom_Window]
series_c_zoom <- as.numeric(series_c)[Zoom_Window]
y_lin_zoom <- y_lin[Zoom_Window]
c_lin_zoom <- c_lin[Zoom_Window]

# ==========================================================
# Plot 1: Output Dynamics (Nonlinear vs Linear)
# ==========================================================
y_range_out <- range(c(series_y_zoom, y_lin_zoom), na.rm = TRUE)
y_range_out[1] <- y_range_out[1] - diff(y_range_out) * 0.30 

plot(Zoom_Window,series_y_zoom, type="l", col="blue", lwd=1, ylim=y_range_out,
     ylab = "Output Value", xlab = "Time", main = "Output (Y): True DGP vs. Linear VAR(1)")
lines(Zoom_Window,y_lin_zoom, col="red", lwd=1, lty=2)

legend("bottom", legend=c("Nonlinear Dynamics", "Linear Dynamics"), 
       col=c("blue", "red"), lty=c(1, 2), lwd=1, horiz=TRUE, bty="n", cex=1)

# ==========================================================
# Plot 2: Consumption Dynamics (Nonlinear vs Linear)
# ==========================================================
y_range_cons <- range(c(series_c_zoom, c_lin_zoom), na.rm = TRUE)
y_range_cons[1] <- y_range_cons[1] - diff(y_range_cons) * 0.30

plot(Zoom_Window,series_c_zoom, type="l", col="blue", lwd=1, ylim=y_range_cons,
     ylab = "Consumption Value", xlab = "Time", main = "Consumption (C): True DGP vs. Linear VAR(1)")
lines(Zoom_Window,c_lin_zoom, col="red", lwd=1, lty=2) 

legend("bottom", legend=c("Nonlinear Dynamics", "Linear Dynamics"), 
       col=c("darkgreen", "orange"), lty=c(1, 2), lwd=1, horiz=TRUE, bty="n", cex=1)

dev.off()
