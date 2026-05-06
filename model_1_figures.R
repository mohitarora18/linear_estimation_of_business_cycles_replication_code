#To clear the global environment
rm(list=ls())

#Beaudry 2017 (stochastic univariate system)
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
# --- 2. SINGLE REALIZATION & PLOTTING ---
# =====================================================================
set.seed(42) 
myeps <- rnorm(Total_Sim)
x_dgp_full <- numeric(Total_Sim + p)

for(i in (p+1):(Total_Sim+p)){
  x_dgp_full[i] <- a1*x_dgp_full[i-1] + a2*x_dgp_full[i-2] + a3*x_dgp_full[i-3] + a_cub*(x_dgp_full[i-1]^3) + sigma_eps*myeps[i-p]
}

# 1. Capture the series starting after the burn-in period
series <- ts(x_dgp_full[(Burn_In+p + 1):(Total_Sim + p)])

d1 <- ts.intersect(series, series.L1=lag(series,-1), series.L2=lag(series,-2), series.L3=lag(series,-3))
lin_mod <- lm(series ~ series.L1 + series.L2 + series.L3, data=d1)

c_est <- coef(lin_mod)["(Intercept)"]
a1_est <- coef(lin_mod)["series.L1"]
a2_est <- coef(lin_mod)["series.L2"]
a3_est <- coef(lin_mod)["series.L3"]

# 2. Linear forecast
x_lin <- numeric(Num_Obs)
x_lin[1:p] <- series[1:p] 

for (t in p+1:Num_Obs) {
  x_lin[t] <- c_est + a1_est*x_lin[t-1] + a2_est*x_lin[t-2] + a3_est*x_lin[t-3] + sigma_eps*myeps[Burn_In+ t]
}

y_min <- min(c(series, x_lin), na.rm = TRUE)
y_max <- max(c(series, x_lin), na.rm = TRUE)

# =====================================================================
# --- 5. Plotting (HD Specifications) ---
# =====================================================================

png("Beaudry_Univariate_Plot_HD.png", width = 2500, height = 1500, res = 300)

par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex = 1.1, cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1) 

# 3. Plot results
Zoom_Window <- 1:150

series_zoom <- as.numeric(series)[Zoom_Window]
x_lin_zoom <- as.numeric(x_lin)[Zoom_Window]

# ==========================================================
# Plot: True Nonlinear DGP vs Estimated Linear AR(3)
# ==========================================================
y_range <- range(c(series_zoom, x_lin_zoom), na.rm = TRUE)
# Expand bottom by 30% for the legend
y_range[1] <- y_range[1] - diff(y_range) * 0.30 

# Plot the nonlinear model
plot(Zoom_Window,series_zoom, type="l", col="blue", lwd=1, ylim=y_range,
     ylab = "Value of x", xlab = "Time")

# Plot the linear forecast model
lines(Zoom_Window,x_lin_zoom, col="red", lwd=1, lty=2)

# Horizontal legend at the bottom with no bounding box
legend("bottomleft", legend=c("Nonlinear Dynamics", "Linear Dynamics"), 
       col=c("blue", "red"), lty=c(1, 2), lwd=1, horiz=TRUE, bty="n", cex=1)

dev.off()
