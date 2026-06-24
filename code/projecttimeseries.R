###############################################################################
# MATH1318 Final Project
# Modelling Australian Real Residential Property Prices with ARIMAX
# Data: FRED/BIS, Q1 1978 – Q4 2025, quarterly
###############################################################################

## 0. Packages and Helpers ----------------------------------------------------

# Install once if needed:
#install.packages(c("TSA","forecast","tseries","lmtest","FitAR"))
install.packages(c("fUnitRoots"))
library(TSA)
library(forecast)
library(tseries)
library(lmtest)
library(FitAR)
library(fUnitRoots)
setwd("/Users/kritiyadav/Downloads")

# Source the course helper functions (sort.score, residual.analysis, etc.)
source("Overview_RCodes_EndofSem2026.R")   # Make sure this file is in your working dir

## 1. Data Loading and Preparation --------------------------------------------

# Read FRED CSVs (must be in the working directory)
rppi_raw  <- read.csv("QAUR628BIS.csv",
                      header = TRUE, col.names = c("date","rppi"))
bond_raw  <- read.csv("IRLTLT01AUQ156N.csv",
                      header = TRUE, col.names = c("date","bond10"))
unemp_raw <- read.csv("LRHUTTTTAUM156S.csv",
                      header = TRUE, col.names = c("date","unemp"))

# Convert dates
rppi_raw$date  <- as.Date(rppi_raw$date)
bond_raw$date  <- as.Date(bond_raw$date)
unemp_raw$date <- as.Date(unemp_raw$date)

# Convert monthly unemployment to quarterly average (start of quarter)
unemp_raw$quarter <- as.Date(paste0(
  format(unemp_raw$date, "%Y-"),
  sprintf("%02d",
          ceiling(as.numeric(format(unemp_raw$date, "%m")) / 3) * 3 - 2),
  "-01"
))
unemp_q <- aggregate(unemp ~ quarter, data = unemp_raw, FUN = mean)
names(unemp_q) <- c("date","unemp")

# Merge to RPPI dates and keep complete cases
df <- merge(rppi_raw, bond_raw, by = "date", all.x = TRUE)
df <- merge(df,       unemp_q,  by = "date", all.x = TRUE)
df <- df[complete.cases(df), ]
df <- df[order(df$date), ]

cat("Dataset:", nrow(df), "quarterly observations\n")
cat("Range  :", format(min(df$date)), "to", format(max(df$date)), "\n\n")

# Create time-series objects (quarterly, start 1978 Q1)
rppi_ts  <- ts(df$rppi,   start = c(1978,1), frequency = 4)
bond_ts  <- ts(df$bond10, start = c(1978,1), frequency = 4)
unemp_ts <- ts(df$unemp,  start = c(1978,1), frequency = 4)

# Regressor matrix for ARIMAX (levels)
xreg_full <- cbind(bond10 = df$bond10,
                   unemp  = df$unemp)

## 2. Descriptive Analysis ----------------------------------------------------

cat("=== Descriptive statistics ===\n")
print(summary(df[, c("rppi","bond10","unemp")]))
cat("\n")

# Time-series plots
par(mfrow = c(3,1), mar = c(3,4,3,1))

plot(rppi_ts, type = "o", pch = 16, cex = 0.4, col = "#2166AC",
     main = "Figure 1. Real RPPI, Australia (Q1 1978–Q4 2025)",
     ylab = "Index (2010 = 100)", xlab = "Year")
abline(v = c(1990,2008,2020), lty = 2, col = "grey50")

plot(bond_ts, type = "o", pch = 16, cex = 0.4, col = "#D6604D",
     main = "Figure 2. 10-year Government Bond Yield (Q1 1978–Q4 2025)",
     ylab = "Yield (%)", xlab = "Year")

plot(unemp_ts, type = "o", pch = 16, cex = 0.4, col = "#1A9850",
     main = "Figure 3. Unemployment Rate (Q1 1978–Q4 2025)",
     ylab = "Rate (%)", xlab = "Year")

par(mfrow = c(1,1))

# ACF/PACF of RPPI
par(mfrow = c(1,2))
acf(rppi_ts,  lag.max = 40, main = "Figure 4. ACF of Real RPPI")
pacf(rppi_ts, lag.max = 40, main = "Figure 5. PACF of Real RPPI")
par(mfrow = c(1,1))

# Scatterplots vs regressors
par(mfrow = c(1,2))
plot(df$bond10, df$rppi,
     main = "Figure 6. RPPI vs Bond Yield",
     xlab = "Bond Yield (%)", ylab = "RPPI (2010 = 100)",
     pch = 16, col = "#D6604D", cex = 0.7)
abline(lm(rppi ~ bond10, data = df), col = "darkred", lwd = 2)

plot(df$unemp, df$rppi,
     main = "Figure 7. RPPI vs Unemployment",
     xlab = "Unemployment (%)", ylab = "RPPI (2010 = 100)",
     pch = 16, col = "#1A9850", cex = 0.7)
abline(lm(rppi ~ unemp, data = df), col = "darkgreen", lwd = 2)
par(mfrow = c(1,1))

## 3. Stationarity and Transformation -----------------------------------------

cat("=== Stationarity tests: RPPI in levels ===\n")
print(adf.test(rppi_ts))
print(pp.test(rppi_ts))
print(kpss.test(rppi_ts))
cat("\n")

# Box-Cox transformation
lambda_bc <- BoxCox.lambda(rppi_ts)
cat("Optimal Box-Cox lambda:", round(lambda_bc,4), "\n\n")
rppi_bc <- BoxCox(rppi_ts, lambda = lambda_bc)

par(mfrow = c(1,2))
plot(rppi_ts, main = "Figure 8a. Original RPPI",
     ylab = "Index", xlab = "Year", col = "#2166AC", type = "l")
plot(rppi_bc, main = "Figure 8b. Box–Cox Transformed RPPI",
     ylab = "Transformed", xlab = "Year", col = "#4D9221", type = "l")
par(mfrow = c(1,1))

# Difference transformed series
rppi_diff <- diff(rppi_bc, differences = 1)

plot(rppi_diff, type = "o", pch = 16, cex = 0.4, col = "#2166AC",
     main = "Figure 9. First-Differenced BC-RPPI",
     ylab = "First difference", xlab = "Year")
abline(h = 0, lty = 2, col = "grey50")

cat("=== Stationarity tests: first-differenced BC-RPPI ===\n")
print(adf.test(rppi_diff))
print(pp.test(rppi_diff))
print(kpss.test(rppi_diff))
cat("\n")

# ACF/PACF of stationary series
par(mfrow = c(1,2))
acf(rppi_diff,  lag.max = 40, main = "Figure 10. ACF of Diff BC-RPPI")
pacf(rppi_diff, lag.max = 40, main = "Figure 11. PACF of Diff BC-RPPI")
par(mfrow = c(1,1))

## 4. Model Specification (ACF/PACF, EACF, BIC table) ------------------------

cat("=== EACF for differenced BC-RPPI ===\n")
eacf(rppi_diff, ar.max = 6, ma.max = 6)

res_arma <- armasubsets(y = rppi_diff, nar = 5, nma = 5,
                        y.name = "p", ar.method = "ols")
plot(res_arma, main = "Figure 12. BIC table for ARMA on Diff BC-RPPI")

# Candidate ARIMAX orders (p,d,q) based on EACF + BIC + ACF/PACF
# Example set:
#   (0,1,1), (1,1,1), (2,1,1), (1,1,2), (2,1,2)

## 5. Model Fitting (ARIMAX) --------------------------------------------------

fit_arimax <- function(p, d, q, ts_data, xreg_mat, bc_lambda) {
  m_css <- tryCatch(
    Arima(ts_data, order = c(p,d,q), xreg = xreg_mat,
          lambda = bc_lambda, method = "CSS"),
    error = function(e) NULL
  )
  m_ml  <- tryCatch(
    Arima(ts_data, order = c(p,d,q), xreg = xreg_mat,
          lambda = bc_lambda, method = "ML"),
    error = function(e) NULL
  )
  list(css = m_css, ml = m_ml)
}

cat("=== Fitting candidate ARIMAX models ===\n")
m011 <- fit_arimax(0,1,1, rppi_ts, xreg_full, lambda_bc)
m111 <- fit_arimax(1,1,1, rppi_ts, xreg_full, lambda_bc)
m211 <- fit_arimax(2,1,1, rppi_ts, xreg_full, lambda_bc)
m112 <- fit_arimax(1,1,2, rppi_ts, xreg_full, lambda_bc)
m212 <- fit_arimax(2,1,2, rppi_ts, xreg_full, lambda_bc)

models_ml <- list(m011$ml, m111$ml, m211$ml, m112$ml, m212$ml)
names(models_ml) <- c("ARIMAX(0,1,1)","ARIMAX(1,1,1)",
                      "ARIMAX(2,1,1)","ARIMAX(1,1,2)","ARIMAX(2,1,2)")

# AIC/BIC tables (using sort.score helper)
aic_tab <- AIC(models_ml[[1]], models_ml[[2]], models_ml[[3]],
               models_ml[[4]], models_ml[[5]])
bic_tab <- AIC(models_ml[[1]], models_ml[[2]], models_ml[[3]],
               models_ml[[4]], models_ml[[5]], k = log(length(rppi_ts)))

cat("\nTable 1. Models sorted by AIC\n")
print(sort.score(aic_tab, score = "aic"))

cat("\nTable 2. Models sorted by BIC\n")
print(sort.score(bic_tab, score = "bic"))

aic_vals <- sapply(models_ml, function(m) if (!is.null(m)) AIC(m) else NA)
bic_vals <- sapply(models_ml, function(m) if (!is.null(m)) BIC(m) else NA)

best_idx  <- which.min(bic_vals)
best_name <- names(models_ml)[best_idx]
best_model <- models_ml[[best_idx]]

cat("\nSelected final model (by BIC):", best_name, "\n\n")

cat("=== Parameter estimates (coeftest) for", best_name, "===\n")
print(coeftest(best_model))

## 6. Diagnostics (using course residual.analysis helper) ---------------------

cat("\n=== Residual diagnostics for", best_name, "===\n")
residual.analysis(model = best_model,
                  std   = TRUE,
                  lag.max = 30,
                  class = "ARIMA")

# Optional: diagnostics for a second model to compare, e.g. ARIMAX(1,1,1)
if (!is.null(m111$ml)) {
  cat("\n=== Residual diagnostics for ARIMAX(1,1,1) ===\n")
  residual.analysis(model = m111$ml,
                    std   = TRUE,
                    lag.max = 30,
                    class = "ARIMA")
}

## 7. Forecasting (10 quarters ahead) -----------------------------------------

# Hold regressors constant at last observed values (simple scenario)
last_bond  <- tail(df$bond10, 1)
last_unemp <- tail(df$unemp,  1)

xreg_future <- matrix(rep(c(last_bond, last_unemp), 10),
                      nrow = 10, ncol = 2, byrow = TRUE,
                      dimnames = list(NULL, c("bond10","unemp")))

cat("\nForecast regressors (held constant):\n")
cat("  Bond yield last value:", round(last_bond,3), "%\n")
cat("  Unemployment last value:", round(last_unemp,3), "%\n\n")

frc <- forecast(best_model, h = 10, xreg = xreg_future)

cat("Table 3. 10-quarter ahead forecasts (Q1 2026 – Q2 2028)\n")
print(frc)

cat("\nTable 4. Forecast accuracy metrics (in-sample)\n")
print(accuracy(frc))

plot(frc,
     main = paste("Figure 13. 10-Quarter Forecast –", best_name),
     ylab = "Real RPPI (2010 = 100)", xlab = "Year",
     col  = "#2166AC", fcol = "#D6604D", flwd = 2,
     shadecols = c("#FDDBC7","#F4A582"))
lines(fitted(best_model), col = "#4DAC26", lty = 2, lwd = 1.5)
legend("topleft",
       legend = c("Observed","Fitted","Forecast","80% CI","95% CI"),
       lty    = c(1,2,1,NA,NA),
       pch    = c(NA,NA,NA,15,15),
       col    = c("#2166AC","#4DAC26","#D6604D","#F4A582","#FDDBC7"),
       bty    = "n", cex = 0.8)

cat("\n=== Final model summary ===\n")
cat("Model:", best_name, "\n")
cat("AIC:", round(AIC(best_model),2),
    "BIC:", round(BIC(best_model),2), "\n\n")
print(summary(best_model))