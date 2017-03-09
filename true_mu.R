# -----------------------------------------------
# code to find the true value of mu*(1, 0) via
# numerical integration for a given mean vector
# and covariance matrix used in data generation
# -----------------------------------------------

library(R2Cuba)

source("functions.R")

# correlation matrix for C', X, Y
m <- matrix(c(1.000, -0.405, -0.080,
             -0.405,  1.000,  0.700,
             -0.080,  0.700,  1.000), 3, 3)

p <- 0.20                      # Pr(C = 1)
mu <- c(0, 10, 16)             # mean vector
s <- sqrt(c(1, 2, 3))          # standard deviations
V <- diag(s) %*% m %*% diag(s) # covariance matrix

x.coef <- xprobit(V, mu, p)   # Coefficents for Pr(C = 1 | X)
xy.coef <- xyprobit(V, mu, p) # Coefficients for Pr(C = 1 | X, Y)

# integrate over f(Y | X, C = 1) * f(X) to find E{Y*(1, 0)}
mu.high <- cuhre(ndim = 2, 
  ncomp = 1, 
  integrand = ey, 
  lower = c(2.93, 7.34), 
  upper = c(17.07, 24.66),
  rel.tol = 1e-6,
  mean = mu,
  V = V,
  x.coef = x.coef,
  xy.coef = xy.coef)

# E{Y*(1, 0)}
mu.high$val
# 16.5682
