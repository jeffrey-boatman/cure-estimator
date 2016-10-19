# -----------------------------------------------
# code to find the true value of mu*(1, 0) via
# numerical integration for a given mean vector
# and covariance matrix used in data generation
# -----------------------------------------------

library(mvtnorm)
library(R2Cuba)

source("functions.R")

# correlation matrix for low discrimination case
rho1 <- -0.405
rho2 <- -0.37
rho3 <- -0.06
rho4 <- 0.9135802
rho5 <- 0.7
rho6 <- 0.8476062
(m.lo <- matrix(c(1,    rho1, rho2, rho3,
                  rho1, 1,    rho4, rho5,
                  rho2, rho4, 1,    rho6,
                  rho3, rho5, rho6, 1), 4, 4))

# correlation for high discrimination case
rho1 <- -0.405 
rho2 <- -0.33 
rho3 <- -0.08 
rho4 <- 0.8148148  
rho5 <- 0.7 
rho6 <- 0.9414324
(m.hi <- matrix(c(1,    rho1, rho2, rho3,
                  rho1, 1,    rho4, rho5,
                  rho2, rho4, 1,    rho6,
                  rho3, rho5, rho6, 1), 4, 4))

# p = Pr(C = 1)
# meanvec is the mean vector for data generation.
# if m is the correlation matrics, the covariance 
# matrix is diag(s) %*% m %*% diag(s)
p <- 0.18
meanvec <- c(0, 10, 1.5, 16) 
s <- sqrt(c(1, 2, 2, 3))


fy <- function(xy, mean, s, m, p) {
  # the joint density f(Y | X, C = 1) * f(X)
  # for meanvec, s, m, the order of variables is 
  # C', X, B, Y
  # meanvec: the vector of means.
  # s: the vector of standard deviations of C', X, B, Y
  # m: the correlation matrix for C', X, B, Y
  # p: Pr(C = 1)
  x <- xy[1]
  y <- xy[2]
  sigma <- diag(s) %*% m %*% diag(s)
  # 'num' = Pr(C = 1 | Y, X)
  ss <- arrange(sigma, c(1, 2, 4))
  s11 <- ss[1, 1, drop = FALSE]
  s12 <- ss[1, 2:3, drop = FALSE]
  s22 <- ss[2:3, 2:3, drop = FALSE]
  s21 <- t(s12)
  eta <- mean[1] + as.numeric(s12 %*% solve(s22) %*% 
    cbind(c(x, y) - mean[c(2, 4)]))
  sigsq <- as.numeric(s11 - s12 %*% solve(s22) %*% s21)
  psi <- qnorm(1 - p)
  num <- pnorm((eta - psi) / sqrt(sigsq))
  # 'den' = Pr(C = 1 | X)
  ss <- arrange(sigma, c(1, 2))
  s11 <- ss[1, 1, drop = FALSE]
  s12 <- ss[1, 2, drop = FALSE]
  s22 <- ss[2, 2, drop = FALSE]
  s21 <- t(s12)
  eta <- mean[1] + as.numeric(s12 %*% solve(s22) %*% 
    cbind(x - mean[2]))
  sigsq <- as.numeric(s11 - s12 %*% solve(s22) %*% s21)
  psi <- qnorm(1 - p)
  den <- pnorm((eta - psi) / sqrt(sigsq))
  num / den * dmvnorm(xy, mean[c(2, 4)], arrange(sigma, c(2, 4)))
}

ey <- function(xy, mean, s, m, p) {
  # expectation of y*(1, 0). see function 'fy' for a 
  # description of the arguments.
  y <- xy[2]
  y * fy(xy, mean, s, m, p)
}

# --------------------------------------------
# find E{Y*(1,0 )} by integrating ey.
# for limits integration, use mean +/- 5 sd.
# does not work with +/- Inf
# --------------------------------------------
meanvec + cbind(-5 * s, 5 * s)

# -----------------------------------------
# E{Y*(1, 0)} for low dsicrmination case. 
# -----------------------------------------
mu.lo <- cuhre(ndim = 2, 
  ncomp = 1, 
  integrand = ey, 
  lower = c(2.93, 7.34), 
  upper = c(17.07, 24.66),
  rel.tol = 1e-6,
  mean = meanvec,
  s = s,
  m = m.lo,
  p = p)

(mu.lo <- mu.lo$val)
# 16.65094


# -----------------------------------------
# E{Y*(1, 0)} for high dsicrmination case. 
# -----------------------------------------
mu.hi <- cuhre(ndim = 2, 
  ncomp = 1, 
  integrand = ey, 
  lower = c(2.93, 7.34), 
  upper = c(17.07, 24.66),
  rel.tol = 1e-6,
  mean = meanvec,
  s = s,
  m = m.hi,
  p = p)
(mu.hi <- mu.hi$val)
# 16.5927
