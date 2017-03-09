mixfit <- function (bio, y, x, known, tol = 1e-06, gamma.start, alpha.start, 
  p.start) 
{		
  # Estimates the parameters of the mixture distribution.
  # Arguments:
  #  - bio: vector with the biomarker
  #  - y: vector with the outcome
  #  - x: vector with the confounder
  #  - known: vector with indicator for whether compliance is known
  #  - tol: tolerance for changes in negative log likelihood.
  #  - gamma.start: starting coefficients for f(B | Y, C)
  #  - alpha.start: starting coefficients for Pr(C = 1 | X, Y)
  #  - p.start: starting estimates for Pr(C = 1 | B = b, Y = y, X = x, D = 1)
  # Returns a list with the following elements:
  #  - gamma.hat: estimated gamma coefficients
  #  - alpha.hat: estimated alpha coefficients
  #  - prob.complaiance: estimated Pr(C = 1 | D = 1)
  #  - post.prob: Estimates of Pr(C = 1, B = b, X = x, Y = y, D = 1)
  #  - num.its: number of iterations until convergence
  #  - neg.log.lik: negative log likelihood

  known <- as.numeric(known)
  p.hat <- p.start
  alpha <- alpha.start
  gamma <- gamma.start  

  # function to evaluate the density
  dens <- function(gamma){
    temp <- matrix(nrow = length(bio), ncol = 3)
    (gamma.c <- gamma[1:2])
    (gamma.nc <- gamma[3:4])
    (b.c.sd <- gamma[5])
    (b.nc.sd <- gamma[5])

    temp[, 1] <- temp[, 2] <- dnorm(bio, c(cbind(1, y) %*% gamma.c), b.c.sd)
    temp[, 3] <- dnorm(bio, c(cbind(1, y) %*% gamma.nc), b.nc.sd)
    temp
  }

  # function to evaluate the incomplete negative log likelihood.
  incomp.nll <- function(gamma, p.hat){
    # p.hat <- pnorm(c(cbind(1, x, x2) %*% gamma))
    k <- cbind(known, 
      (1 - known) * p.hat, 
      (1 - known) * (1 - p.hat))
    -sum(log(apply(k * dens(gamma), 1, sum)))
  }

  k <- cbind(known, 
    (1 - known) * p.hat, 
    (1 - known) * (1 - p.hat))
  nll <- -sum(log(apply(k * dens(gamma), 1, sum)))

  # first E-step
  dens_ <- cbind(p.hat, 1- p.hat) * (dens(gamma)[, -1])
  w <- dens_ / apply(dens_, 1, sum)

  nll.diff <- 1
  iter <- 0		
  
  nll.total <- nll

  # weighted sum of squares function
  wss <- function(pars, yy, w1, w2, ...) {
    ll <- list(...)
    xx <- cbind(1, do.call(cbind, ll))
    cc <- w1 * (yy - c(xx %*% pars[-2])) ^ 2
    nc <- w2 * (yy - c(xx %*% pars[-1])) ^ 2
    sum(cc, nc)
  }

  while(nll.diff > tol) {
    gg <- numeric(length(gamma))
    aa <- numeric(length(alpha))
    comp_model <- suppressWarnings(
      glm(w[, 1] ~ x + y, weights = (1 - known),
      family = binomial(link = "probit")))
    aa <- coef(comp_model)
  
    bio.lm <- optim(gamma[c(1, 3, 2)], 
      fn = wss, 
      yy = bio,
      w1 = known + (1 - known) * w[, 1],
      w2 = (1 - known) * w[, 2],
      y = y)

    gg[1:2] <- bio.lm$par[-2]
    gg[3:4] <- bio.lm$par[-1]
    gg[5:6] <- weighted.sd(bio, 
      mean.1 = c(cbind(1, y) %*% bio.lm$par[-2]),
      mean.2 = c(cbind(1, y) %*% bio.lm$par[-1]),
      weight.1 = known + (1 - known) * w[, 1],
      weight.2 = (1 - known) * w[, 2])

    gamma <- gg
    alpha <- unname(aa)

    p.hat <- pnorm(c(cbind(1, x, y) %*% alpha))
  
    k <- cbind(known,
      (1 - known) * p.hat,
      (1 - known) * (1 - p.hat))
  
    # E-step
    dens_ <- cbind(p.hat, 1- p.hat) * (dens(gamma)[, -1])
    w <- dens_ / apply(dens_, 1, sum)

    nll.update <- incomp.nll(gamma, p.hat)
    nll.diff <- nll - nll.update

    nll <- nll.update
    iter <- iter + 1
    nll.total <- c(nll.total, nll)
  }

  dens_ <- cbind(p.hat, 1 - p.hat) * dens(gamma)[, -1]
  post.prob <- dens_[, 1]/apply(dens_, 1, sum)

  list(gamma.hat = gamma, 
    alpha.hat = alpha,
    prob.compliance = mean(p.hat[!known]), 
    post.prob = post.prob,
    num.its = length(nll.total),
    neg.log.lik = nll.total[length(nll.total)])
}

failwith <- function(default = NULL, f, quiet = TRUE) {
  # borrowed from dplyr package. 
  function(...) {
    out <- default
    try(out <- f(...), silent = quiet)
    out
  }
}

weighted.sd <- function(x, mean.1, mean.2, weight.1, weight.2) {
  # returns weighted sd
  sum.w <- sum(weight.1) + sum(weight.2)
  sqrt(sum(weight.1 * (x - mean.1) ^ 2 + weight.2 * (x - mean.2) ^ 2) / sum.w)
}

xyprobit <- function(V, mu, pp) {
  # Returns the coefficients of probit regression of C on X, Y.
  # V is the variance matrix of C', x, ..., Xn, and Y, in that order.
  # mu is the mean vector of C', x, ...Xn, Y, in that order.
  # pp = Pr(C = 1)
  (s11 <- V[1, 1, drop = FALSE])
  (s12 <- V[1, -1, drop = FALSE])
  (s22 <- V[-1, -1])
  # mu[1] = 0 but including it for sake of coding clarity
  int <- mu[1] - s12 %*% solve(s22) %*% cbind(mu[-1])
  slope <- s12 %*% solve(s22)
  (variance <- s11 - s12 %*% solve(s22) %*% t(s12))
  offset <- -qnorm(1 - pp)
  c(int + offset, slope) / sqrt(variance)
}

xprobit <- function(V, mu, pp) {
  # Returns the coefficients of probit regression of C on X.
  # V is the variance matrix of C', x, ..., Xn, and Y, in that order.
  # mu is the mean vector of C', x, ..., Xn, Y, in that order.
  # pp = Pr(C = 1)
  nn <- nrow(V)
  V <- V[-nn, -nn]
  mu <- mu[-nn]
  (s11 <- V[1, 1, drop = FALSE])
  (s12 <- V[1, -1, drop = FALSE])
  (s22 <- V[-1, -1])
  # mu[1] = 0 but including it for sake of coding clarity
  int <- mu[1] - s12 %*% solve(s22) %*% cbind(mu[-1])
  slope <- s12 %*% solve(s22)
  (variance <- s11 - s12 %*% solve(s22) %*% t(s12))
  offset <- -qnorm(1 - pp)
  c(int + offset, slope) / sqrt(variance)
}

get.starting.vals <- function(bio, x, y, comp, known) {
  # gets starting vals for the EM algorithm assuming we're estimating
  # the numerator of the weights as Pr(C = 1 | B, X, Y) * f(B | Y, C = 1) / sum(...)
  bio.lm <- lm(bio ~ comp + y)
  bio.coefs <- coef(bio.lm)
  bio.rmse <- sqrt(mean(residuals(bio.lm) ^ 2))
  gamma.start <- unname(c(sum(bio.coefs[1:2]), bio.coefs[c(3, 1, 3)], rep(bio.rmse, 2)))
  alpha.glm <- glm(comp ~ x + y, 
    family = binomial(link = "probit"),
    weights = 1 - known)
  alpha.start <- unname(coef(alpha.glm))
  p.start <- pnorm(c(cbind(1, x, y) %*% alpha.start))
  list(gamma.start = gamma.start,
    alpha.start = alpha.start,
    p.start = p.start)
}

fy <- function(xy, mean, V, x.coef, xy.coef) {
  # The joint density f(Y | X1, X2, C = 1) * f(X1, X2).
  # xy: the vector(X1, X2, Y)
  # V: the covariance matrix for C', X1, X2, Y, in that order.
  # x.coef: the coefficients for the probit model Pr(C = 1 | X1, X2).
  # xy.coef: the coefficients for the probit model Pr(C = 1 | X1, X2, Y).
  # --------------------------------------------------------------------
  # 'num' = Pr(C = 1 | Y, X)
  num <- pnorm(c(rbind(c(1, xy)) %*% xy.coef))
  # 'den' = Pr(C = 1 | X)
  den <- pnorm(c(rbind(c(1, xy[-length(xy)])) %*% x.coef))
  num / den * mvtnorm::dmvnorm(xy, mean[-1], V[-1, -1])
}

ey <- function(xy, mean, V, x.coef, xy.coef) {
  # expectation of y*(1, 0). see function 'fy' for a 
  # description of the arguments.
  y <- xy[length(xy)]
  y * fy(xy, mean, V, x.coef, xy.coef)
}

colSd <- function(x) apply(x, 2, sd, na.rm = TRUE)

summarize <- function(df) {
  # summary function used in creating output table
  keepers <- c("pp", "src.ipw", "ipw", "cutoff.ipw.8",
    "cure.8", "cutoff.ipw.9", "cure.9")
  bias <- colMeans(df[paste0(keepers, ".bias")], 
    na.rm = TRUE)
  mcsd <- colSd(df[paste0(keepers, ".bias")])
  mean.se <- colMeans(df[paste0(keepers, ".se")], na.rm = TRUE)
  cover.prob <- colMeans(df[paste0(keepers, ".covers")],
    na.rm = TRUE)
  mse <- colMeans(df[paste0(keepers, ".sqe")], 
    na.rm = TRUE)
  out <- cbind(bias, mcsd, mean.se, cover.prob, mse)
  out
}

`%!in%` <- Negate(`%in%`)
