# ---------------------------------------------------------
# Example analysis using simulated data set. 
# Compares ITT, per protocol, IPW based on self-report, 
# IPW using "cutoff" estimator, and the CURE estimator.
# Creates a table similar as in the application in the
# manuscript.
# ---------------------------------------------------------


source("functions.R")

n <- 1000 # number of observations in each treatment group.

set.seed(333)

# generate data for treatment group 2. --------------------------------
# For this group, mu(2, 0) = 20. All participants are assumed to be 
# compliant for simplicity.
y.a2 <- rnorm(n, 20, 3)
(y.a2.mean <- mean(y.a2))

# generate data from treatment group 2. -------------------------------
# For this group, mu(1, 0) = 16.5927. Some participants are 
# noncompliant.

# mu(2, 0) - mu(1, 0) = 3.4073 is the target of inference. 

# Pr(C = 1) = 0.18
p <- 0.18 

# correlation matrix for data generation. 
rho1 <- -0.405 
rho2 <- -0.33 
rho3 <- -0.08 
rho4 <- 0.8148148  
rho5 <- 0.7 
rho6 <- 0.9414324
m <- matrix(c(1,    rho1, rho2, rho3,
              rho1, 1,    rho4, rho5,
              rho2, rho4, 1,    rho6,
              rho3, rho5, rho6, 1), 4, 4)

n.k <- 100 # number of participants known to be compliant.
nn <- n + n.k

# indicator for whether subject reports noncompliance without error. 
hh <- runif(nn) < 0.7

# generate data for main sample
meanvec <- c(0, 10, 1.5, 16)
s <- sqrt(c(1, 2, 2, 3))
sigma <- diag(s) %*% m %*% diag(s)
ts <- t(chol(sigma)) 
dd <- t(replicate(n, meanvec + ts %*% rnorm(4), simplify = TRUE))
cstar <- dd[, 1]
x <- dd[, 2]
bio <- dd[, 3]
y <- dd[, 4]
comp <- cstar > qnorm(1 - p)

# generate data for known compliers
ff <- t(replicate(1e4, meanvec + ts %*% rnorm(4), simplify = TRUE))
comp_k <- ff[, 1] > qnorm(1 - 0.18)
x_k <- ff[comp_k, 2][seq_len(n.k)]
bio_k <- ff[comp_k, 3][seq_len(n.k)]
y_k <- ff[comp_k, 4][seq_len(n.k)]

# indicator for whether compliance status is known
k <- rep(c(TRUE, FALSE), c(n.k, n))

# combine data for entire sample.
k <- rep(c(TRUE, FALSE), c(n.k, n))
bio <- c(bio_k, bio)
y <- c(y_k, y)
x <- c(x_k, x)
comp <- c(rep(TRUE, n.k), comp)
comp[k] <- TRUE 

# indicator for self-reported compliance.
src <- !k & (comp | !hh)

# indicator for whether subject contributes to the 
# mixture distribution estimation
inmix <- k | src

# starting values for estimating the mixture distribution.
ymeanstart <- c(mean(y[comp & !k]), mean(y[!comp & inmix]))
ysdstart <- weighted.sd(y[!k & inmix], 
  ymeanstart[1], 
  ymeanstart[2],
  as.numeric(comp[!k & inmix]),
  as.numeric(1 - comp[!k & inmix]))              
cmod <- lm(bio ~ y, weights = as.numeric(comp))
ncmod <- lm(bio ~ y, weights = as.numeric(!comp & inmix))
rr <- weighted.sd(bio, fitted(cmod), fitted(ncmod), 
  as.numeric(comp),
  as.numeric(!comp & inmix))
theta.start <- c(
  ymeanstart,
  rep(ysdstart, 2),
  coef(cmod)[1],
  coef(ncmod)[1],
  coef(cmod)[2],
  coef(ncmod)[2],
  rep(rr, 2))

# estimating the mixture distribution.
fit <- mixfit(bio = bio[inmix], 
  y = y[inmix], 
  known = k[inmix],
  start = theta.start, 
  p_start = mean(comp[src]))

# numerator for CURE. 
prob.compliant <- rep(0, nn)
prob.compliant[inmix] <- fit$post.prob

# denominator for CURE
cure.den.model <- suppressWarnings(glm(prob.compliant[!k] ~ x[!k], 
  family = binomial(link = "probit")))	
cure.den.weight <- c(pnorm(cbind(1, x) %*% coef(cure.den.model)))

# CURE
(cure <- y.a2.mean - weighted.mean(y, prob.compliant * src / cure.den.weight))

# numerator for the cutoff IPW. 
comp.hat <- (prob.compliant * src) > 0.5

# denominator for the cutoff IPW.
cutoff.den.model <- suppressWarnings(glm(comp.hat[!k] ~ x[!k], 
  family = binomial(link = "probit")))
cutoff.den.weight <- c(pnorm(cbind(1, x) %*% coef(cutoff.den.model)))

# cutoff IPW
(cutoff.ipw <- y.a2.mean - weighted.mean(y, (comp.hat)/cutoff.den.weight))

# per protocol
(pp <- y.a2.mean - mean(y[src]))

# denominator for IPW based on self-reported compliance. 
src.den.model <- glm(src[!k] ~ x[!k], family = binomial(link = "probit"))
src.den.weight <- c(pnorm(cbind(1, x) %*% coef(src.den.model)))

# IPW based on self-reported compliant. 
(src.ipw <- y.a2.mean - weighted.mean(y, src  / src.den.weight))

# ITT
(itt <- y.a2.mean - mean(y[!k]))

# bootstrap
n.boot <- 1000
boot <- array(dim = c(n.boot, 5))
for(b in seq_len(n.boot)){
  out <- rep(NA, 5)
  try({
    y.a2.mean_t <- mean(sample(y.a2, replace = TRUE))
    known.index <- sample(1:n.k, replace = TRUE)
    unknown.index <- sample(n.k + 1:n, replace = TRUE)
    index <- sort(c(known.index, unknown.index))
    bio_t <- bio[index]
    k_t <- k[index]
    x_t <- x[index]
    y_t <- y[index]
    src_t <- src[index]
    inmix_t <- inmix[index]
    hh_t <- hh[index]
    comp_t <- comp[index]

    fit_t <- mixfit(bio = bio_t[inmix_t], 
      y = y_t[inmix_t], 
      known = k_t[inmix_t], 
      start = c(fit$theta),
      p_start = fit$prob)	

    prob.compliant_t <- rep(0, nn)
    prob.compliant_t[inmix_t] <- fit_t$post.prob
    cure.den.model_t <- suppressWarnings(glm(prob.compliant_t[!k_t] ~ x[!k_t], 
      family = binomial(link = "probit")))	
    cure.den.weight_t <- c(pnorm(cbind(1, x_t) %*% coef(cure.den.model_t)))
    (cure_t <- y.a2.mean_t - weighted.mean(y_t, prob.compliant_t * 
      src_t / cure.den.weight_t))


    # 5th estimator. Same as 5, with cutoff for compliance. 
    comp.hat_t <- (prob.compliant_t * src_t) > 0.5
    cutoff.den.model_t <- suppressWarnings(glm(comp.hat_t[!k_t] ~ x_t[!k_t], 
      family = binomial(link = "probit")))
    cutoff.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
      coef(cutoff.den.model_t)))
    (cutoff.ipw_t <- y.a2.mean_t - weighted.mean(y_t, comp.hat_t / 
      cutoff.den.weight_t))

    (pp_t <- y.a2.mean_t - mean(y_t[src_t]))

    src.den.model_t <- glm(src_t[!k_t] ~ x_t[!k_t], 
      family = binomial(link = "probit"))
    src.den.weight_t <- c(pnorm(cbind(1, x_t) %*% coef(src.den.model_t)))
    (src.ipw_t <- y.a2.mean_t - weighted.mean(y_t, src_t  / src.den.weight_t))

    # ITT
    (itt_t <- y.a2.mean_t - mean(y_t[!k_t]))

    out[1] <- itt_t
    out[2] <- pp_t
    out[3] <- src.ipw_t
    out[4] <- cutoff.ipw_t
    out[5] <- cure_t

  }, silent = TRUE)
  boot[b, ] <- out
}

# estimated SEs
boot.se <- apply(boot, 2, sd, na.rm = TRUE)

# percentile CI
boot.conf.ints <- t(apply(boot, 2, quantile, c(0.025, 0.975), na.rm = TRUE))

# create output table
summary.table <- cbind(y.a2.mean,
  y.a2.mean - c(itt, pp, src.ipw, cutoff.ipw, cure),
  c(itt, pp, src.ipw, cutoff.ipw, cure),
  boot.se,
  boot.conf.ints)

estimators <- c("ITT", "Per Protocol", "Self-Report IPW", "Cutoff IPW", "CURE")
dimnames(summary.table) <- list("Estimator" = estimators)
colnames(summary.table) <- c("mu(2, 0)", "mu(1, 0)", "mu(2, 0) - mu(1, 0)",
  "SE", "2.5%", "97.5%")

round(summary.table, 2)
