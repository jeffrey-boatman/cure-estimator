# ---------------------------------------------------------
# Example analysis using simulated data set. 
# Compares ITT, per protocol, IPW based on self-report, 
# IPW using "cutoff" estimator, and the CURE estimator.
# Creates a table similar as in the application in the
# manuscript.
# ---------------------------------------------------------


source("functions.R")

example.data <- read.table("example_data.txt", header = TRUE, as.is = TRUE)

# x is the confounder
# bio is the biomarker
# y is the outcome, src is the indicator for self-reported compliance
# src is the indicator for self-reported compliance. Note that 
#  participants with known = 1 have src = 0
# known is the indicator for whether ppnt's compliance status is known
# treatment is the treatment group. For simplicity, all particpants in
#  group B are assumed to be compliant.
str(example.data)

group.a <- subset(example.data, treatment == "A")
group.b <- subset(example.data, treatment == "B")

(y.a2.mean <- mean(group.b$y))

# get estimated mean for group A.
x <- group.a$x
bio <- group.a$bio
y <- group.a$y
src <- group.a$src
known <- k <- group.a$known

n <- sum(1 - k)
n_k <- sum(k)

# indicator for whether subject contributes to the 
# mixture distribution estimation
inmix <- k | src

# guess on who compliant participants are
bio.cut <- quantile(bio, 0.2)
comp.guess <- bio <= bio.cut

# get starting values for em algorithm
starting.vals <- get.starting.vals(bio = bio[inmix],
  x = x[inmix],
  y = y[inmix],
  comp = comp.guess[inmix],
  known = k[inmix])

# CURE estimator
fit <- mixfit(bio = bio[inmix], 
  y = y[inmix], 
  x = x[inmix],
  known = known[inmix], 
  gamma.start = starting.vals$gamma.start, 
  alpha.start = starting.vals$alpha.start, 
  p.start = starting.vals$p.start)
# numerator for CURE
prob.compliant <- rep(0, n + n_k)
prob.compliant[inmix] <- fit$post.prob
# denominator for CURE
cure.den.model <- suppressWarnings(glm(prob.compliant ~ x, 
  family = binomial(link = "probit"),
  weights = 1 - k))	
cure.den.weight <- c(pnorm(cbind(1, x) %*% coef(cure.den.model)))
(cure <- y.a2.mean - 
  weighted.mean(y, (1 - k) * prob.compliant / cure.den.weight))

# cutoff IPW
comp.hat <- ((1 - k) * prob.compliant) > 0.5
# denominator
cutoff.den.model <- suppressWarnings(glm(comp.hat ~ x, 
  family = binomial(link = "probit"),
  weights = 1 - k))
cutoff.den.weight <- c(pnorm(cbind(1, x) %*% coef(cutoff.den.model)))
(cutoff.ipw <- y.a2.mean - weighted.mean(y, comp.hat / cutoff.den.weight))

# per protocol
(pp <- y.a2.mean - mean(y[src]))

# IPW based on self-reported compliance. 
src.den.model <- glm(src ~ x, 
  family = binomial(link = "probit"),
  weights = 1 - k)
src.den.weight <- c(pnorm(cbind(1, x) %*% coef(src.den.model)))
(src.ipw <- y.a2.mean - weighted.mean(y, as.numeric(src)  / src.den.weight))

# ITT
itt <- (y.a2.mean - mean(y[!k]))

n.boot <- 1000
boot <- array(dim = c(n.boot, 5))
set.seed(322)
for(b in seq_len(n.boot)){
  out <- rep(NA, 5)
  try({
    known.index <- numeric(0)
    if(n_k > 0) {
      known.index <- sample(1:n_k, replace = TRUE)
    }
    unknown.index <- sample(n_k + 1:n, replace = TRUE)
    index <- sort(c(known.index, unknown.index))
    bio_t <- bio[index]
    k_t <- k[index]
    x_t <- x[index]
    y_t <- y[index]
    src_t <- src[index]
    inmix_t <- inmix[index]
    comp.guess_t <- comp.guess[index]
    known_t <- known[index]
    y.a2.mean_t <- mean(sample(group.b$y, replace = TRUE))

    starting.vals_t <- get.starting.vals(bio = bio_t[inmix_t],
      x = x_t[inmix_t],
      y = y_t[inmix_t],
      comp = comp.guess_t[inmix_t],
      known = k_t[inmix_t])
    fit_t <- mixfit(bio = bio_t[inmix_t], 
      y = y_t[inmix_t], 
      x = x_t[inmix_t],
      known = known_t[inmix_t], 
      gamma.start = starting.vals_t$gamma.start, 
      alpha.start = starting.vals_t$alpha.start, 
      p.start = starting.vals_t$p.start)
    prob.compliant_t <- rep(0, n + n_k)
    prob.compliant_t[inmix_t] <- fit_t$post.prob
    cure.den.model_t <- suppressWarnings(glm(prob.compliant_t ~ x_t, 
      family = binomial(link = "probit"),
      weights = 1 - k_t))	
    cure.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
      coef(cure.den.model_t)))
    (cure_t <- y.a2.mean_t - weighted.mean(y_t, 
      (1 - k_t) * prob.compliant_t / cure.den.weight_t))

    comp.hat_t <- ((1 - k_t) * prob.compliant_t) > 0.5
    cutoff.den.model_t <- suppressWarnings(glm(comp.hat_t ~ x_t, 
      family = binomial(link = "probit"),
      weights = 1 - k_t))
    cutoff.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
      coef(cutoff.den.model_t)))
    (cutoff.ipw_t <- y.a2.mean_t - weighted.mean(y_t, 
      comp.hat_t / cutoff.den.weight_t))

    (pp_t <- y.a2.mean_t - mean(y_t[src_t]))

    src.den.model_t <- glm(src_t ~ x_t, 
      family = binomial(link = "probit"),
      weights = 1 - k_t)
    src.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
      coef(src.den.model_t)))
    (src.ipw_t <- y.a2.mean_t - weighted.mean(y_t, 
      as.numeric(src_t)  / src.den.weight_t))

    itt_t <- (y.a2.mean_t - mean(y_t[!k_t]))

    out[1] <- itt_t
    out[2] <- pp_t
    out[3] <- src.ipw_t
    out[4] <- cutoff.ipw_t
    out[5] <- cure_t

  }, silent = TRUE)
  boot[b, ] <- out
}

# percentile CI
boot.conf.ints <- t(apply(boot, 2, quantile, c(0.025, 0.975), 
  na.rm = TRUE))

# estimated SEs
boot.se <- apply(boot, 2, sd, na.rm = TRUE)

# create output table
summary.table <- cbind(y.a2.mean,
  y.a2.mean - c(itt, pp, src.ipw, cutoff.ipw, cure),
  c(itt, pp, src.ipw, cutoff.ipw, cure),
  boot.se,
  boot.conf.ints)

estimators <- c("ITT", "Per Protocol", "Self-Report IPW", "Cutoff IPW", "CURE")
dimnames(summary.table) <- list("Estimator" = estimators)
colnames(summary.table) <- c("mu(A, 0)", "mu(B, 0)", "mu(A, 0) - mu(B, 0)",
  "SE", "2.5%", "97.5%")

round(summary.table, 2)	
