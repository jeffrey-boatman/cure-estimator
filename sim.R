# Simulation parameters ---------------------------------------------------
n.mc <- 1000             # number of Monte Carlo Iterations
n.boot <- 1000           # number of bootstrap iterations
n.cores <- 12            # number of cores to use for parallel processing
# -------------------------------------------------------------------------

options(digits = 3)
library(parallel)

# designed to run on 5 servers
host <- system2("hostname", stdout = TRUE)
hosts <- paste0(c("carbon", "cesium", "chromium", 
  "potassium", "silicon"), 
  ".ccbr.umn.edu")
n.s <- length(hosts)
j <- match(host, hosts)

source("functions.R")

# # for debugging:
# mc.it <- 1
# include.se <- TRUE
# n.boot <- 10
# include.known <- TRUE
# n <- 1000

sim <- function(mc.it, include.se, n.boot = 1000, include.known){
  out.big <- NULL # for storing results
  for(n in c(225, 500, 1000)) {
    seed <- mc.it + 1090
    set.seed(seed)

    # number of participants who are known to be compliant. 
    n_k <- 0
    if(include.known) n_k <- ceiling(n * 0.1)

    ii <- runif(n + n_k) < 0.7

    # E{Y*(1, 0)}
    true.mu <- 16.5682 

    # correlation matrix for data
    m <- matrix(c(1.000, -0.4050000, -0.3300000, -0.0800000,
                 -0.405,  1.0000000,  0.8148148,  0.7000000,
                 -0.330,  0.8148148,  1.0000000,  0.9414324,
                 -0.080,  0.7000000,  0.9414324,  1.0000000), 4, 4)

    p <- 0.20                         # Pr(C = 1)
    mu <- c(0, 10, 1.5, 16)           # mean vector
    s <- sqrt(c(1, 2, 2, 3))          # standard deviations
    V <- diag(s) %*% m %*% diag(s)    # covariance matrix

    # generate data for main sample
    ts <- t(chol(V)) 
    dd <- t(replicate(n, mu + ts %*% rnorm(4), simplify = TRUE))
    cstar <- dd[, 1]
    x <- dd[, 2]
    y <- dd[, 4]
    comp <- cstar > qnorm(1 - p)

    # generate data for known compliers
    ff <- t(replicate(1e4, mu + ts %*% rnorm(4), simplify = TRUE))
    comp_k <- ff[, 1] > qnorm(1 - 0.18)
    x_k <- ff[comp_k, 2][seq_len(n_k)]
    y_k <- ff[comp_k, 4][seq_len(n_k)]

    # combine data for entire sample.
    k <- rep(c(TRUE, FALSE), c(n_k, n))
    y <- c(y_k, y)
    x <- c(x_k, x)
    comp <- c(rep(TRUE, n_k), comp)
    comp[k] <- TRUE 

    # generate biomarker data
    gamma <- c(-9.3, -0.8, 0.7)    
    bio.lp <- c(cbind(1, as.numeric(comp), y) %*% gamma)
    biost <- rnorm(n + n_k)
    b.sd.8 <- 0.66875 # For AUC = 0.8
    b.sd.9 <- 0.44    # For AUC = 0.9
    bio.8 <- bio.lp + b.sd.8 * biost
    bio.9 <- bio.lp + b.sd.9 * biost


    # indicator for whether subject reports noncompliance without error. 
    hh <- runif(n + n_k) < 0.7

    # indicator for self-reported compliance
    src <- !k & (comp | !hh)

    # indicator for contributing to estimation of mixture dist'n
    inmix <- k | src

    # indicator for known to be compliant
    known <- as.numeric(k)

    # IPW estimator
    ipw.den.model <- glm(comp ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k)
    pi.x <- c(pnorm(cbind(1, x) %*% coef(ipw.den.model)))           
    (ipw <- weighted.mean(y, as.numeric(!k & comp) / pi.x))


    # CURE estimators --------------------------------------
    # AUC 0.8 first, then AUC 0.9
    # starting vals for EM algorithm
    starting.vals.8 <- get.starting.vals(bio = bio.8[inmix],
      x = x[inmix],
      y = y[inmix],
      comp = comp[inmix],
      known = k[inmix])
    # mixture dist'n
    fit.8 <- mixfit(bio = bio.8[inmix], 
      y = y[inmix], 
      x = x[inmix],
      known = known[inmix], 
      gamma.start = starting.vals.8$gamma.start, 
      alpha.start = starting.vals.8$alpha.start, 
      p.start = starting.vals.8$p.start)
    # numerator for CURE
    prob.compliant.8 <- rep(0, n + n_k)
    prob.compliant.8[inmix] <- fit.8$post.prob
    # denominator for CURE
    cure.den.model.8 <- suppressWarnings(glm(prob.compliant.8 ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k))	
    cure.den.weight.8 <- c(pnorm(cbind(1, x) %*% coef(cure.den.model.8)))
    (cure.8 <- weighted.mean(y, (1 - k) * prob.compliant.8 / cure.den.weight.8))

    starting.vals.9 <- get.starting.vals(bio = bio.9[inmix],
      x = x[inmix],
      y = y[inmix],
      comp = comp[inmix],
      known = k[inmix])
    fit.9 <- mixfit(bio = bio.9[inmix], 
      y = y[inmix], 
      x = x[inmix],
      known = known[inmix], 
      gamma.start = starting.vals.9$gamma.start, 
      alpha.start = starting.vals.9$alpha.start, 
      p.start = starting.vals.9$p.start)
    prob.compliant.9 <- rep(0, n + n_k)
    prob.compliant.9[inmix] <- fit.9$post.prob
    cure.den.model.9 <- suppressWarnings(glm(prob.compliant.9 ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k))	
    cure.den.weight.9 <- c(pnorm(cbind(1, x) %*% coef(cure.den.model.9)))
    (cure.9 <- weighted.mean(y, (1 - k) * prob.compliant.9 / cure.den.weight.9))

    # Cutoff IPW estimators -----------------------------------
    # AUC 0.8 first, then AUC 0.9
    # numerator
    comp.hat.8 <- ((1 - k) * prob.compliant.8) > 0.5
    # denominator
    cutoff.den.model.8 <- suppressWarnings(glm(comp.hat.8 ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k))
    cutoff.den.weight.8 <- c(pnorm(cbind(1, x) %*% coef(cutoff.den.model.8)))
    (cutoff.ipw.8 <- weighted.mean(y, comp.hat.8 / cutoff.den.weight.8))

    comp.hat.9 <- ((1 - k) * prob.compliant.9) > 0.5
    cutoff.den.model.9 <- suppressWarnings(glm(comp.hat.9 ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k))
    cutoff.den.weight.9 <- c(pnorm(cbind(1, x) %*% coef(cutoff.den.model.9)))
    (cutoff.ipw.9 <- weighted.mean(y, comp.hat.9 / cutoff.den.weight.9))

    # per protocol
    (pp <- mean(y[src]))

    # IPW based on self-reported compliance. 
    src.den.model <- glm(src ~ x, 
      family = binomial(link = "probit"),
      weights = 1 - k)
    src.den.weight <- c(pnorm(cbind(1, x) %*% coef(src.den.model)))
    (src.ipw <- weighted.mean(y, as.numeric(src)  / src.den.weight))

    # bootstrap
    if(include.se) {
      boot <- array(dim = c(n.boot, 7))
      colnames(boot) <- c("ipw", "cure.8", "cure.9", "cutoff.8", "cutoff.9",
        "pp", "src.ipw")
      for(b in seq_len(n.boot)){
        boot.seed <- 122578 + mc.it * n.boot + b
        set.seed(boot.seed)
        out <- rep(NA, 6)
        try({
          known.index <- numeric(0)
          if(n_k > 0) {
            known.index <- sample(1:n_k, replace = TRUE)
          }
          unknown.index <- sample(n_k + 1:n, replace = TRUE)
          index <- sort(c(known.index, unknown.index))
          bio.8_t <- bio.8[index]
          bio.9_t <- bio.9[index]
          k_t <- k[index]
          x_t <- x[index]
          y_t <- y[index]
          src_t <- src[index]
          inmix_t <- inmix[index]
          hh_t <- hh[index]
          comp_t <- comp[index]
          known_t <- as.numeric(k_t) 

          ipw.den.model_t <- glm(comp_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t)
          pi.x_t <- c(pnorm(cbind(1, x_t) %*% coef(ipw.den.model_t)))           
          (ipw_t <- weighted.mean(y_t, as.numeric(!k_t & comp_t) / pi.x_t))

          (out[1] <- ipw_t)

          starting.vals.8_t <- get.starting.vals(bio = bio.8_t[inmix_t],
            x = x_t[inmix_t],
            y = y_t[inmix_t],
            comp = comp_t[inmix_t],
            known = k_t[inmix_t])
          fit.8_t <- mixfit(bio = bio.8_t[inmix_t], 
            y = y_t[inmix_t], 
            x = x_t[inmix_t],
            known = known_t[inmix_t], 
            gamma.start = starting.vals.8_t$gamma.start, 
            alpha.start = starting.vals.8_t$alpha.start, 
            p.start = starting.vals.8_t$p.start)
          prob.compliant.8_t <- rep(0, n + n_k)
          prob.compliant.8_t[inmix_t] <- fit.8_t$post.prob
          cure.den.model.8_t <- suppressWarnings(glm(prob.compliant.8_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t))	
          cure.den.weight.8_t <- c(pnorm(cbind(1, x_t) %*% 
            coef(cure.den.model.8_t)))
          (cure.8_t <- weighted.mean(y_t, 
            (1 - k_t) * prob.compliant.8_t / cure.den.weight.8_t))
          out[2] <- cure.8_t

          starting.vals.9_t <- get.starting.vals(bio = bio.9_t[inmix_t],
            x = x_t[inmix_t],
            y = y_t[inmix_t],
            comp = comp_t[inmix_t],
            known = k_t[inmix_t])
          fit.9_t <- mixfit(bio = bio.9_t[inmix_t], 
            y = y_t[inmix_t], 
            x = x_t[inmix_t],
            known = known_t[inmix_t], 
            gamma.start = starting.vals.9_t$gamma.start, 
            alpha.start = starting.vals.9_t$alpha.start, 
            p.start = starting.vals.9_t$p.start)
          prob.compliant.9_t <- rep(0, n + n_k)
          prob.compliant.9_t[inmix_t] <- fit.9_t$post.prob
          cure.den.model.9_t <- suppressWarnings(glm(prob.compliant.9_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t))	
          cure.den.weight.9_t <- c(pnorm(cbind(1, x_t) %*% 
            coef(cure.den.model.9_t)))
          (cure.9_t <- weighted.mean(y_t, 
            (1 - k_t) * prob.compliant.9_t / cure.den.weight.9_t))
          out[3] <- cure.9_t

          comp.hat.8_t <- ((1 - k_t) * prob.compliant.8_t) > 0.5
          cutoff.den.model.8_t <- suppressWarnings(glm(comp.hat.8_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t))
          cutoff.den.weight.8_t <- c(pnorm(cbind(1, x_t) %*% 
            coef(cutoff.den.model.8_t)))
          (cutoff.ipw.8_t <- weighted.mean(y_t, 
            comp.hat.8_t / cutoff.den.weight.8_t))
          out[4] <- cutoff.ipw.8_t

          comp.hat.9_t <- ((1 - k_t) * prob.compliant.9_t) > 0.5
          cutoff.den.model.9_t <- suppressWarnings(glm(comp.hat.9_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t))
          cutoff.den.weight.9_t <- c(pnorm(cbind(1, x_t) %*% 
            coef(cutoff.den.model.9_t)))
          (cutoff.ipw.9_t <- weighted.mean(y_t, 
            comp.hat.9_t / cutoff.den.weight.9_t))
          out[5] <- cutoff.ipw.9_t

          (pp_t <- mean(y_t[src_t]))
          out[6] <- pp_t

          src.den.model_t <- glm(src_t ~ x_t, 
            family = binomial(link = "probit"),
            weights = 1 - k_t)
          src.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
            coef(src.den.model_t)))
          (src.ipw_t <- weighted.mean(y_t, 
            as.numeric(src_t)  / src.den.weight_t))
          out[7] <- src.ipw_t
        }, silent = TRUE)
        boot[b, ] <- out
      }

      # percentile CI
      boot.conf.ints <- t(apply(boot, 2, quantile, c(0.025, 0.975), 
        na.rm = TRUE))

      # indicator for whether CI includes true.mu
      boot.covers <- boot.conf.ints[, 1] < true.mu & 
        true.mu < boot.conf.ints[, 2]
      boot.covers <- as.numeric(boot.covers)

      boot.se <- apply(boot, 2, sd, na.rm = TRUE)
    }

    if(!include.se) {
      boot.covers <- boot.se <- rep(NA, 7)
    }

    out.little <- data.frame(n,
      mc.it,
      seed,
      ipw - true.mu,
      cure.8 - true.mu,
      cure.9 - true.mu,
      cutoff.ipw.8 - true.mu,
      cutoff.ipw.9 - true.mu,
      pp - true.mu,
      src.ipw - true.mu,
      t(boot.covers),
      t(boot.se)
    )
    out.big <- rbind(out.big, out.little) 
  }
  write.table(out.big,
    file = output.file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE)
  out.big
}


# output column names
cols <- c(
  'n', 
  'mc.it',
  'seed',
  paste0(c("ipw", "cure.8", "cure.9", "cutoff.ipw.8", "cutoff.ipw.9", "pp", "src.ipw"), ".bias"),
  paste0(c("ipw", "cure.8", "cure.9", "cutoff.ipw.8", "cutoff.ipw.9", "pp", "src.ipw"), ".covers"),
  paste0(c("ipw", "cure.8", "cure.9", "cutoff.ipw.8", "cutoff.ipw.9", "pp", "src.ipw"), ".se")
)

# create output data files
output.file <- paste0("outfiles/out_", j, ".txt")
write.table(t(cols), 
  file = output.file, 
  row.names = FALSE,   
  col.names = FALSE,
  append = FALSE, 
  quote = FALSE)


# run simulation in parallel
sims <- ((j * n.mc / n.s) - (n.mc / n.s - 1)):(j * n.mc / n.s)
out.list <- mclapply(sims,
  failwith(NULL, sim),
  mc.cores = n.cores,
  include.se = TRUE,
  n.boot = n.boot,
  include.known = TRUE)
