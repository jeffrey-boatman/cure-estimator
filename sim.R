# Simulation parameters --------------------------------------
n.mc <- 1000             # number of Monte Carlo Iterations
n.boot <- 1000           # number of bootstrap iterations
n.cores <- 4             # number of cores to use for parallel processing
output.file <- "out.txt" # where to write results

library(parallel)

source("functions.R")

sim <- function(mc.it, n.boot = 1000){
  # main simulation function. 
  # mc.it: the monte carlo iteration. 
  # n.boot: the number of bootstrap iterations to use.
  # Writes to the file specified by output.file.
  # The output includes:
  # n: sample size
  # dis: low or high for low discrmination or high discrimination
  # mc.it: the monte carlo iteration number
  # seed: the random number seed
  # the bias of the CURE, cutoff IPW, per protocol, and self-report IPW, 
  #  in that order
  # indicators for whether the 95% CI includes the true value for CURE, cutoff IPW,
  #  per protocol, and self-report IPW, in that order.
  # the estimated SE of the estimator for CURE, cutoff IPW, per protocol, and self-
  #   report IPW, in that order.  

  out.big <- NULL # for storing results
  for(n in c(225, 500, 1000)) {
    for(dis in c('low', 'high')) {

        set.seed(mc.it + 1090)

        # number of participants who are known to be compliant. 
        n_k <- ceiling(n * 0.1)

        if(dis == 'low') {
          # for low discrimination case
          true.mu <- 16.65094
          rho1 <- -0.405
          rho2 <- -0.37
          rho3 <- -0.06
          rho4 <- 0.9135802
          rho5 <- 0.7
          rho6 <- 0.8476062
        }
        if(dis == 'high') {
          # for high discrimination case
          true.mu <- 16.5927   
          rho1 <- -0.405 
          rho2 <- -0.33 
          rho3 <- -0.08 
          rho4 <- 0.8148148  
          rho5 <- 0.7 
          rho6 <- 0.9414324
        }

        # correlation matrix for data generation. 
        m <- matrix(c(1,    rho1, rho2, rho3,
                      rho1, 1,    rho4, rho5,
                      rho2, rho4, 1,    rho6,
                      rho3, rho5, rho6, 1), 4, 4)

 
        # Pr(C = 1) = 0.18
        p <- 0.18 


        # indicator for whether subject reports noncompliance without error. 
        hh <- runif(n + n_k) < 0.7

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
        x_k <- ff[comp_k, 2][seq_len(n_k)]
        bio_k <- ff[comp_k, 3][seq_len(n_k)]
        y_k <- ff[comp_k, 4][seq_len(n_k)]

        # indicator for whether compliance status is known
        k <- rep(c(TRUE, FALSE), c(n_k, n))

        # combine data for entire sample.
        k <- rep(c(TRUE, FALSE), c(n_k, n))
        bio <- c(bio_k, bio)
        y <- c(y_k, y)
        x <- c(x_k, x)
        comp <- c(rep(TRUE, n_k), comp)
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

        # estimating the mixturedistribution.
        fit <- mixfit(bio = bio[inmix], 
          y = y[inmix], 
          known = k[inmix],
          start = theta.start, 
          p_start = mean(comp[src]))

        # numerator for the CURE estimator. 
        prob.compliant <- rep(0, n + n_k)
        prob.compliant[inmix] <- fit$post.prob

        # denominator for the CURE estimator
        cure.den.model <- suppressWarnings(glm(prob.compliant[!k] ~ x[!k], 
          family = binomial(link = "probit")))	
        cure.den.weight <- c(pnorm(cbind(1, x) %*% coef(cure.den.model)))

        # CURE
        (cure <- weighted.mean(y, prob.compliant * src / cure.den.weight))

        # numerator for the cutoff IPW estimator. 
        comp.hat <- (prob.compliant * src) > 0.5

        # denominator for the cutoff IPW estimator.
        cutoff.den.model <- suppressWarnings(glm(comp.hat[!k] ~ x[!k], 
          family = binomial(link = "probit")))
        cutoff.den.weight <- c(pnorm(cbind(1, x) %*% coef(cutoff.den.model)))

        # cutoff IPW
        (cutoff.ipw <- weighted.mean(y, (comp.hat)/cutoff.den.weight))

        # per protocol
        (pp <- mean(y[src]))

        # denominator for IPW based on self-reported compliance. 
        src.den.model <- glm(src[!k] ~ x[!k], family = binomial(link = "probit"))
        src.den.weight <- c(pnorm(cbind(1, x) %*% coef(src.den.model)))
        (src.ipw <- weighted.mean(y, src  / src.den.weight))

        # bootstrap
        boot <- array(dim = c(n.boot, 4))
        set.seed(122578 + mc.it * 1000)
        for(b in seq_len(n.boot)){
          out <- rep(NA, 4)
          try({
            known.index <- sample(1:n_k, replace = TRUE)
            unknown.index <- sample(n_k + 1:n, replace = TRUE)
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

            prob.compliant_t <- rep(0, n + n_k)
            prob.compliant_t[inmix_t] <- fit_t$post.prob
            cure.den.model_t <- suppressWarnings(glm(prob.compliant_t[!k_t] ~ x[!k_t], 
              family = binomial(link = "probit")))	
            cure.den.weight_t <- c(pnorm(cbind(1, x_t) %*% coef(cure.den.model_t)))
            (cure_t <- weighted.mean(y_t, prob.compliant_t * 
              src_t / cure.den.weight_t))
            out[1] <- cure_t

            # 5th estimator. Same as 5, with cutoff for compliance. 
            comp.hat_t <- (prob.compliant_t * src_t) > 0.5
            cutoff.den.model_t <- suppressWarnings(glm(comp.hat_t[!k_t] ~ x_t[!k_t], 
              family = binomial(link = "probit")))
            cutoff.den.weight_t <- c(pnorm(cbind(1, x_t) %*% 
              coef(cutoff.den.model_t)))
            (cutoff.ipw_t <- weighted.mean(y_t, comp.hat_t / cutoff.den.weight_t))
            out[2] <- cutoff.ipw_t

            # mean of self-reported compliers
            # (naive1 <- mean(y[!k & comp | !hh]))
            (pp_t <- mean(y_t[src_t]))
            out[3] <- pp_t

            # IPW w/ self-reported compliers. 
            src.den.model_t <- glm(src_t[!k_t] ~ x_t[!k_t], 
              family = binomial(link = "probit"))
            src.den.weight_t <- c(pnorm(cbind(1, x_t) %*% coef(src.den.model_t)))
            (src.ipw_t <- weighted.mean(y_t, src_t  / src.den.weight_t))
            out[4] <- src.ipw_t

          }, silent = TRUE)
          boot[b, ] <- out
        }

        # estimated SEs
        boot.se <- apply(boot, 2, sd, na.rm = TRUE)

        # percentile CI
        boot.conf.ints <- t(apply(boot, 2, quantile, c(0.025, 0.975), na.rm = TRUE))

        # indicator for whether CI includes true.mu
        boot.covers <- boot.conf.ints[, 1] < true.mu & 
          true.mu < boot.conf.ints[, 2]
        boot.covers <- as.numeric(boot.covers)

        out.little <- c(n,
          dis, 
          mc.it,
          mc.it + 1090,
          cure - true.mu,
          cutoff.ipw - true.mu,
          pp - true.mu,
          src.ipw - true.mu,
          boot.covers,
          boot.se
        )
        out.big <- rbind(out.big, out.little)
    }
  }
  out.big
}

# output column names
cols <- c(
  'n', 
  'dis',
  'mc.it',
  'seed',
  paste0(c("cure", "cutoff.ipw", "pp", "src.ipw"), ".bias"),
  paste0(c("cure", "cutoff.ipw", "pp", "src.ipw"), ".covers"),
  paste0(c("cure", "cutoff.ipw", "pp", "src.ipw"), ".se")
)

# run simulation in parallel
out.list <- mclapply(1:n.mc, 
  failwith(NULL, sim),
  mc.cores = n.cores,
  n.boot = n.boot)

out <- do.call(rbind, out.list)
colnames(out) <- cols

write.table(out,
  file = output.file,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  append = TRUE)