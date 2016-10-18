mixfit <- function (bio, y, known, tol = 1e-06, start = NULL, p_start = NULL) 
{		
  # Estimates the parameters of the mixture distribution.
  # Arguments:
  #  - bio: vector with the biomarker
  #  - y: vector with the outcome
  #  - tol: tolerance for changes in negative log likelihood.
  #  - start: starting values for parameters for EM algorithm.
  #  - p_start: starting value for Pr(C = 1 | B, Y, D = 1)
  # Returns a list with the following elements:
  #  - theta.hat: the estimated parameters of the mixture distribution
  #  - prob.compliance: estimated Pr(C = 1 | D = 1)
  #  - post.prob: vector of estimated Pr(C = 1 | B, Y, D = 1). Note that
  #    this returns estimates even for observations who are  
  #    known to be compliant. 

  known <- as.numeric(known)
  p.hat <- p_start
  theta <- start  

  # function to evaluate the density
  dens <- function(theta){
    temp <- matrix(nrow = length(bio), ncol = 3)
    ymean <- theta[1:2]
    ysd <- theta[3:4]
    bint <- theta[5:6]
    bslope <- theta[7:8]
    bse <- theta[9:10]
    temp[, 1] <- dnorm(bio, mean = bint[1] + bslope[1] * y, sd = bse[1])
    temp[, 2] <- 
      dnorm(bio, mean = bint[1] + bslope[1] * y, sd = bse[1]) * 
        dnorm(y, mean = ymean[1], sd = ysd[1])
    temp[, 3] <-
      dnorm(bio, mean = bint[2] + bslope[2] * y, sd = bse[2]) * 
        dnorm(y, mean = ymean[2], sd = ysd[2])
    temp
  }

  # function to evaluate the incomplete negative log likelihood.
  incomp.nll <- function(pars){
    theta <- pars[-length(pars)]
    p.hat <- pars[length(pars)]
    k <- cbind(known, 
      (1 - known) * p.hat, 
      (1 - known) * (1 - p.hat))
   -sum(log(apply(k * dens(theta), 1, sum)))
  }

  k <- cbind(known, 
    (1 - known) * p.hat, 
    (1 - known) * (1 - p.hat))
    nll <- -sum(log(apply(k * dens(theta), 1, sum)))

  # first E-step
  dens_ <- t(c(p.hat, 1- p.hat) * t(dens(theta)[, -1]))
  w <- dens_ / apply(dens_, 1, sum)

  nll.diff <- 1
  iter <- 0		

  # EM algorithm
  while(nll.diff > tol){
    tt <- numeric(length(theta))
    # M-step
    tt[1] <- weighted.mean(y, (1 - known) * w[, 1]) 
    tt[2] <- weighted.mean(y, (1 - known) * w[, 2])
    tt[3:4] <- weighted.sd(y, tt[1], tt[2], 
      (1 - known) * w[, 1], 
      (1 - known) * w[, 2])
    clm <- lm(bio ~ y, weights = known + (1 - known) * w[, 1])
    nclm <- lm(bio ~ y, weights = (1 - known) * w[, 2])
    tt[5] <- coef(clm)[1]
    tt[6] <- coef(nclm)[1]
    tt[7] <- coef(clm)[2]
    tt[8] <- coef(nclm)[2]
    tt[9:10] <- weighted.sd(bio, 
      fitted(clm), 
      fitted(nclm), 
      known + (1 - known) * w[, 1],
      (1 - known) * w[, 2])
    theta <- tt

    if (sum(known) > 0) {
      p.hat <- mean(w[-(1:sum(known)), 1])
    } else {
      p.hat <- mean(w[, 1])
    } 

    k <- cbind(known,
      (1 - known) * p.hat,
      (1 - known) * (1 - p.hat))

    # E-step
    dens_ <- t(c(p.hat, 1- p.hat) * t(dens(theta)[, -1]))
    w <- dens_ / apply(dens_, 1, sum)

    nll.update <- incomp.nll(c(theta, p.hat))
    nll.diff <- nll - nll.update

    nll <- nll.update
    iter <- iter + 1
  }

  out <- theta

  theta.hat <- matrix(out, nrow = 2)
  rownames(theta.hat) <- c("compliant", "non.compliant")
  nn <- c("ymean", "ysd", "bint", "bslope", "bsd")
  colnames(theta.hat) <- nn

  #  posterior probability of compliance
  dens_ <- t(c(p.hat, 1- p.hat) * t(dens(theta)[, -1]))
  post.prob <- dens_[,1]/apply(dens_, 1, sum)

  list(theta.hat = theta.hat, 
    prob.compliance = p.hat[1], 
    post.prob = post.prob)
}

failwith <- function(default = NULL, f, quiet = TRUE) {
  # borrowed from dplyr package. 
  function(...) {
    out <- default
    try(out <- f(...), silent = TRUE)
    out
  }
}

weighted.sd <- function(x, mean.1, mean.2, weight.1, weight.2) {
  # returns weighted sd
  sum.w <- sum(weight.1) + sum(weight.2)
  sqrt(sum(weight.1 * (x - mean.1) ^ 2 + weight.2 * (x - mean.2) ^ 2) / sum.w)
}
