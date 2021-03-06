source("try-4.r")
set.seed(10000)

if (!file.exists("simres-4.rdata")) {
  simres <- simExpGammaTrial(1000, 200, 
                             list(rbar = -3.5, sigmasq_r = 0.5, 
                                  alphabar = 3.5, sigmasq_alpha = 0.5, 
                                  betabar = matrix(0.5, 1, 1), Sigma_beta = 0.5 * diag(1)))
                                  #betabar = matrix(c(0.5, 0.3), 2, 1), Sigma_beta = 0.5 * diag(2)))
  simres2 <- simExpGammaTrial(1000, 200,
                              list(rbar = -3.5, sigmasq_r = 0.5,
                                   alphabar = 3.5, sigmasq_alpha = 0.5,
                                   betabar = matrix(c(0.5, 0.3), 2, 1), Sigma_beta = 0.5 * diag(2)))  
  save(simres, simres2, file = "simres-4.rdata")
} else {
  load("simres-4.rdata")
}

# Test the likelihood function

objfn <- function(params, y, X) {
  params <- exp(params)
  params <- c(as.list(params[1:2]), list(params[-(1:2)]))
  names(params) <- c("r", "alpha", "beta")
  if (!identical(length(params$beta), ncol(X)))
    stop("Dimension mismatch between beta and X")
  ll <- loglikelihood(params, y, X, 1)
  #print(ll)
  ll
}

fitplot <- function(x, y, ...) {
  plot(x, y, xlab = "True", ylab = "Estimated", pch = ".", cex = 4, ...)
  abline(a = 0, b = 1, col = 2)
}

test_oneid <- function(idval, sim, print = TRUE, ...) {
  r <- sim$indPars$r[idval]
  alpha <- sim$indPars$alpha[idval]
  beta <- sim$indPars$beta[idval,]
  y <- sim$data[id == idval, y]
  time <- sim$data[id == idval, time]
  Xsel <- grep("^X", names(sim$data), value = TRUE)
  X <- sim$data[id == idval, Xsel, with = FALSE]
  X <- as.matrix(X)
  
  # Calculate likelihood at 'true' parameter values
  ll <- loglikelihood(list(r = r, alpha = alpha, beta = beta), y = y, X = X, 1)
  
  # Run MLE 
  mlstart <- c(rep(-1, 2), rep(-1, length(beta)))
  mle <- optim(mlstart, objfn, y = y, X = X, control = list(fnscale = -1, maxit = 20000), ...)
  
  # Create list of vectors with comparisons
  complist <- list(loglik = c(ll, mle$value), 
                   r = c(r, exp(mle$par[1])), 
                   alpha = c(alpha, exp(mle$par[2])),
                   beta = c(beta, exp(mle$par[-(1:2)])))
  
  # Print comparisons
  if (print) {
    oopt <- options()
    options(digits = 3)
    on.exit(options(oopt))
    compmat <- do.call("rbind", complist[-length(complist)])
    compmat2 <- matrix(complist[[length(complist)]], nrow = length(beta), ncol = 2)
    rownames(compmat2) <- paste0("beta", seq_along(beta))
    colnames(compmat) <- c("True", "Estimate")
    cat("==== Comparison =====\n")
    print(rbind(compmat, compmat2))
    cat("MLE Convergence: ", mle$convergence, "\n")
  }
  invisible(list(complist = complist, mle = mle))
}

test_allids <- function(sim, ...){  
  res <- lapply(seq(sim$sample_size), test_oneid, sim = sim, print = FALSE, ...)
  loglik <- do.call("rbind", lapply(res, function(x) x$complist$loglik))
  r <- do.call("rbind", lapply(res, function(x) x$complist$r))
  alpha <- do.call("rbind", lapply(res, function(x) x$complist$alpha))
  beta <- do.call("rbind", lapply(res, function(x) x$complist$beta))
  beta <- array(beta, dim = c(nrow(beta), ncol(beta)/2, 2))
  mlestatus <- do.call("c", lapply(res, function(x) x$mle$convergence))
  nr <- ceiling((3 + dim(beta)[2])/2)
  par(mfrow = c(nr,2))
  fitplot(loglik[,1], loglik[,2], main = "Loglikelihood")
  fitplot(r[,1], r[,2], main = "r")
  fitplot(alpha[,1], alpha[,2], main = "alpha")
  for (j in seq(dim(beta)[2]))
    fitplot(beta[,j,1], beta[,j,2], main = paste0("beta",j))
  return(list(loglik = loglik, r = r, alpha = alpha, beta = beta, mlestatus = mlestatus))
  par(mfrow = c(1,1))
}

test_oneid(1, simres)
testres1 <- test_allids(simres, method = "BFGS")
testres2 <- test_allids(simres2, method = "BFGS")

save(testres1, testres2, file = "simtests-4.rdata")

#### NOT YET DONE ####

if (FALSE) {
mle.test <- mle1(simres)
conv <- do.call("c", lapply(mle.test, function(x) x$convergence))
estpars <- do.call("rbind", lapply(mle.test, function(x) x$par))
estpars[,-1] <- exp(estpars[,-1])
estpars[,1] <- invlogit(estpars[,1])
colnames(estpars) <- c("p0.hat", "r.hat", "alpha.hat")

fitplot <- function(x, y, param) { 
  plot(x, y, main = paste("Fit of", param), xlab = "Obs", ylab = "Est", pch = ".", cex = 4)
  abline(a = 0, b = 1, col = 2)
}

opar <- par(no.readonly = TRUE)
par(mfrow = c(3,1), mar = c(4.1, 4.1, 4.1, 2.1)) 
fitplot(simres$indPars$p0, estpars[,"p0.hat"], "p0")
fitplot(simres$indPars$r, estpars[,"r.hat"], "r")
fitplot(simres$indPars$alpha, estpars[,"alpha.hat"], "alpha")
par(opar)
}