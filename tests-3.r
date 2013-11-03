source("try-3.r")
set.seed(10000)

if (!file.exists("simres-3.rdata")) {
  simres <- simExpGammaTrial(1000, 200, 
                             list(p0bar = -2, sigmasq_p0 = 0.5, 
                                  rbar = -3.5, sigmasq_r = 0.5, 
                                  alphabar = 3.5, sigmasq_alpha = 0.5, 
                                  betabar = matrix(0.5, 1, 1), Sigma_beta = 0.5 * diag(1)))
                                  #betabar = matrix(c(0.5, 0.3), 2, 1), Sigma_beta = 0.5 * diag(2)))
  simres2 <- simExpGammaTrial(1000, 200,
                              list(p0bar = -2, sigmasq_p0 = 0.5,
                                   rbar = -3.5, sigmasq_r = 0.5,
                                   alphabar = 3.5, sigmasq_alpha = 0.5,
                                   betabar = matrix(c(0.5, 0.3), 2, 1), Sigma_beta = 0.5 * diag(2)))  
  save(simres, simres2, file = "simres-3.rdata")
} else {
  load("simres-3.rdata")
}

# Test the likelihood function

objfn <- function(params, y, X) {
  params[1] <- invlogit(params[1])
  params[-1] <- exp(params[-1])
  params <- c(as.list(params[1:3]), list(params[-(1:3)]))
  names(params) <- c("p0", "r", "alpha", "beta")
  if (!identical(length(params$beta), ncol(X)))
    stop("Dimension mismatch between beta and X")
  ll <- loglikelihood(params, y, X, 1)
  #print(ll)
  ll
}


mle1 <- function(lst){
  res <- vector("list", lst$sample_size)
  for (i in seq(lst$sample_size)) {
    if (i %% 10 == 0) cat(i, " done \n")
    yi <- lst$data[id == i, y]
    sel <- grep("^X", names(lst$data), value = TRUE)
    Xi <- lst$data[id == i, sel, with = FALSE]
    Xi <- as.matrix(Xi)
    timei <- lst$data[id == i, time]
    res[[i]] <- try(optim(c(-1, -1, -1, -1), objfn, y = yi, X = Xi, control = list(fnscale = -1)))
  }
  res
}

reportpars <- function(mle) {
  par <- mle$par
  p0 <- invlogit(par[1])
  r <- exp(par[2])
  alpha <- exp(par[3])
  beta <- exp(par[-(1:3)])
  list(p0 = p0, r = r, alpha = alpha, beta = beta)
}

curvefit <- function(mle, y, X){
  res <- reportpars(mle)
  fit <- cdf(res, X)
  data.table(ObsTrialPct = y, PredTrialPct = fit)
}

test_oneid <- function(sim, idval, ...) {
  p0 <- sim$indPars$p0[idval]
  r <- sim$indPars$r[idval]
  alpha <- sim$indPars$alpha[idval]
  beta <- sim$indPars$beta[idval,]
  y <- sim$data[id == idval, y]
  time <- sim$data[id == idval, time]
  Xsel <- grep("^X", names(sim$data), value = TRUE)
  X <- sim$data[id == idval, Xsel, with = FALSE]
  X <- as.matrix(X)
  
  # Calculate likelihood at 'true' parameter values
  ll <- loglikelihood(list(p0 = p0, r = r, alpha = alpha, beta = beta), y = y, X = X, 1)
  
  # Run MLE 
  mlstart <- c(rep(-1, 3), rep(-1, length(beta)))
  mle <- optim(mlstart, objfn, y = y, X = X, control = list(fnscale = -1, maxit = 20000), ...)
  
  # Print comparisons
  cat("Statistic \tTrue \tEstimate\n")
  cat("---------------------------------------\n")
  cat("LogLik:\t", ll, "\t", mle$value, "\n")
  cat("p0:\t", p0, "\t", invlogit(mle$par[1]), "\n")
  cat("r:\t", r, "\t", exp(mle$par[2]), "\n")
  cat("alpha:\t", alpha, "\t", exp(mle$par[3]), "\n")
  cat("beta:\t", beta, "\t", exp(mle$par[-(1:3)]), "\n")
  invisible(list(ll = ll, mle = mle))
}


test_oneid <- function(idval, sim, print = TRUE, ...) {
  p0 <- sim$indPars$p0[idval]
  r <- sim$indPars$r[idval]
  alpha <- sim$indPars$alpha[idval]
  beta <- sim$indPars$beta[idval,]
  y <- sim$data[id == idval, y]
  time <- sim$data[id == idval, time]
  Xsel <- grep("^X", names(sim$data), value = TRUE)
  X <- sim$data[id == idval, Xsel, with = FALSE]
  X <- as.matrix(X)
  
  # Calculate likelihood at 'true' parameter values
  ll <- loglikelihood(list(p0 = p0, r = r, alpha = alpha, beta = beta), y = y, X = X, 1)
  
  # Run MLE 
  mlstart <- c(rep(-1, 3), rep(-1, length(beta)))
  mle <- optim(mlstart, objfn, y = y, X = X, control = list(fnscale = -1, maxit = 20000), ...)
  
  # Create list of vectors with comparisons
  complist <- list(loglik = c(ll, mle$value), 
                   p0 = c(p0, invlogit(mle$par[1])), 
                   r = c(r, exp(mle$par[2])), 
                   alpha = c(alpha, exp(mle$par[3])),
                   beta = c(beta, exp(mle$par[-(1:3)])))
  
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

fitplot <- function(x, y, ...) {
  plot(x, y, xlab = "True", ylab = "Estimated", pch = ".", cex = 4, ...)
  abline(a = 0, b = 1, col = 2)
}

test_allids <- function(sim, ...){  
  res <- lapply(seq(sim$sample_size), test_oneid, sim = sim, print = FALSE, ...)
  loglik <- do.call("rbind", lapply(res, function(x) x$complist$loglik))
  p0 <- do.call("rbind", lapply(res, function(x) x$complist$p0))
  r <- do.call("rbind", lapply(res, function(x) x$complist$r))
  alpha <- do.call("rbind", lapply(res, function(x) x$complist$alpha))
  beta <- do.call("rbind", lapply(res, function(x) x$complist$beta))
  beta <- array(beta, dim = c(nrow(beta), ncol(beta)/2, 2))
  mlestatus <- do.call("c", lapply(res, function(x) x$mle$convergence))
  nr <- ceiling((3 + dim(beta)[2])/2)
  par(mfrow = c(nr,2))
  fitplot(loglik[,1], loglik[,2], main = "Loglikelihood")
  fitplot(p0[,1], p0[,2],main = "p0")
  fitplot(r[,1], r[,2], main = "r")
  fitplot(alpha[,1], alpha[,2], main = "alpha")
  for (j in seq(dim(beta)[2]))
    fitplot(beta[,j,1], beta[,j,2], main = paste0("beta",j))
  return(list(loglik = loglik, p0 = p0, r = r, alpha = alpha, beta = beta, mlestatus = mlestatus))
  par(mfrow = c(1,1))
}

test_oneid(1, simres)
testres <- test_allids(simres, method = "BFGS")

testres1 <- test_allids(simres)
testres2 <- test_allids(simres2)

save(testres1, testres2, file = "simtests-3.rdata")

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