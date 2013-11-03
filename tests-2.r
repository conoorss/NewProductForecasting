source("try-2.r")
set.seed(10000)

if (!file.exists("simres-2.rdata")) {
	simres <- simExpGammaTrial(1000, 200, 
                             list(p0bar = -2, sigmasq_p0 = 0.5, 
                                  rbar = -4, sigmasq_r = 0.5, 
                                  alphabar = 3, sigmasq_alpha = 0.5))
	save(simres, file = "simres-2.rdata")
} else {
	load("simres-2.rdata")
}

# Test the likelihood function

objfn <- function(params, y, time) {
	params[-1] <- exp(params[-1])
  params[1] <- invlogit(params[1])
	names(params) <- c("p0", "r", "alpha")
	ll <- loglikelihood(params, y, time, 1)
	#print(ll)
	ll
}

mle1 <- function(lst){
	res <- vector("list", lst$sample_size)
	for (i in seq(lst$sample_size)) {
		if (i %% 10 == 0) cat(i, " done \n")
		yi <- lst$data[id == i, y]
		timei <- lst$data[id == i, time]
		res[[i]] <- try(optim(c(-1, -1, -1), objfn, y = yi, time = timei, control = list(fnscale = -1)))
	}
	res
}

test_oneid <- function(idval, sim, print = TRUE, ...) {
  p0 <- sim$indPars$p0[idval]
  r <- sim$indPars$r[idval]
  alpha <- sim$indPars$alpha[idval]
  y <- sim$data[id == idval, y]
  time <- sim$data[id == idval, time]
  
  # Calculate likelihood at 'true' parameter values
  ll <- loglikelihood(list(p0 = p0, r = r, alpha = alpha, beta = beta), y = y, time = time, 1)
  
  # Run MLE 
  mlstart <- rep(-1, 3)
  mle <- optim(mlstart, objfn, y = y, time = time, control = list(fnscale = -1, maxit = 20000), ...)
  
  # Create list of vectors with comparisons
  complist <- list(loglik = c(ll, mle$value), 
                   p0 = c(p0, invlogit(mle$par[1])), 
                   r = c(r, exp(mle$par[2])), 
                   alpha = c(alpha, exp(mle$par[3])))
      
  # Print comparisons
  if (print) {
    oopt <- options()
    options(digits = 3)
    on.exit(options(oopt))
    compmat <- do.call("rbind", complist)
    colnames(compmat) <- c("True", "Estimate")
    cat("==== Comparison =====\n")
    print(compmat)
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
  mlestatus <- do.call("c", lapply(res, function(x) x$mle$convergence))
  
  par(mfrow = c(2,2))
  fitplot(loglik[,1], loglik[,2], main = "Loglikelihood")
  fitplot(p0[,1], p0[,2],main = "p0")
  fitplot(r[,1], r[,2], main = "r")
  fitplot(alpha[,1], alpha[,2], main = "alpha")
  return(list(loglik = loglik, p0 = p0, r = r, alpha = alpha, mlestatus = mlestatus))
  par(mfrow = c(1,1))
}

test_oneid(1, simres)
testres <- test_allids(simres, method = "BFGS")


if (FALSE) {
tcurve <- function(r, alpha, p0 = 1) p0 * (1 - (alpha / (alpha + seq(200)))^r)
t1 <- tcurve(r = simres$indPars$r[1], alpha = simres$indPars$alpha[1])
t2 <- do.call(tcurve, as.list(exp(mle.test$par)))
matplot(cbind(t1, t2), type = "l", lwd = c(1,2))

#plotSimData(simres)

#res <- gibbsSampler(simres, Mcmc = list(burnin = 100, samples = 100, thin = 10, printThin = 10, rwscale_r = 0.02, rwscale_alpha = 0.01, 
#										startVals = list(r = rep(1, 1000), alpha = rep(1, 1000), rbar = 0, alphabar = 0, sigmasq_r = 1, sigmasq_alpha = 1)))

res <- gibbsSampler1(simres, 
					 Mcmc = list(burnin = 50000, samples = 10000, thin = 10, printThin = 1000, 
								 rwscale_r = 4, rwscale_alpha = 1, 
								 startVals = list(r = r1, alpha = 1),
								 rFixed = TRUE, alphaFixed = FALSE))


# Test mle with alpha fixed to true value
loglikelihood(c(r = r1, alpha = alpha1), y = y1, time = time1, 1)
mle.test <- optim(c(-1, -1), objfn, y = y1, time = time1, control = list(fnscale = -1))
mle.test2 <- optim(c(-1), objfn, y = y1, time = time1, alpha = alpha1, control = list(fnscale = -1))
}
