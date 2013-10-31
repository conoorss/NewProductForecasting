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


# Tests for first id
p01 <- simres$indPars$p0[1]
r1 <- simres$indPars$r[1]
alpha1 <- simres$indPars$alpha[1]

y1 <- simres$data[id == 1, y]
time1 <- simres$data[id == 1, time]

loglikelihood(c(p0 = p01, r = r1, alpha = alpha1), y = y1, time = time1, 1)
mle.test <- optim(c(-1, -1, -1), objfn, y = y1, time = time1, control = list(fnscale = -1))

mle.test <- mle1(simres)
conv <- do.call("c", lapply(mle.test, function(x) x$convergence))
estpars <- do.call("rbind", lapply(mle.test, function(x) x$par))
estpars[,-1] <- exp(estpars[,-1])
estpars[,1] <- invlogit(estpars[,1])
colnames(estpars) <- c("p0.hat", "r.hat", "alpha.hat")

fitplot <- function(x, y, param) { 
  plot(x, y, main = paste("Fit of", param), xlab = "Obs", ylab = "Est", pch = ".", cex = 3)
  abline(a = 0, b = 1, col = 2)
}

par(mfrow = c(3,1)) 
fitplot(simres$indPars$p0, estpars[,"p0.hat"], "p0")
fitplot(simres$indPars$r, estpars[,"r.hat"], "r")
fitplot(simres$indPars$alpha, estpars[,"alpha.hat"], "alpha")
par(mfrow = c(1,1))





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
