source("stewart.R")

sigmasq_r = 0.5 
sigmasq_alpha = 1 
Sigma.prior = matrix(c(sigmasq_r, 0, 0, sigmasq_alpha), ncol = 2)
simres <- simExpGammaTrial(1, 200, list(rbar = -2, sigmasq_r = sigmasq_r, alphabar = -4, sigmasq_alpha = sigmasq_alpha))

#if (!file.exists("simres-1.rdata")) {
	#simres <- simExpGammaTrial(1, 200, list(rbar = -2, sigmasq_r = sigmasq_r, alphabar = -4, sigmasq_alpha = sigmasq_alpha))
#	save(simres, file = "simres-1.rdata")
#} else {
	#load("simres-1.rdata")
#}

# Test the likelihood function

objfn <- function(params, y, time) {
	params <- exp(params)
	names(params) <- c("r", "alpha")
	ll <- loglikelihood(params, y, time, 1)
	print(ll)
	ll
}

c(simres$indPars$r[1], simres$indPars$alpha[1])
c(log(simres$indPars$r[1]),log(simres$indPars$alpha[1]))

U(c(log(simres$indPars$r[1]),log(simres$indPars$alpha[1])), y = simres$data[id == 1, y], time = simres$data[id == 1, time], 1)
loglikelihood(c(r = simres$indPars$r[1], alpha = simres$indPars$alpha[1]), y = simres$data[id == 1, y], time = simres$data[id == 1, time], 1)

mle.test <- optim(c(-1, -1), objfn, y = simres$data[id == 1, y], time = simres$data[id == 1, time], control = list(fnscale = -1))
mle.test <- optim(c(-1, -1), objfn, y = simres$data[id == 1, y], time = simres$data[id == 1, time], control = list(fnscale = -1), method = "BFGS")

tcurve <- function(r, alpha, p0 = 1) p0 * (1 - (alpha / (alpha + seq(200)))^r)
t1 <- tcurve(r = simres$indPars$r[1], alpha = simres$indPars$alpha[1])
t2 <- do.call(tcurve, as.list(exp(mle.test$par)))
matplot(cbind(t1, t2), type = "l", lwd = c(1,2))

#plotSimData(simres)

#res <- gibbsSampler(simres, Mcmc = list(burnin = 100, samples = 100, thin = 10, printThin = 10, rwscale_r = 0.02, rwscale_alpha = 0.01, 
#										startVals = list(r = rep(1, 1000), alpha = rep(1, 1000), rbar = 0, alphabar = 0, sigmasq_r = 1, sigmasq_alpha = 1)))

source("stewart.R")
Sigma = matrix(c(100,0,0,100), ncol = 2)
res <- HMC(simres, 
					 Mcmc = list(burnin = 2000, samples = 2000, thin = 1, printThin = 1, 
                                epsilon = 1, L = 20, Sigma = Sigma,
								 startVals = list(r = -20, alpha = -20) ) )
								# startVals = list(r = log( simres$indPars$r[1]), alpha = log( simres$indPars$alpha[1]))) )


