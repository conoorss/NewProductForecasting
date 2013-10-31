#
# Model Spec:
#	y_{jt} = (1 - (alpha_{j}/(alpha_{j} + t))^r_{j})
#
#	log(alpha) ~ N(alphabar, sigmasq_alpha)
#	log(r) ~ N(rbar, sigmasq_r)
#

source("runireg.r")

simExpGammaTrial <- function(sample_size, horizon, simPars, seed = 1234) {
	require(data.table)
	r <- exp(rnorm(sample_size, mean = simPars$rbar, sd = sqrt(simPars$sigmasq_r)))
	alpha <- exp(rnorm(sample_size, mean = simPars$alphabar, sd = sqrt(simPars$sigmasq_alpha)))
	data <- vector("list", sample_size)
	time <- seq(horizon)
	for (i in seq(sample_size)) {
		y <- (1 - (alpha[i] / (alpha[i] + time))^r[i])
		data[[i]] <- data.table(id = i, time = time, y = y)
	}
	data <- rbindlist(data)
	list(sample_size = sample_size, 
       potential = rep(1, sample_size), 
       horizon = horizon, 
       simPars = simPars, 
       indPars = list(r = r, alpha = alpha), 
       data = data)
}

plotSimData <- function(lst) {
	mat <- matrix(lst$data$y, nrow = lst$horizon)
	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
	ix <- seq(ncol(mat)) %% 10 == 0
	matplot(mat[,ix], type = "l", col = 1, main = "Trial curves")
	plot(density(lst$indPars$r), xlab = "", ylab = "", main = "Density of r")
	plot(density(lst$indPars$alpha), xlab = "", ylab = "", main = "Density of alpha")
	#invisible()
}

gibbsSampler1 <- function(Data, Priors, Mcmc) {
	# Sampler for one unit
   	require(coda)
	default.prior <- list(beta = c(0), sigmasq = c(10))
	if (missing(Priors))
		Priors <- list(r = default.prior, alpha = default.prior)
	if (is.null(Priors$r))
		Priors$r <- default.prior
	if (is.null(Priors$alpha))
		Priors$alpha <- default.prior

	stored <- Mcmc$samples / Mcmc$thin
	output <- list(r = mcmc(matrix(NA, stored, 1)), alpha = mcmc(matrix(NA, stored, 1)))
	r <- Mcmc$startVals$r
	alpha <- Mcmc$startVals$alpha

	alpha.rejCount <- r.rejCount <- 0

	printStatus <- function() {
		curr.time <- proc.time()[3]
		cat("Iter: ", iter, "-- Time Elapsed: ", round((curr.time - start.time) / 60, 2), "mins", "\n")
		cat("r RejRate: ", round(100 * r.rejCount / iter), "\n")
		cat("alpha RejRate: ", round(100 * alpha.rejCount / iter), "\n")
		cat("\n")
	}

	start.time <- proc.time()[3]
	
	curr_id <- 1
	dati <- Data$data[id ==curr_id]
	for (iter in seq(Mcmc$burnin + Mcmc$samples)){
		if (!Mcmc$rFixed) {
			rres <- sample.r(r, alpha, Priors$r$beta, Priors$r$sigmasq, Mcmc$rwscale_r, dati$y, dati$time, Data$potential[curr_id])
			r <- rres$r
			r.rejCount <- r.rejCount + rres$reject
		}
		
		if (!Mcmc$alphaFixed) {
			alphares <- sample.alpha(alpha, r, Priors$alpha$beta, Priors$alpha$sigmasq, Mcmc$rwscale_alpha, dati$y, dati$time, Data$potential[curr_id])
			alpha <- alphares$alpha
			alpha.rejCount <- alpha.rejCount + alphares$reject
		}

		if (iter %% Mcmc$printThin == 0) printStatus()
		#cat(iter, r, alpha, "\n")
		if (iter > Mcmc$burnin && ((iter - Mcmc$burnin)  %% Mcmc$thin == 0)) {
			j <- (iter - Mcmc$burnin) / Mcmc$thin
			output$r[j, 1] <- r
			output$alpha[j, 1] <- alpha
		}
	}

	return(output)
}

gibbsSampler <- function(Data, Priors, Mcmc) {
	require(coda)
	#
	# Data: list with components sample_size, horizon, data
	# 		data must have columns y and time
	#
	# Priors: list with components alphabar, sigmasq_alpha, rbar, sigmasq_r
	#		Each of the above is a list with parameters for the corresponding model parameter
	# 
	# Mcmc: list with components burnin, stored, thin, printThin, alphaTuning, rTuning, startVals
	#		burnin, samples, thin, and printThin are integers
	#		alphaTuning and rTuning are numeric
	#		startVals is a list with named components matching the parameters r, alpha, rbar, sigmasq_r, alphabar, sigmasq_alpha
	#		seed - random seed

	# Check Data inputs
	ids <- Data$data[,unique(id)]
	# Check Prior inputs
	default.prior <- list(betabar = c(0), A = 0.01 * diag(1L), nu = 3, ssq = 10)
	if (missing(Priors))
		Priors <- list(r = default.prior, alpha = default.prior)
	if (is.null(Priors$r))
		Priors$r <- default.prior
	if (is.null(Priors$alpha))
		Priors$alpha <- default.prior

	# Check MCMC inputs

	# Return an object of class mcmc defined in coda
	stored <- Mcmc$samples / Mcmc$thin
	output <- list(r = mcmc(matrix(NA, stored, Data$sample_size)), alpha = mcmc(matrix(NA, stored, Data$sample_size)), 
				   hyperPars = mcmc(matrix(NA, stored, 4, dimnames = list(NULL, c("rbar", "sigmasq_r", "alphabar", "sigmasq_alpha")))))

	r <- Mcmc$startVals$r
	alpha <- Mcmc$startVals$alpha
	rbar <- Mcmc$startVals$rbar
	alphabar <- Mcmc$startVals$alphabar
	sigmasq_r <- Mcmc$startVals$sigmasq_r
	sigmasq_alpha <- Mcmc$startVals$sigmasq_alpha

	alpha.rejCount <- r.rejCount <- rep(0, length(r))

	printStatus <- function() {
		curr.time <- proc.time()[3]
		cat("Iter: ", iter, "-- Time Elapsed: ", round((curr.time - start.time) / 60, 2), "mins", "\n")
		cat("r RejRate: \n")
		print(summary(100 * r.rejCount / iter))
		cat("alpha RejRate: \n")
		print(summary(100 * alpha.rejCount / iter))			
		cat("\n")
	}

	start.time <- proc.time()[3]

	for (iter in seq(Mcmc$burnin + Mcmc$samples)){
		# Sample individual level parameters
		### Need to make this parallel
		for (i in seq(Data$sample_size)) {
			dati <- Data$data[id == ids[i]]
			if (!Mcmc$rFixed) {
				rres <- sample.r(r[i], alpha[i], rbar, sigmasq_r, Mcmc$rwscale_r, dati$y, dati$time, Data$potential[i])
				r[i] <- rres$r
				r.rejCount[i] <- r.rejCount[i] + rres$reject
			}

			if (!alphaFixed) {
				alphares <- sample.alpha(alpha[i], r[i], alphabar, sigmasq_alpha, Mcmc$rwscale_alpha, dati$y, dati$time, Data$potential[i])
				alpha[i] <- alphares$alpha
				alpha.rejCount[i] <- alpha.rejCount[i] + alphares$reject
			}
		}
		# Sample hyper parameters
#		if (missing(Priors) || is.null(Priors$r)) 
#			rhres <- sample.hpars(r) 
#		else 
#			rhres <- sample.hpars(r, Priors$r)
		rhres <- sample.hpars(list(y = log(r), X = matrix(1, nrow = length(r), ncol = 1)), Priors$r)
		rbar <- rhres$beta[,]
		sigmasq_r <- rhres$sigmasq[,]

#		if (missing(Priors) || is.null(Priors$alpha)) 
#			alphahres <- sample.hpars(alpha) 
#		else 
#			alphahres <- sample.hpars(alpha, Priors$alpha)
		alphahres <- sample.hpars(list(y = log(alpha), X = matrix(1, nrow = length(alpha), ncol = 1)), Priors$alpha)
		alphabar <- alphahres$beta[,]
		sigmasq_alpha <- alphahres$sigmasq[,]

		# Print/Store results
		#if (iter %% Mcmc$printThin == 0) printStatus()
		cat(iter, rbar, alphabar, sigmasq_r, sigmasq_alpha, "\n")
		cat(iter, sum(is.na(r)), min(r, na.rm = TRUE), median(r, na.rm = TRUE), max(r, na.rm = TRUE),  "\n")
		cat(iter, sum(is.na(alpha)), min(alpha, na.rm = TRUE), median(alpha, na.rm = TRUE), max(alpha, na.rm = TRUE),  "\n")
		if (iter > Mcmc$burnin && ((iter - Mcmc$burnin)  %% Mcmc$thin == 0)) {
			j <- (iter - Mcmc$burnin) / Mcmc$thin
			output$r[j, ] <- r
			output$alpha[j, ] <- alpha
			output$hyperPars[j, "rbar"] <- rbar
			output$hyperPars[j, "alphabar"] <- alphabar
			output$hyperPars[j, "sigmasq_r"] <- sigmasq_r
			output$hyperPars[j, "sigmasq_alpha"] <- sigmasq_alpha
		}
	}
	return(output)
}

sample.r <- function(r, alpha, rbar, sigmasq_r, rwscale_r, ...){
	r.prop <- propose_param(r, rwscale_r)
	loglik <- loglikelihood(c(r = r, alpha = alpha), ...)
	logprior <- dnorm(log(r), mean = rbar, sd = sqrt(sigmasq_r), log = TRUE) - log(r)
	loglik.prop <- loglikelihood(c(r = r.prop, alpha = alpha), ...)
	logprior.prop <- dnorm(log(r.prop), mean = rbar, sd = sqrt(sigmasq_r), log = TRUE) - log(r.prop)
	mhratio <- exp(loglik.prop + logprior.prop - loglik - logprior)
	mhratio <- exp(loglik.prop - loglik)
	threshold <- runif(1)
	if (FALSE) {
	cat("=====================================\n")
	cat(r, r.prop, "\n")
	cat(loglik, loglik.prop, "\n")
	cat(mhratio, threshold, "\n")
	}
	if (mhratio > threshold)
		list(r = r.prop, reject = FALSE)
	else
		list(r = r, reject = TRUE)
}

sample.alpha <- function(alpha, r, alphabar, sigmasq_alpha, rwscale_alpha, ...) {
	alpha.prop <- propose_param(alpha, rwscale_alpha)
	loglik <- loglikelihood(c(r = r, alpha = alpha), ...)
	logprior <- dnorm(log(alpha), mean = alphabar, sd = sqrt(sigmasq_alpha), log = TRUE) - log(alpha)
	loglik.prop <- loglikelihood(c(r = r, alpha = alpha.prop), ...)
	logprior.prop <- dnorm(log(alpha.prop), mean = alphabar, sd = sqrt(sigmasq_alpha), log = TRUE) - log(alpha.prop)
#	mhratio <- exp(loglik.prop + logprior.prop - loglik - logprior)
	mhratio <- exp(loglik.prop - loglik)
	threshold <- runif(1)
	if (FALSE) {
	cat("**************************************\n")
	cat(alpha, alpha.prop, "\n")
	cat(loglik, loglik.prop, "\n")
	cat(mhratio, threshold, "\n")
	}
	if (mhratio > runif(1))
		list(alpha = alpha.prop, reject = FALSE)
	else
		list(alpha = alpha, reject = TRUE)
}

cdf <- function(pars, tau) {
  r <- pars[["r"]]
  alpha <- pars[["alpha"]]
  (1 - (alpha / (alpha + tau))^r)
}


loglikelihood <- function(pars, y, time, potential) {
	incy <- diff(c(0, y))
	cdf1 <- cdf(pars, time)
	cdf2 <- cdf(pars, time - 1)
	cdf3 <- cdf(pars, max(time))

	ll <- sum(incy * log(cdf1 - cdf2)) + (potential - sum(incy)) * log(1 - cdf3)
	ll
}

propose_param <- function(param, scale) exp(log(param) + scale * rnorm(length(param)))

stripMcmcAttributes <- function(x) {
	a <- names(attributes(x))
	a <- setdiff(a, "dim")
	for (j in a) attr(x, j) <- NULL
	x
}


#sample.hpars <- function(indparam, prior){
#	require(bayesm)
#	res <- runireg(list(y = indparam, X = matrix(1, nrow = length(indparam), ncol = 1)), prior, list(R = 1))	
#	res <- lapply(res, stripMcmcAttributes)
#}
#
sample.hpars <- function(Data, Prior) {
	nobs <- length(Data$y)
	nvar <- ncol(Data$X)

	betabar <- Prior$betabar
	nu <- Prior$nu
	ssq <- Prior$ssq
	
	
	RA <- chol(Prior$A)
	W <- rbind(Data$X, RA)
	z <- c(Data$y, as.vector(RA %*% Prior$betabar))
	IR <- backsolve(chol(crossprod(W)), diag(nvar))
	btilde <- crossprod(t(IR)) %*% crossprod(W, z)
	res <- z - W %*% btilde
	s <- crossprod(res)
	sigmasq <- with(Prior, (nu * ssq + s) / rchisq(1, nu + nobs))
	beta <- btilde + as.vector(sqrt(sigmasq)) * IR %*% rnorm(nvar)
	list(beta = beta, sigmasq = sigmasq)
}

