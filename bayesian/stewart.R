#
# Model Spec:
#	y_{jt} = (1 - (alpha <- {j}/(alpha_{j} + t))^r_{j})
#
#	alpha ~ N(alphabar, sigmasq_alpha)
#	r ~ N(rbar, sigmasq_r)
#

source("runireg.r")
library(data.table)

simExpGammaTrial <- function(sample_size, horizon, simPars, seed = 1234) {
	set.seed(1234)
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
	list(sample_size = sample_size, potential = rep(1, sample_size), horizon = horizon, simPars = simPars, indPars = list(r = r, alpha = alpha), seed = seed, data = data)
}

simExpGammaTrial2 <- function(sample_size, horizon, simPars, seed = 1234) {
	set.seed(1234)
	require(data.table)
	r <- exp(rnorm(sample_size, mean = simPars$rbar, sd = sqrt(simPars$sigmasq_r)))
	alpha <- exp(rnorm(sample_size, mean = simPars$alphabar, sd = sqrt(simPars$sigmasq_alpha)))
	data <- vector("list", sample_size)
	time <- seq(horizon)
    num.ind = 5000
	for (i in seq(sample_size)) {
        lambda = rgamma(num.ind, r[i], alpha[i])
        trial.time = rexp(num.ind, lambda)
        trial.dis.time = sapply(time, function(x) sum(trial.time > (x - 1) & trial.time < x ) )
        trial.dis.time[horizon] = sum(trial.time > horizon)
        incy = trial.dis.time / sum(trial.dis.time)
        y = cumsum(incy)
		#y <- (1 - (alpha[i] / (alpha[i] + time))^r[i])
		data[[i]] <- data.table(id = i, time = time, y = y)
	}
	data <- rbindlist(data)
	list(sample_size = sample_size, potential = rep(1, sample_size), horizon = horizon, simPars = simPars, indPars = list(r = r, alpha = alpha), seed = seed, data = data)
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

HMC <- function(Data, Priors, Mcmc) {
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
		cat("RejRate: ", round(100 * r.rejCount / iter), "\n")
		#cat("alpha RejRate: ", round(100 * alpha.rejCount / iter), "\n")
		cat("\n")
	}

	start.time <- proc.time()[3]
	
	curr_id <- 1
	dati <- Data$data[id ==curr_id]

    current_q = c(r, alpha)

##############################
##### HMC tunning parameters
    epsilon = Mcmc$epsilon
    L = Mcmc$L
##############################

    y = dati$y; time = dati$time; potential = Data$potential[curr_id]

    #Sigma = matrix(c(0.01, 0, 0, 0.1), ncol = 2)
    #Sigma = matrix(c(100, 0, 0, 100), ncol = 2)
    Sigma = Mcmc$Sigma
    sig  = t(chol(solve(Sigma)))

	for (iter in seq(Mcmc$burnin + Mcmc$samples)){


          cat("current q is: ", current_q, "\n")
          q = current_q
          #p = rnorm(length(q))  
          p = sig %*% rnorm(length(q))  
          #p = rnorm(length(q), 0 , 0.1)  
          #cat("new p is : ", p , "\n")
          current_p = p

          # Make a half step for momentum at the beginning

          gg =  grad_U(q, y, time, potential)
          #gg = gg + InvCov %*% (q - q_bar)
          p = p - epsilon * gg / 2

          # Alternate full steps for position and momentum

          for (i in 1:L)
          {
            # Make a full step for the position

            #q = q + epsilon * Sigma %*% p
            q = q + epsilon * p

            # Make a full step for the momentum, except at end of trajectory

            gg =  grad_U(q, y, time, potential)
            #gg = gg + InvCov %*% (q - q_bar)
            if (i!=L) p = p - epsilon * gg
          }

          # Make a half step for momentum at the end.
          gg =  grad_U(q, y, time, potential)
################################
##### if you have upper model, uncomment out the folloig line
          #gg = gg + InvCov %*% (q - q_bar)
          p = p - epsilon * gg / 2

          # Negate momentum at end of trajectory to make the proposal symmetric

          p = -p

          # Evaluate potential and kinetic energies at start and end of trajectory

          current_U = U(current_q, y, time, potential)
################################
##### if you have upper model, uncomment out the folloig line
          #current_U2 = current_U + t(current_q - q_bar) %*% InvCov %*% (current_q - q_bar) / 2
          current_K = sum(current_p^2) / 2
          #current_K = t(current_p) %*% Sigma %*% current_p / 2
          proposed_U = U(q, y, time, potential)
################################
##### if you have upper model, uncomment out the folloig line
          #proposed_U2 = proposed_U + t(q - q_bar) %*% InvCov %*% (q - q_bar) / 2
          #proposed_K = t(p) %*% Sigma %*% p / 2
          proposed_K = sum(p^2) / 2

          cat(current_U, current_K, proposed_U, proposed_K, "\n")
################################
##### if you have upper model, uncomment out the folloig line
          #cat(current_U2, current_K, proposed_U2, proposed_K, "\n")

          # Accept or reject the state at end of trajectory, returning either
          # the position at the end of the trajectory or the initial position

          threshold = runif(1) 
          #ratio = exp(current_U2-proposed_U2+current_K-proposed_K)
          ratio = exp(current_U-proposed_U+current_K-proposed_K)

          cat("===", ratio, threshold, "===", "\n")

          #if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
          if( is.na(ratio ) ) { cat("ratio is not a number, continue....\n") 
                r.rejCount <- r.rejCount + 1
              } else {
                  if ( ratio > threshold)
                  {
                    current_q = q  # accept
                  } else {
                    r.rejCount <- r.rejCount + 1
                  }
          }

		if (iter %% Mcmc$printThin == 0) printStatus()
		#cat(iter, r, alpha, "\n")
		if (iter > Mcmc$burnin && ((iter - Mcmc$burnin)  %% Mcmc$thin == 0)) {
			j <- (iter - Mcmc$burnin) / Mcmc$thin
			output$r[j, 1] <- current_q[1] 
			output$alpha[j, 1] <-  current_q[2]
		}
	}

	return(output)
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

sample.r <- function(r, alpha, rbar, sigmasq_r, rwscale_r, ...){
	r.prop <- propose_param(r, rwscale_r)
	loglik <- loglikelihood(c(r = r, alpha = alpha), ...)
	logprior <- dnorm(log(r), mean = rbar, sd = sqrt(sigmasq_r), log = TRUE) - log(r)
	loglik.prop <- loglikelihood(c(r = r.prop, alpha = alpha), ...)
	logprior.prop <- dnorm(log(r.prop), mean = rbar, sd = sqrt(sigmasq_r), log = TRUE) - log(r.prop)
	mhratio <- exp(loglik.prop + logprior.prop - loglik - logprior)
	if (mhratio > runif(1))
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
	mhratio <- exp(loglik.prop + logprior.prop - loglik - logprior)
	if (mhratio > runif(1))
		list(alpha = alpha.prop, reject = FALSE)
	else
		list(alpha = alpha, reject = TRUE)
}

loglikelihood <- function(pars, y, time, potential) {
	cumcurve <- function(tau) (1 - (alpha / (alpha + tau))^r)
	r <- pars[["r"]]
	alpha <- pars[["alpha"]]
	incy <- diff(c(0, y))
	cdf1 <- cumcurve(time)
	cdf2 <- cumcurve(time - 1)
	cdf3 <- cumcurve(max(time))

	ll <- sum(incy * log(cdf1 - cdf2)) + (potential - sum(incy)) * log(1 - cdf3)
	ll
}


U <- function(pars, y, time, potential) {

	cumcurve <- function(tau){ 
            ratio = exp(alpha) / (exp(alpha) + tau)
            res = 1 - ratio^exp(r)
            res
    }

	r <- pars[1]
	alpha <- pars[2]
	incy <- diff(c(0, y))
	cdf1 <- cumcurve(time)
	cdf2 <- cumcurve(time - 1)
	cdf3 <- cumcurve(max(time))

	ll <- sum(incy * log(cdf1 - cdf2)) + (potential - sum(incy)) * log(1 - cdf3)
	-ll
}

grad_U <- function(pars, y, time, potential) {

	cumcurve <- function(tau){ 
            ratio = exp(alpha) / (exp(alpha) + tau)
            res = 1 - ratio^exp(r)
            res
    }

	gradient.r <- function(tau){ 
            ratio = exp(alpha) / (exp(alpha) + tau)
            res = - log(ratio) * exp(r) * ratio^exp(r)
            res
    }

	gradient.alpha <- function(tau){ 
            ratio = exp(alpha) / (exp(alpha) + tau)
            res = - exp(r) * exp(alpha) * (tau / (exp(alpha) + tau)^2) * ratio^(exp(r) - 1)
            res
    }

	r <- pars[1]
	alpha <- pars[2]
	incy <- diff(c(0, y))

	cdf1 <- cumcurve(time)
	cdf2 <- cumcurve(time - 1)
	cdf3 <- cumcurve(max(time))

#### wrt to r
	gradient1 <- gradient.r(time)
	gradient2 <- gradient.r(time - 1)
	gradient3 <- gradient.r(max(time))

	ll.grad.r <- sum(incy * ( 1 / (cdf1 - cdf2)) * (gradient1 - gradient2) )  +
            (potential - sum(incy)) * (1 / (1 - cdf3)) * ( - gradient3 )

#### wrt to alpha 
	gradient1 <- gradient.alpha(time)
	gradient2 <- gradient.alpha(time - 1)
	gradient3 <- gradient.alpha(max(time))

	ll.grad.alpha <- sum(incy * ( 1 / (cdf1 - cdf2)) * (gradient1 - gradient2) )  + 
            (potential - sum(incy)) * (1 / (1 - cdf3)) * ( - gradient3 )

    c(-ll.grad.r, -ll.grad.alpha)
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

