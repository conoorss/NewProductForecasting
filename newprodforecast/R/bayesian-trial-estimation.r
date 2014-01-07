#
#

loglik_nls <- function(params, tranforms, sigsq, y, x) {
	# Likelihood for the nonlinear regression y = zeta(x, params) + epsilon; epsilon ~ N(0, sigsq)
	# zeta is a nonlinear function of parameters and covariates
	# The call to cumulative_curve computes zeta
	#
	# params is a named list of constrained parameters
	# sigsq is a numeric scalar (NOTE: this is the LOG of sigsq)
	# y is a numeric vector 
	# x is a numeric matrix

	sigsq <- exp(sigsq)
	params <- mapply(function(x, fun) fun(x), params, transforms)
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	tpredictor <- cumulative_curve(params, predictor, acv = 1)	
	-sum((y - tpredictor)^2) / (2 * exp(sigsq))
	dnorm(y, mean = tpredictor, sd = sqrt(sigsq), log = TRUE)
}

logprior_normal <- function(params, paramName, paramIx, mu, Sigma) { 
	# Calculate the conditional normal density for the component given by paramName
	# params is a named list of unconstrained parameters
	# paramName is a string
	# paramIx is an integer vector representing the indices of the component paramName in the atomic version of param (i.e., unlist(param))
	# mu is numeric vector
	# Sigma is a numeric Matrix
	x <- params[[paramName]]
	mat1 <- Sigma[paramIx, -paramIx, drop = FALSE] %*% chol2inv(chol(Sigma[-paramIx, -paramIx, drop = FALSE]))
	mu_x <- mu[paramIx] - mat1 %*% (x - mu[-paramIx])
	Sigma_x <- Sigma[paramIx, paramIx, drop = FALSE] - mat1 %*% Sigma[-paramIx, paramIx, drop = FALSE]
	dmvnorm(x, mu_x, Sigma_x, log = TRUE)	
}

sigsq_sample <- function(sigsq, params, transforms, y, x, loglik, logprior, rejections, rwscale, sigsq_mean, sigsq_sd) {
	sigsqProp <- sigsq + rwscale * rnorm(1)          # NOTE: sigsq is the log of sigma^2 (error variance)
	loglikProp <- loglik_nls(params, transforms, sigsqProp, y, x)
	logpriorProp <- dnorm(sigsqProp, sigsq_mean, sigsq_sd, log = TRUE)
	posteriorRatio <- exp(sum(loglikProp) + logpriorProp - sum(loglik) - logprior)

	msg <- sprintf("ll prop: %f\nlp prop: %f\nll: %f\nlp: %f\n", loglikProp, logpriorProp, loglik, logprior)
	if (is.nan(posteriorRatio) || is.na(posteriorRatio))
		stop(paste(sprintf("%f generated\n\n", posteriorRatio), msg))

	if (posteriorRatio > runif(1)) {
		sigsq <- sigsqProp
		loglik <- loglikProp
		logprior <- logpriorProp
		rejections <- rejections + 1L
	}
	list(sigsq = sigsq, loglik = loglik, logprior = logprior, rejections = rejections)
}

metropolis_sample <- function(params, paramName, paramIx, transforms, y, x, loglik, logprior, rejections, rwscale, sigsq, mu, Sigma) {
	# Metropolis Hastings update to the component paramName of the list params
	#
	# params is a named list of unconstrained parameters
	# paramName is a string 
	# paramIx is an integer vector representing the indices of the component paramName in the atomic version of param (i.e., unlist(param))
	# y is a numeric vector
	# x is a numeric matrix
	# loglik is a numeric scalar
	# logprior is a numeric scalar
	# rejections is an integer scalar
	# sigsq is a numeric scalar
	# mu is a numeric vector
	# Sigma is a numeric matrix

	paramsProp <- params
	k <- length(params[[paramName]])
	paramsProp[[paramName]] <- params[[paramName]] + rwscale * rnorm(k)
	loglikProp <- loglik_nls(paramsProp, transforms, sigsq, y, x)
	logpriorProp <- logprior_normal(params, paramName, paramIx, mu, Sigma)
	posteriorRatio <- exp(loglikProp + logpriorProp - loglik - logprior)

	msg <- sprintf("ll prop: %f\nlp prop: %f\nll: %f\nlp: %f\n", loglikProp, logpriorProp, loglik, logprior)
	if (is.nan(posteriorRatio) || is.na(posteriorRatio))
		stop(paste(sprintf("%f generated\n\n", posteriorRatio), msg))
	
	if (posteriorRatio > runif(1)) {
		params <- paramsProp
		loglik <- loglikProp
		logprior <- logpriorProp
		rejections <- rejections + 1L
	}
	list(params = params, loglik = loglik, logprior = logprior, rejections = rejections)
}

#' @title Hierarchical Bayesian Estimation of the Trial Model
#' 
#' @description
#' Fits a hierarchical bayesian trial model
#' @param Data is a list with components yList, xList, and z. yList and xList are lists of vectors and matrices respectively, representing the dependent variable and the covariates. Each element of the list represents a group in the lower level of the hierarchical model. z is a matrix of covariates for the upper level of the hierarchy.
#' @param Priors is a list with components Deltabar, A, nu, V, sigsq_mean, and sigsq_sd. The components Deltabar, A, nu, and V are parameters of the hyperprior in the multivariate regression for the upper level model (LINK TO bayesm::rmultireg). sigsq_mean and sigsq_sd represent the mean and standard deviation of the univariate normal prior on the log of the variance of the observation error term (sigma^2). 
#' @param Mcmc is a list with starting values of model parameters, scaling parameters for the Random Walk Metropolis algorithm, and options controlling various aspects of the MCMC sampler. The following are components of the Mcmc list: 
#' @return A list of lists with components param_samples, prior_samples, rejections, loglikelihood, logprior. param_samples is a list of arrays/matrices storing samples of parameters from the lower model. prior_samples is a list of arrays/matrices storing parameters from the upper model. loglikelihood is a matrix of log likelihood values. logprior is a list of matrices storing the log priors for lower model parameters. 
#'
#'
#' @export
hb_trialmodel <- function(Data, Priors, Mcmc) {
	# 
	# Data is a list with the following components
	# 	yList list of dependent variable vector
	# 	xList list of covariate matrices
	# 	z is matrix of upper model covariates
	#
	# Priors is a list of prior parameters
	#	Deltabar prior parameter for Delta
	#	A prior parameter for Delta
	#	nu prior parameter for Sigma
	#	V prior parameter for Sigma
	#	sigsq_mean prior parameter for sigsq
	#	sigsq_sd parameter for sigsq
	#
	# Mcmc is a list of Mcmc parameters
	# 	The following elements are starting values
	#		paramsList list of named lists. Each named list contains group-level unconstrained parameters
	#		sigsq is a scalar (NOTE: this is the log of error variance)
	#		Delta is a matrix
	#		Sigma is a matrix
	#	The following elements are random walk scaling parameters
	#		rscale
	#		alphascale
	#		betascale	
	#		sigsqscale
	#	The following elements are options for length of run, printing and saving
	#		burn
	#		samples
	#		thin
	#		printThin
	#
	#
	# Return value is a list with the following components
	#	params_samples is a list with components:
	#		r 		- matrix of samples x groups
	#		alpha 	- matrix of samples x groups
	#		beta 	- array of samples x groups x covariates			
	#		sigsq	- vector
	#	NOTE: ALL PARAMETERS ARE REPORTED UNTRANSFORMED
	#	prior_samples is a list with components:
	#		Delta 	- array of samples x upper covariates x lower parameters
	#		Sigma 	- array of samples x upper covariates x upper covariates
	#	rejections is a list with components
	#		r 		- matrix of samples x groups
	#		alpha 	- matrix of samples x groups
	#		beta 	- matrix of samples x groups
	#	loglikelihood is matrix of samples x groups
	#	logprior is a list with components
	#		r 		- matrix of samples x groups
	#		alpha	- matrix of samples x groups
	#		beta	- matrix of samples x groups

	numGroups <- length(Data$yList)
	numPars <- length(unlist(Mcmc$paramList[[1]]))
	numCovariates <- length(Mcmc$paramList[[1]]$beta)
	numUpperCovariates <- NCOL(Data$z)
	numIter <- Mcmc$samples + Mcmc$burn
	numSaved <- Mcmc$samples %/% Mcmc$thin

	# Define model specification (Needed for generating parameter transformations)
	tmSpec <- list(family = "exponential-gamma", p0 = FALSE, numCovariates = NCOL(Data$xList[[1]])) # HARD CODED FOR NOW. NEEDS TO BE AN INPUT
	# Generate functions that are used to transform the parameters while computing likelihood
	paramTransforms <- generate_transforms(tmSpec)
	paramIx <- generate_vec_splitter(tmSpec)(seq(numPars))

	paramsList <- Mcmc$paramsList
	sigsq <- Mcmc$sigsq

	loglik <- mapply(loglik_nls, paramsList, Data$yList, Data$xList, MoreArgs = list(sigsq = Mcmc$sigsq, transforms = paramTransforms))
	mu <- Data$z %*% Mcmc$Delta
	Sigma <- Mcmc$Sigma
	r_logprior <- lapply(paramsList, logprior_normal, paramName = "r", paramIx = paramIx$r, mu = mu, Sigma = Mcmc$Sigma)
	alpha_logprior <- lapply(paramsList, logprior_normal, paramName = "alpha", paramIx = paramIx$alpha, mu = mu, Sigma = Sigma)
	beta_logprior <- lapply(paramsList, logprior_normal, paramName = "beta", paramIx = paramIx$beta, mu = mu, Sigma = Sigma)

	# Initialize rejection counters
	r_rejections <- alpha_rejections <- beta_rejections <- rep(0, numGroups)

	# Initialize storage for outputs
	params_samples <- list(r = matrix(NA, numSaved, numGroups),
						   alpha = matrix(NA, numSaved, numGroups),
						   beta = array(NA, c(numSaved, numGroups, numCovariates)),
						   sigsq = rep(NA, numSaved))
	prior_samples <- list(Delta = array(NA, c(numSaved, numUpperCovariates, numPars)),
						  Sigma = array(NA, c(numSaved, numUpperCovariates, numUpperCovariates)))
	rejections <- list(r = matrix(NA, numSaved, numGroups),
					   alpha = matrix(NA, numSaved, numGroups),
					   beta = matrix(NA, numSaved, numGroups))
	loglikelihood <- matrix(NA, numSaved, numGroups)
	logprior <- list(r = matrix(NA, numSaved, numGroups),
					 alpha = matrix(NA, numSaved, numGroups),
					 beta = matrix(NA, numSaved, numGroups)) # Initialize storage for log prior

	save_iter <- 0L
	for (iter in seq(numIter)) {
		cat("Starting burn-in phase ...", fill = TRUE)
		# Sample r
		r_mhres <- mapply(metropolis_sample, paramsList, Data$yList, Data$xList, loglik, r_logprior, r_rejections,
						  MoreArgs = list(paramName = "r", paramIx = paramIx$r, sigsq = sigsq, rwscale = Mcmc$rscale, mu = mu, Sigma = Sigma))
		paramsList <- lapply(r_mhres, "[[", i = "params")
		loglik <- sapply(r_mhres, "[[", i = "loglik")
		r_logprior <- sapply(r_mhres, "[[", i = "logprior")
		r_rejections <- sapply(r_mhres, "[[", i = "rejections")

		# Sample alpha
		alpha_mhres <- mapply(metropolis_sample, paramsList, yList, xList, loglik, alpha_logprior, alpha_rejections,
							  MoreArgs = list(paramName = "alpha", paramIx = paramIx$alpha, sigsq = sigsq, rwscale = Mcmc$alphascale, mu = mu, Sigma = Sigma))
		loglik <- sapply(alpha_mhres, "[[", i = "loglik")
		alpha_logprior <- sapply(alpha_mhres, "[[", i = "logprior")
		alpha_rejections <- sapply(alpha_mhres, "[[", i = "rejections")

		# Sample beta
		beta_mhres <- mapply(metropolis_sample, paramsList, yList, xList, loglik, beta_logprior, beta_rejections, 
							 MoreArgs = list(paramName = "beta", paramIx = paramIx$beta, sigsq = sigsq, rwscale = Mcmc$betascale, mu = mu, Sigma = Sigma))
		loglik <- sapply(beta_mhres, "[[", i = "loglik")
		beta_logprior <- sapply(beta_mhres, "[[", i = "logprior")
		beta_rejections <- sapply(beta_mhres, "[[", i = "rejections")

		# Sample sigsq
		sigsq_mhres <- sigsq_sample(sigsq, )
		loglik <- do.call("c", sigsq_mhres$loglik)
		sigsq_logprior <- do.call("c", sigsq_mhres$logprior)
		sigsq_rejections <- do.call("c", sigsq_mhres$rejections)

		# Sample upper model parameters (Delta, Sigma)
		paramsMatrix <- do.call("rbind", unlist(paramsList))
		upper_res <- rmultireg(paramsMatrix, z, Priors$Deltabar, Priors$A, Priors$nu, Priors$V)
		Delta <- upper_res$B
		mu <- z %*% Delta
		Sigma <- upper_res$Sigma

		# Print diagnostics
		if (iter %% Mcmc$printThin == 0) {
			cat("Iterations:\t, iter, \n")
			cat("Rejection Rates:\n")
			rrfmt <- "%6s %.2f%%"
			msg <- sprintf(paste(rep(rrfmt, 3), collapse = "\n"), 
						   "r", mean(100 * r_rejections / numIter), 
						   "alpha", mean(100 * alpha_rejections / numIter), 
						   "beta", mean(100 * beta_rejections / numIter))
			cat(msg, fill = TRUE)
			linesep()
		}
		# Save
		if (iter > burn) cat("Burn-in Done. Starting sampling ...", fill = TRUE)
		if (iter > burn && (iter - burn) %% thin == 0L) {
			save_iter <- save_iter + 1L
			params_samples$r[save_iter,] <- sapply(paramsList, "[[", i = "r")
			params_samples$alpha[save_iter,] <- sapply(paramsList, "[[", i = "alpha")
			params_samples$beta[save_iter,,] <- do.call("rbind", lapply(paramsList, "[[", i = "beta"))
			prior_samples$Delta[save_iter,,] <- Delta
			prior_samples$Sigma[save_iter,,] <- Sigma
			rejections$r[save_iter,] <- r_rejections
			rejections$alpha[save_iter,] <- alpha_rejections
			rejections$beta[save_iter,] <- beta_rejections
			loglikelihood[save_iter,] <- loglik
			logprior$r[save_iter,] <- r_logprior
			logprior$alpha[save_iter,] <- alpha_logprior
			logprior$beta[save_iter,] <- beta_logprior
		}
	}
	list(params_samples = params_samples, prior_samples = prior_samples, rejections = rejections, loglikelihood = loglikelihood, logprior = logprior)
}

# All code below this is deprecated for now

if (FALSE) {
	extract_xy <- function(dat, mfcall) {
		m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
		mfcall <- mfcall[c(1L, m)]
		mfcall[[1L]] <- as.name("model.frame")
		mfcall$data <- quote(dat)
		mf <- as.data.table(eval(mfcall))
		mt <- attr(mf, "terms")
		y <- model.response(mf, "numeric")
		x <- model.matrix(mt, mf)
		list(y = y, x = x)
	}

	hb_trialmodel <- function(formula, 
							  data, 
							  group,
							  startvals = numeric(),
							  tmSpec = list(),
							  priors = list(),
							  mcmcControl = list()) {

		cl <- match.call()
		mf <- match.call(expand.dots = FALSE)
		cl$formula <- mf$formula <- autocorrect_formula(formula)

		if (missing(formula) || !inherits(formula, "formula") || length(formula) < 3L) # Check for two sided formula
			stop("Require a formula argument")
		if (missing(data) || !inherits(data, "data.frame"))   # Check if data arg is a data frame
			stop("Require a data.frame or data.table")
		if (!inherits(data, "data.table"))   # Coerce to data.table if needed
			data <- as.data.table(data)
		if (!missing(group) && !is.character(group))
			stop("The group argument must be a string giving the name of a column in the data argument")

		# Specification of the functional form for trial model
		tmSpec <- set_spec_trialmodel(tmSpec)
		tmSpec$numCovariates <- length(attr(terms(mf$formula), "term.labels"))
		if (!is.null(tmSpec$acvMultiplier) && tmSpec$acvMultiplier %notin% names(data))
			stop("The acv multiplier variable in the model specification was not found in the data")


		dataList <- dat[, extract_xy(.SD, mf), by = group]

		numGroups <- length(datlst)
		#L <- mcmcControl$L
		#epsilon <- mcmcControl$epsilon

		loglik <- mapply(loglik_trialmodel, paramList, , )
		numIter <- mcmcControl$samples + mcmcControl$burn

		for (iter in seq(numIter)) {
			for (i in seq(numGroups)) {
				paramList[[i]]$lambda <- metropolis_sample(paramList[[i]], dataList[[i]], loglik[i], logprior[i], rejections[i], mu, Sigma, rwscale, "lambda", 1L)
				z <- extract_lambda(paramList)
				uppermodel <- prior_sample(z, mu, Sigma)
				mu <- uppermodel$mu
				Sigma <- uppermodel$Sigma
				if (i > mcmcControl$burn && i %% mcmcControl$thin) {
					# Store draws in appropriate container
				}
			}
		}
	}

	sample_lambda <- function(params, xylist, loglik, logprior, rejections, mu, Sigma, rwscale, tmSpec){
		# Propose a new value of lambda
		paramsProp <- params
		paramsProp$lambda <- params$lambda + rwscale * rnorm(1)
		loglikProp <- loglik_trialmodel(paramsProp, xylist$y, xylist$x, acv = 1, tmSpec)
		logpriorProp <- calcPrior(paramsProp$lambda, mu, sigma)
		mhRatio <- exp(loglikProp + logpriorProp - loglik - logprior)
		if (is.na(mhRatio) || mhratio > runif(1))
			list(params = params, loglik = loglikProp, logprior = logpriorProp, rejections = rejections)
		else
			list(params = params, loglik = loglik, logprior = logprior, rejections = rejections + 1)
	}


	metropolis_sample <- function(params, xylist, loglik, logprior, rejections, mu, Sigma, rwscale, paramName, paramIx) {
		#
		# params is a named list of parameters
		# xylist is a named lists with a vector y and matrix x
		# loglik is a scalar
		# logprior is a scalar
		# rejections is a scalar
		# mu is a numeric vector representing the mean of the Normal prior  
		# Sigma is a matrix representing the covariance of the Normal prior
		# rwscale is a scalar
		# paramName is a string
		# paramIx is an integer vector
		#
		paramsProp <- params
		k <- length(params[[paramName]])
		paramsProp[[paramName]] <- params[[paramName]] + rwscale * rnorm(k)
		loglikProp <- logloglik_trialmodel(paramsProp, xylist$y, xylist$x, acv = 1)
		# This section needs testing for lower dimension cases and compatibility with dmvnorm
		# Need to include dependency for dmvnorm function
		if (length(mu) == 1L) {
			condMean <- mu
			condVar <- Sigma
		} else {
			condPars <- unlist(paramsProp)[-paramIx]
			condMean <- mu[paramIx] - Sigma[paramIx, -paramIx] %*% chol2inv(chol(Sigma[-paramIx, -paramIx])) %*% (condPars - mu[-paramIx])
			condVar <- Sigma[-paramIx, -paramIx] - Sigma[paramIx, -paramIx] %*% chol2inv(chol(Sigma[-paramIx, -paramIx])) %*% Sigma[-paramIx, paramIx]
		}
		logpriorProp <- logPrior(param, condMean, condVar)
		mhRatio <- exp(loglikProp + logpriorProp - loglik - logprior)
		if (is.na(mhRatio) || mhRatio > runif(1)) {
			params <- paramsProp
			loglik <- loglikProp
			logprior <- logpriorProp
		}
		list(params = params, loglik = loglik, loglik = loglik)
	}

	logPrior <- function(x, mu, Sigma) dmvnorm(x, mu, Sigma, log = TRUE)

}

#		cat("Current q is: ", current_q, "\n")
#		q <- curren_q
#		p <- sig %*% rnorm(length(q))
#		current_p <- p
#
#		# Make a half step for momentum in the beginning
#		gg <- grad_U(q, y, time, potential)
#		p <- p - epsilon * gg / 2
#
#		# Alternate full steps for position and momentum
#		for (i in 1:L) {
#			q <- q + epsilon * p
#			gg <- grad_U(q, y, time, potential)
#			if (i != L) 
#				p <- p - epsilon * gg			
#		}
#		gg <- grad_U(q, y, time, potential)
#		p <- p - epsilon * gg / 2
#		# Negate momentum at the end of the trajectory
#		p <- -p
#		current_U <- U(current_q, y, time, potential)
#		current_K <- sum(current_p^2) / 2
#		proposed_U <- U(q, y, time, potential)
#		proposed_K <- sum(p^2) / 2
#		cat(current_U, current_K, proposed_U, proposed_K, "\n")
#	}
#

