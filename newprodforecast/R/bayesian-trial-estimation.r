#
#

#' @title Convenience function to create data inputs for Hierarchical Bayesian trial model
#'
#' @description
#' Splits an input data.frame or data.table by a group variable and returns a list with two components ylist and xlist containing the dependent variable and covariates for each group respectively
#' @param data a data.frame or data.table containing the relevant variables
#' @param groupvar is a string that specifies the name of the group variable
#' @param yvar is a string that specifies the name of the dependent variable
#' @parma xvar is a character vector specifying the covariates
#'
#' @return A list with two components ylist and xlist which are lists containing group-level vectors/matrices of the dependent variable and covariates respectively

get_xylist <- function(data, groupvar, yvar, xvars) {
	data <- as.data.table(data)
	y <- split(data[, yvar, with = FALSE], data[, groupvar, with = FALSE])
	x <- split(data[, xvars, with = FALSE], data[, groupvar, with = FALSE])
	x <- lapply(x, as.matrix)
	list(ylist = y, xlist = x)
}


loglik_nls <- function(params, y, x, transforms, sigsq) {
	# Likelihood for the nonlinear regression y = zeta(x, params) + epsilon; epsilon ~ N(0, sigsq)
	# zeta is a nonlinear function of parameters and covariates
	# The call to cumulative_curve computes zeta
	#
	# params is a named list of constrained parameters
	# sigsq is a numeric scalar (NOTE: this is the LOG of sigsq)
	# y is a numeric vector 
	# x is a numeric matrix

	sigsq <- exp(sigsq)
	params <- apply_transforms(params, transforms)
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	tpredictor <- cumulative_curve(params, predictor, acv = 1)	
	sum(dnorm(y, mean = tpredictor, sd = sqrt(sigsq), log = TRUE))
}

logprior_normal <- function(params, paramIx, mu, Sigma) { 
	# Calculate the conditional normal density for a partitioned Multivariate Normal distribution
	# params is a named list of unconstrained parameters which are ~ MVN(mu, Sigma)
	# paramIx is an integer vector of indices of the subvector whose conditional density is required
	# mu is numeric vector
	# Sigma is a numeric Matrix
	x <- unlist(params)
	x2 <- x[paramIx]
	x1 <- x[-paramIx]
	mu2 <- mu[paramIx]
	mu1 <- mu[-paramIx]
	mat1 <- Sigma[paramIx, -paramIx, drop = FALSE] %*% chol2inv(chol(Sigma[-paramIx, -paramIx, drop = FALSE]))
    # Conditional mean and variance
	mu_2_1 <- mu2 - mat1 %*% as.matrix(x1 - mu1)
	Sigma_2_1 <- Sigma[paramIx, paramIx, drop = FALSE] - mat1 %*% Sigma[-paramIx, paramIx, drop = FALSE]
	dmvnorm(x2, mu_2_1, Sigma_2_1, log = TRUE)
}

sigsq_sample <- function(sigsq, params, transforms, y, x, loglik, logprior, rejections, rwscale, sigsq_mean, sigsq_sd) {
	sigsqProp <- sigsq + rwscale * rnorm(1)          # NOTE: sigsq is the log of sigma^2 (error variance)
	loglikProp <- mapply(loglik_nls, params, y, x, MoreArgs = list(sigsq = sigsq, transforms = transforms))
	logpriorProp <- dnorm(sigsqProp, sigsq_mean, sigsq_sd, log = TRUE)
	posteriorRatio <- exp(sum(loglikProp) + logpriorProp - sum(loglik) - logprior)

	msg <- sprintf("ll prop: %f\nlp prop: %f\nll: %f\nlp: %f\n", loglikProp, logpriorProp, loglik, logprior)
	if (is.nan(posteriorRatio) || is.na(posteriorRatio))
		stop(paste(sprintf("%f generated\n\n", posteriorRatio), msg))

	if (posteriorRatio > runif(1)) {
		sigsq <- sigsqProp
		loglik <- loglikProp
		logprior <- logpriorProp
	} else {
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

	if (!is.relistable(params))
		stop("params must be a relistable object")

	paramsProp <- unlist(params)
	paramsProp[paramIx] <- paramsProp[paramIx] + rwscale * rnorm(length(paramIx))
	paramsProp <- relist(paramsProp)
	
	#paramsProp <- params
	#k <- length(params[[paramName]])
	#paramsProp[[paramName]] <- params[[paramName]] + rwscale * rnorm(k)
	loglikProp <- loglik_nls(paramsProp, y, x, transforms, sigsq)
	logpriorProp <- logprior_normal(paramsProp, paramIx, mu, Sigma)
	posteriorRatio <- exp(loglikProp + logpriorProp - loglik - logprior)

	msg <- sprintf("ll prop: %f\nlp prop: %f\nll: %f\nlp: %f\n", loglikProp, logpriorProp, loglik, logprior)
	if (is.nan(posteriorRatio) || is.na(posteriorRatio))
		stop(paste(sprintf("%f generated\n\n", posteriorRatio), msg))
	
	if (posteriorRatio > runif(1)) {
		params <- paramsProp
		loglik <- loglikProp
		logprior <- logpriorProp
	} else {
		rejections <- rejections + 1L
	}
	list(params = params, loglik = loglik, logprior = logprior, rejections = rejections)
}

#' @title Hierarchical Bayesian Estimation of the Trial Model
#' 
#' @description
#' Fits a hierarchical bayesian trial model
#' @param Data a list with components ylist, xlist, and z. ylist and xlist are lists of vectors and matrices respectively, representing the dependent variable and the covariates. Each element of the list represents a group in the lower level of the hierarchical model. z is a matrix of covariates for the upper level of the hierarchy.
#' @param Priors a list with components Deltabar, A, nu, V, sigsq_mean, and sigsq_sd. The components Deltabar, A, nu, and V are parameters of the hyperprior in the multivariate regression for the upper level model (LINK TO bayesm::rmultireg). sigsq_mean and sigsq_sd represent the mean and standard deviation of the univariate normal prior on the log of the variance of the observation error term (sigma^2). 
#' @param Mcmc a list with starting values of model parameters, scaling parameters for the Random Walk Metropolis algorithm, and options controlling various aspects of the MCMC sampler. 
#' The components of Mcmc are as follows:
#' paramsList: a list of lists. Each sublist is a named list of unconstrained parameters 
#' sigsq: is a scalar representing the log of the variance of the error term
#' Delta: is a matrix of parameters for the mean of the multivariate normal prior distribution. (TODO: Add reference to Rossi, Allenby, McCulloch)
#' Sigma: is the covariance matrix for the multivariate normal prior distribution
#' rscale, alphascale, betascale, sigsqscale: are scale parameters for the random walk Metropolis for (unconstrained) parameters r, alpha, beta, and log(sigma^2) respectively
#' burn: is the number of burnin iterations prior to sampling
#' samples: is the number of iterations during the sampling phase
#' thin: controls the frequency at which samples are stored
#' printThin: controls the frequency of printing diagnostic information
#'
#' @return A list of lists with components param_samples, prior_samples, rejections, loglikelihood, logprior. param_samples is a list of arrays/matrices storing samples of parameters from the lower model. prior_samples is a list of arrays/matrices storing parameters from the upper model. loglikelihood is a matrix of log likelihood values. logprior is a list of matrices storing the log priors for lower model parameters. 
#'
#' @export
hb_trialmodel <- function(Data, Priors, Mcmc) {
	numGroups <- length(Data$ylist)
	numPars <- length(unlist(Mcmc$paramsList[[1]]))
	numCovariates <- length(Mcmc$paramsList[[1]]$beta)
	numUpperCovariates <- NCOL(Data$z)
	numIter <- Mcmc$samples + Mcmc$burn
	numSaved <- Mcmc$samples %/% Mcmc$thin

	# Define model specification (Needed for generating parameter transformations)
	tmSpec <- list(family = "exponential-gamma", p0 = FALSE, numCovariates = NCOL(Data$xlist[[1]])) # HARD CODED FOR NOW. NEEDS TO BE AN INPUT
	# Generate functions that are used to transform the parameters while computing likelihood
	paramTransforms <- generate_transforms(tmSpec)
	paramIx <- generate_vec_splitter(tmSpec)(seq(numPars))

	paramsList <- Mcmc$paramsList
	sigsq <- Mcmc$sigsq

	loglik <- mapply(loglik_nls, paramsList, Data$ylist, Data$xlist, 
					 MoreArgs = list(sigsq = Mcmc$sigsq, transforms = paramTransforms))
	mu <- Data$z %*% Mcmc$Delta
	muList <- split(mu, seq(numGroups))
	Sigma <- Mcmc$Sigma
	r_logprior <- mapply(logprior_normal, paramsList, muList,
						 MoreArgs = list(paramIx = paramIx$r, Sigma = Sigma))
	alpha_logprior <- mapply(logprior_normal, paramsList, muList, 
							 MoreArgs = list(paramIx = paramIx$alpha, Sigma = Sigma))
	beta_logprior <- mapply(logprior_normal, paramsList, muList,
							 MoreArgs = list(paramIx = paramIx$beta, Sigma = Sigma))
    sigsq_logprior <- dnorm(sigsq, Priors$sigsq_mean, Priors$sigsq_sd, log = TRUE)

	# Initialize rejection counters
	r_rejections <- alpha_rejections <- beta_rejections <- sigsq_rejections <- rep(0, numGroups)

	# Initialize storage for outputs
	params_samples <- list(r = matrix(NA, numSaved, numGroups),
						   alpha = matrix(NA, numSaved, numGroups),
						   beta = array(NA, c(numSaved, numGroups, numCovariates)),
						   sigsq = rep(NA, numSaved))
	prior_samples <- list(Delta = array(NA, c(numSaved, numUpperCovariates, numPars)),
						  Sigma = array(NA, c(numSaved, numPars, numPars)))
	rejections <- list(r = matrix(NA, numSaved, numGroups),
					   alpha = matrix(NA, numSaved, numGroups),
					   beta = matrix(NA, numSaved, numGroups))
	loglikelihood <- matrix(NA, numSaved, numGroups)
	logprior <- list(r = matrix(NA, numSaved, numGroups),
					 alpha = matrix(NA, numSaved, numGroups),
					 beta = matrix(NA, numSaved, numGroups)) # Initialize storage for log prior

	save_iter <- 0L
    start_time <- proc.time()[3]
	cat("Starting burn-in phase ...", fill = TRUE)
	for (iter in seq(numIter)) {		
		# Sample r
		r_mhres <- mapply(metropolis_sample, paramsList, Data$ylist, Data$xlist, loglik, r_logprior, r_rejections, muList,
						  MoreArgs = list(paramName = "r", transforms = paramTransforms, paramIx = paramIx$r, sigsq = sigsq, 
                                          rwscale = Mcmc$rscale, Sigma = Sigma), SIMPLIFY = FALSE)
		paramsList <- lapply(r_mhres, "[[", i = "params")
		loglik <- sapply(r_mhres, "[[", i = "loglik")
		r_logprior <- sapply(r_mhres, "[[", i = "logprior")
		r_rejections <- sapply(r_mhres, "[[", i = "rejections")

		# Sample alpha
		alpha_mhres <- mapply(metropolis_sample, paramsList, Data$ylist, Data$xlist, loglik, alpha_logprior, alpha_rejections, muList,
							  MoreArgs = list(paramName = "alpha", transforms = paramTransforms, paramIx = paramIx$alpha, sigsq = sigsq, 
                                              rwscale = Mcmc$alphascale, Sigma = Sigma), SIMPLIFY = FALSE)
		paramsList <- lapply(alpha_mhres, "[[", i = "params")
        loglik <- sapply(alpha_mhres, "[[", i = "loglik")
		alpha_logprior <- sapply(alpha_mhres, "[[", i = "logprior")
		alpha_rejections <- sapply(alpha_mhres, "[[", i = "rejections")

		# Sample beta
		beta_mhres <- mapply(metropolis_sample, paramsList, Data$ylist, Data$xlist, loglik, beta_logprior, beta_rejections, muList,
							 MoreArgs = list(paramName = "beta", transforms = paramTransforms, paramIx = paramIx$beta, sigsq = sigsq, 
                                             rwscale = Mcmc$betascale, Sigma = Sigma), SIMPLIFY = FALSE)
		paramsList <- lapply(beta_mhres, "[[", i = "params")
		loglik <- sapply(beta_mhres, "[[", i = "loglik")
		beta_logprior <- sapply(beta_mhres, "[[", i = "logprior")
		beta_rejections <- sapply(beta_mhres, "[[", i = "rejections")

		# Sample sigsq
		sigsq_mhres <- sigsq_sample(sigsq, paramsList, paramTransforms, Data$ylist, Data$xlist, loglik, sigsq_logprior, 
                                    sigsq_rejections, Mcmc$sigsqscale, Priors$sigsq_mean, Priors$sigsq_sd)
        sigsq <- sigsq_mhres$sigsq
		loglik <- sigsq_mhres$loglik
		sigsq_logprior <- sigsq_mhres$logprior
		sigsq_rejections <- sigsq_mhres$rejections

		# Sample upper model parameters (Delta, Sigma)
		paramsMatrix <- do.call("rbind", lapply(paramsList, unlist))
		upper_res <- rmultireg(paramsMatrix, Data$z, Priors$Deltabar, Priors$A, Priors$nu, Priors$V)
		Delta <- upper_res$B
		mu <- Data$z %*% Delta
		muList <- split(mu, seq(numGroups))
		Sigma <- upper_res$Sigma

		# Print diagnostics
		if (iter %% Mcmc$printThin == 0) {            
            now <- proc.time()[3]
            sumfuns <- list(min, median, mean, max)
            fmt <- paste("%6s", paste(rep("%6.2f%%", length(sumfuns)), collapse = " "))
            #print(fmt)            
			cat("Iterations:", iter, round(now - start_time), " secs", fill = TRUE)
			cat("Rejection Rates:", fill = TRUE)
            cat(sprintf("%12s %7s %7s %7s", "Min", "Med", "Mean", "Max"), fill = TRUE)
            msg <- do.call(sprintf, c(fmt, "r", lapply(sumfuns, function(f, x) f(x) * 100 / numIter, x = r_rejections)))
			cat(msg, fill = TRUE)
            msg <- do.call(sprintf, c(fmt, "alpha", lapply(sumfuns, function(f, x) f(x) * 100 / numIter, x = alpha_rejections)))
            cat(msg, fill = TRUE)
            msg <- do.call(sprintf, c(fmt, "beta", lapply(sumfuns, function(f, x) f(x) * 100 / numIter, x = beta_rejections)))
            cat(msg, fill = TRUE)
            msg <- do.call(sprintf, c(fmt, "sigsq", lapply(sumfuns, function(f, x) f(x) * 100 / numIter, x = sigsq_rejections)))
            cat(msg, fill = TRUE)       
			linesep()
		}
		# Save
		if (iter == (Mcmc$burn + 1L)) {
            now <- proc.time()[3]
            cat("Burn-in done in ", round(now - start_time), " secs. Starting sampling ...", fill = TRUE)
		}
		if (iter > Mcmc$burn && (iter - Mcmc$burn) %% Mcmc$thin == 0L) {
			save_iter <- save_iter + 1L
			params_samples$r[save_iter,] <- sapply(paramsList, "[[", i = "r")
			params_samples$alpha[save_iter,] <- sapply(paramsList, "[[", i = "alpha")
			params_samples$beta[save_iter,,] <- do.call("rbind", lapply(paramsList, "[[", i = "beta"))
            params_samples$sigsq[save_iter] <- sigsq
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
	list(params_samples = params_samples, 
         prior_samples = prior_samples, 
         rejections = rejections, 
         loglikelihood = loglikelihood, 
         logprior = logprior)
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

