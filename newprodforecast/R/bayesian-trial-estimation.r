#
#
hb_trialmodel <- function(formula, 
						  data, 
						  group,
						  startvals = numeric(),
						  tmSpec = list(),
						  priors = list(),
						  mcmcControl = list()) {

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	cl$formula <- mf$formula <- autocorrectModelFormula(formula)

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

	getXY <- function(dat, mfcall) {
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

	datalst <- dat[, getXY(.SD, mf), by = group]

	sampleSize <- length(datlst)
	burnin <- mcmcControl$burnin
	samples <- mcmcControl$samples
	numIter <- burnin + samples
	L <- mcmcControl$L
	epsilon <- mcmcControl$epsilon

	if (tmSpec$family == "exponential") {
		for (iter in seq(samples)) {
			#lambda <- mapply(datalst, params, 	
      lambda <- sampleLambda()
			# sample beta
			beta <- sampleBeta()
			# sample Delta, Sigma
		}
	} else {
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


metropolis_sample <- function(params, xylist, loglik, logprior, rejections, mu, Sigma, rwscale, paramName, paramIx, tmSpec) {
	paramsProp <- params
	k <- length(params[[paramName]])
	paramsProp[[paramName]] <- params[[paramName]] + rwscale * rnorm(k)
	loglikProp <- logloglik_trialmodel(paramsProp, xylist$y, xylist$x, acv = 1, tmSpec)
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

