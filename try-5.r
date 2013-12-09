require(data.table)
source("utils.r")
	
mle_allids <- function(lst, useSF = FALSE, parallel = TRUE, cpus = 4, ...) {
	if (useSF) 
		require(snowfall)
	# Wrapper to apply MLE to all products
	# lst : list of data.tables, each element corresponds to a product
	# useSF : flag to use snowfall
	# parallel, cpus : arguments to be passed to snowfall
	# ... : arguments to be passed to mle_oneid
	call <- match.call()
	if (useSF) {
		sfInit(parallel = parallel, cpus = cpus, type = 'SOCK')
		sfExport("mle_oneid", "objfn", "cdf", "loglikelihood", "%notin%", "invlogit")
		res <- sfLapply(lst, mle_oneid, ...)
		sfStop()
	} else {
		res <- lapply(lst, mle_oneid, ...)
	}
	structure(res, call = call)
}

mle_oneid <- function(dat, yvar, xvars = NULL, ACV = NULL, p0 = FALSE, startvals = NULL, ...) {
	# Perform maximum likelihood estimation for one product 
	# dat : data.table containing data for model estimation
	# yvar : string with name of dependent variable in dat
	# xvars : vector of strings with names of covariates. Default value of NULL indicates no covariates
	# ACV : string with name of variable containing monotone transform of distribution. Default value of NULL indicates that this term is not included in the model
	# p0 : logical. If TRUE it indicates that the asymptote parameter is to be included in the model
	# startvals : numeric vector of starting values for parameters. Default value is NULL
	# Value : list produced by optim. If optim fails returns a try-error object
	call <- match.call()
	npars <- 2L + p0 + length(xvars) + !is.null(ACV)
	if (is.null(startvals))
		startvals <- rep(-1, npars)
	else if (length(startvals) != npars)
		stop("Incorrect number of starting parameters")

	y <- dat[[yvar]]
	if (is.null(xvars))
		X <- matrix(seq_along(y), ncol = 1)
	else 
		X <- as.matrix(dat[, xvars, with = FALSE])

	if (!is.null(ACV)) 
		ACV <- dat[[ACV]]

	res <- try(optim(startvals, objfn, y = y, X = X, hasP0 = p0, ACV = ACV, hasCov = !is.null(xvars),
					 control = list(fnscale = -1, maxit = 20000), method = "BFGS"))
	if (!inherits(res, "try-error")) res <- structure(res, call = call)
	res
}	

objfn <- function(params, y, X, hasP0, ACV, hasCov) {
	# Convenience function to handle different model specifications before calling likelihood function
	# params : numeric vector of unconstrained parameters
	# y : vector with dependent variable
	# X : matrix with covariates. If no covariates present it is a single column with time periods
	# hasP0 : logical. If TRUE the asymptote parameter is included in the model
	# ACV : numeric vector with monotone transformation of distribution. If this term is not included in the model then it is NULL
	# hasCov : logical. If TRUE then covariates are included in the model
	# Value : numeric value of likelihood
	boffset <- 2L                         # offset from start of param vector to position of first beta (covariate effect)
	nms <- c("r", "alpha")                # minimal vector of param names

	if (hasP0) {
		params[1] <- invlogit(params[1])     # Inverse logit transform for p0 parameter
		params[-1] <- exp(params[-1])        # Exponential transform for other parameters
		boffset <- boffset + 1L
		nms <- c("p0", nms)
	} else {
	   params <- exp(params)
	}	
	if (!is.null(ACV)) {
		boffset <- boffset + 1L
		if (hasP0) 
			nms <- c(nms[1], "gamma", nms[-1]) 
		else 
			nms <- c("gamma", nms)
	}
	if (hasCov) {
		params <- c(as.list(params[seq(boffset)]), list(params[-seq(boffset)]))
		nms <- c(nms, "beta")
	} else {
		params <- as.list(params)
	}
	names(params) <- nms
	ll <- loglikelihood(params, y, X, ACV, 1)	
	ll
}

loglikelihood <- function(pars, y, X, ACV, potential) {
	incy <- diff(c(0, y))
	cdf1 <- cdf(pars, X, ACV)
	cdf2 <- cdf(pars, X, ACV, "lag")
	cdf3 <- cdf(pars, X, ACV, "max")
	ll <- sum(incy * log(cdf1 - cdf2)) + (potential - sum(incy)) * log(1 - cdf3)
	ll
}

cdf <- function(pars, X, ACV, mod = "") {
	# Calculate the trial curve as a function of parameters and covariates
	# pars : named list of parameters
	# X : matrix of covariate
	# ACV : numeric vector with montone transform of distribution
	# mod : string representing one of three scenarios. Default is "". Other scenarios are "lag" and "max"
	r <- pars$r
	alpha <- pars$alpha
	# Set defaults for optional parameters
	if (is.null(pars$p0)) p0 <- 1 else p0 <- pars$p0
	if (is.null(pars$beta)) beta <- matrix(0, 1, 1) else beta <- pars$beta
	if (is.null(pars$gamma)) gamma <- 0 else gamma <- pars$gamma
	if (is.null(ACV)) ACV <- 1

	A <- exp(X %*% beta)
	A <- switch(mod,
				"lag" = cumsum(c(0, head(A, -1))),
				"max" = tail(cumsum(A), 1),
				cumsum(A))
	ACV <- switch(mod, 
				  "lag" = 1,
				  "max" = max(ACV),
				  ACV)
	p0 * ((ACV/100) ^ gamma) * (1 - (alpha / (alpha + A))^r)
}

get_params <- function(mle_res) {
	if (inherits(mle_res, "try-error")) {
		return(list())
	} else {
		cl <- as.list(attr(mle_res, "call"))
		params <- mle_res$par
		nms <- c("r", "alpha")
		boffset <- 2L
		hasP0 <- ("p0" %in% names(cl)) && cl$p0
		hasCov <- ("xvars" %in% names(cl)) && !is.null(cl$xvars)

		# Check the call to see if p0 is 
		if (hasP0) {
			params[1] <- invlogit(params[1])
			params[-1] <- exp(params[-1])
			boffset <- boffset + 1L
			nms <- c("p0", nms)		
		} else {
			params <- exp(params)			
		}
		if ("ACV" %in% names(cl) && !is.null(cl$ACV)) {
			boffset <- boffset + 1L
			if (hasP0)
				nms <- c(nms[1], "gamma", nms[-1])
			else
				nms <- c("gamma", nms)
		}
		if (hasCov) {
			params <- c(as.list(params[seq(boffset)]), list(params[-seq(boffset)]))
			nms <- c(nms, "beta")
		} else {
			params <- as.list(params)
		}
		names(params) <- nms
		return(params)
	}
}



fit_allids <- function(mle_res, obs_dat, ...) {
	cl <- as.list(attr(mle_res, "call"))
	mlst <- list(yvar = cl$yvar, xvars = cl$xvars, ACV = cl$ACV, p0 = ifelse(is.null(cl$p0), FALSE, cl$p0))
	res <- mapply(fit_oneid, mle_res, obs_dat, MoreArgs = mlst, SIMPLIFY = FALSE)
	res <- rbindlist(res)
	structure(res, call = cl)
}


fit_oneid <- function(mle_res, obs_dat, yvar, xvars = NULL, ACV = NULL, p0 = FALSE) {
	est_params <- get_params(mle_res)
	if (length(est_params)) {
		y <- obs_dat[[yvar]]
		if (is.null(xvars))
			X <- matrix(seq_along(y), ncol = 1)
		else
			X <- as.matrix(obs_dat[, xvars, with = FALSE])
		if (!is.null(ACV))
			ACV <- obs_dat[[ACV]]

		fit <- cdf(est_params, X, ACV)		
		mape <- median(100 * abs(fit / (y + 1e-6) - 1))
		res <- cbind(obs_dat, PredTrialPct = fit, MAPE = mape)
		if (!is.null(est_params$p0))
			res$p0 <- est_params$p0
		if (!is.null(est_params$gamma))
			res$gamma <- est_params$gamma
		res$r <- est_params$r
		res$alpha <- est_params$alpha
		if (!is.null(est_params$beta)) {
			betamat <- tcrossprod(rep(1, nrow(obs_dat)), est_params$beta)
			colnames(betamat) <- xvars
			res <- cbind(res, beta = betamat)
		}
		return(res)
	} else {
		return(data.table())
	}
}

