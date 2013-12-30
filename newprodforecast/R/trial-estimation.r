#' @title Specification for Trial Model
#'
#' @description
#' Used to specify options for the trial model:
#' functional family of the cumulative trial curve
#' presence of asymptote parameter
#' presence of multiplicative term for ACV (distribution)
#'
#' @param family string with value "exponential" or "exponential-gamma". Default value is "exponential-gamma"
#' @param p0 logical. If the value is TRUE then the model includes an asymptote parameter p0. Default value is FALSE.
#' @param acv_multiplier string giving the column name of the variable which is to be used in the multiplicative term for distribution.
#' @return A list with named elements "family", "p0", "acv_multiplier"
#' 
#' @export
tm.spec <- function(family = "exponential-gamma", p0 = FALSE, acv_multiplier = NULL, numCovariates = 0L) {
	if (family %notin% c("exponential", "exponential-gamma"))
		warning(gettextf("family must be either %s or %s. Setting it to exponential-gamma.", dQuote("exponential"), dQuote("exponential-gamma")))
	if (!is.null(acv_multiplier) && !is.character(acv_multiplier)) {
		stop("acv_multiplier must be a string")		
	} else if (length(acv_multiplier) > 1L) {
		warning("acv_multiplier must a single string (not a vector). Ignoring all but the first element of the vector")
		acv_multiplier <- acv_multiplier[[1L]]
	}
	list(family = family, p0 = p0, acv_multiplier = acv_multiplier, numCovariates = 0L)
}

#' @title Specification of control parameters for \code{snowfall}
#'
#' @description
#' Used to specify options for snowfall
#'
#' @param parallel logical
#' @param cpus integer
#' @return A list with named elements "parallel", "cpus"
#' 
#' @export
sf.control <- function(parallel = TRUE, cpus = 2L) {
	maxc <- detectCores(logical = FALSE)
	if (cpus > maxc) {
		warning(gettextf("Number of cpus specified is greater than available cores (%i). Setting cpus at %i.", maxc, maxc))
		cpus <- maxc
	}
	list(parallel = parallel, cpus = cpus)
}


tm.param.transform <- function(rawpar, spec) {
	# Internal function to transform an atomic vector 'rawpar' of unconstrained parameters into an appropriate named list of transformed parameters
	covParOffset <- length(rawpar) - spec$numCovariates # Offset for index of covariate effects
	params <- rawpar
	if (spec$p0) {
		params[1] <- invlogit(params[1])
		params[-1] <- exp(params[-1])
	} else {
		params <- exp(params)
	}
	params <- c(as.list(params[seq(covParOffset)]), if (spec$numCovariates) list(params[-seq(covParOffset)]))
	nms <- switch(spec$family, "exponential" = "lambda", "exponential-gamma" = c("r", "alpha"))
	nms <- if (!is.null(spec$acv_multiplier)) c("gamma", nms) else nms
	nms <- if (spec$p0) c("p0", nms) else nms
	nms <- if (spec$numCovariates) c(nms, "beta") else nms
	names(params) <- nms
	params
}


tm.mle.fit <- function(dat, mfcall, startvals, tmSpec, method, optimControl) {
	# Internal function for setting up model frame and calling optim
	#cat(unique(dat$UPC), fill = TRUE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
	mfcall <- mfcall[c(1L, m)]
	mfcall[[1L]] <- as.name("model.frame")
	mfcall$data <- quote(dat)
	mf <- as.data.table(eval(mfcall))
	mt <- attr(mf, "terms")
	y <- stats::model.response(mf, "numeric")
	x <- stats::model.matrix(mt, mf)
	if (!is.null(tmSpec$acv_multiplier)) 
		acv <- dat[, tmSpec$acv_multiplier, with = FALSE]
	else
		acv <- 1
	res <- try(optim(startvals, tm.likelihood, y = y, x = x, acv = acv, tmSpec = tmSpec, method = method, control = optimControl), silent = TRUE)
	#structure(res, y = y, x = x, terms = mt)
	if (inherits(res, "try-error"))
		return(structure(list(), optres = res))
	else
		return(structure(tm.param.transform(res$par, tmSpec), optres = res))
}

tm.setstart <- function(start, spec) {
	# Internal function to check starting values / set default starting values
	numReqPars <- spec$p0 + (!is.null(spec$acv_multiplier)) + spec$numCovariates + switch(spec$fam, "exponential"=1L, "exponential-gamma"=2L)
	if (length(start) == numReqPars)
		return(start)
	else if (length(start) == 0L)
		return(rep(-1, numReqPars))
	else 
		stop(gettextf("Length of startvals (%i) does not match required length (%i)", length(start), numReqPars))
}

tm.formula.autocorrect <- function(formula) {
	# Internal function to modify a formula to drop intercept if covariates are supplied
	tm <- terms(formula)
	numCovariates <- length(attr(tm, "term.labels"))
	if (numCovariates && attr(tm, "intercept"))
		formula[[3]] <- as.call(c(as.name("-"), list(formula[[3]]), 1))
	formula
}

tm.likelihood <- function(params, y, x, acv, tmSpec) {
	# Internal function to compute negative loglikelihood
	# Name and transform parameter vector as needed
#	parNames <- switch(tmSpec$family, "exponential" = c("lambda"), "exponential-gamma" = c("r", "alpha"))
#	parNames <- if (!is.null(tmSpec$acv_multiplier)) c("gamma", parNames) else parNames
#	parNames <- if (tmSpec$p0) c("p0", parNames) else parNames
#	parNames <- if (tmSpec$numCovariates) c(parNames, "beta") else parNames
#	covParOffset <- length(params) - tmSpec$numCovariates
#	if (tmSpec$p0) {
#		params[1] <- invlogit(params[1])
#		params[-1] <- exp(params[-1])
#	} else {
#		params <- exp(params)
#	}
#	params <- c(as.list(params[seq(covParOffset)]), if (tmSpec$numCovariates) list(params[-seq(covParOffset)]))
#	names(params) <- parNames
	params <- tm.param.transform(params, tmSpec)
	incy <- diff(c(0, y))
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(c(0, predictor))
#	tpredictor <- cumFun(params, tmSpec, predictor, acv)
	tpredictor <- cumFun(params, predictor, acv)
	cdf1 <- tail(tpredictor, -1)
	cdf2 <- head(tpredictor, -1)
	cdf3 <- tail(tpredictor, 1)
	ll <- sum(incy * log(cdf1 - cdf2)) + (1 - sum(incy)) * log(1 - cdf3)
	ll
}

#cumFun <- function(params, tmSpec, time, acv) {	
cumFun <- function(params, time, acv) {	
	# Internal function to calculate cumulative trial curve as a function of time
	y <- if(is.null(params$lambda)) with(params, 1 - (alpha / (alpha + time))^r) else with(params, 1 - exp(-lambda*time))
	if (!is.null(params$p0)) 
		y <- params$p0 * y
	if (!is.null(params$gamma)) 
		y <- (acv ^ params$gamma) * y
	y
}

#' @title Maximum Likelihood Estimation for Trial Model
#'
#' @description
#' runs maximum likelihood estimation for each unit in the data and returns a list of lists 
#'
#' @param model is a list with two or three components. The first component is the dependent variable, the second component is the independent variable 
#' @param useSF logical flag which turns on usage of the \code{snowfall} package if \code{TRUE}. Default is \code{FALSE}
#' @param parallel argument to be passed onto \code{snowfall} if useSF is TRUE. The default is TRUE.
#' @param cpus argument to be passed onto \code{snowfall} if useSF is TRUE. The default is to use 4 processors.
#' @param ... are other arguments to be passed onto maximum likelihood estimation for each unit.
#' @return A list of lists with each sublist containing either the results from the estimation or an object of class \code{try-error} if maximum likelihood fails
#' 
#' @export
trialmodel.mle <- function(formula, data, group, startvals = numeric(), tmSpec = list(), sf = FALSE, sfControl = list(), method = "BFGS", optimControl = list(fnscale = -1, maxit = 20000)) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  cl$formula <- mf$formula <- tm.formula.autocorrect(formula)

  if (!inherits(formula, "formula") || length(formula) < 3L) # Check for two sided formula
	  stop("Require a formula argument")
  if (!inherits(data, "data.frame"))   # Check if data arg is a data frame
	  stop("Require a data.frame or data.table")
  if (!inherits(data, "data.table"))   # Coerce to data.table if needed
	  data <- as.data.table(data)
  if (!missing(group) && !is.character(group))
	  stop("The group argument must be a string giving the name of a column in the data argument")

  # Specification of the functional form for trial model
  tms <- tm.spec()
  if (length(tmSpec)) {
	  if (any(names(tmSpec) %notin% names(tms)))
		  warning(gettextf("Ignoring invalid elements %s of tmSpec.", paste(setdiff(names(tmSpec), names(tms)), collapse = ", ")))
	  if (!is.null(tmSpec$family) && tmSpec$family %notin% c("exponential", "exponential-gamma"))
		  stop(gettextf("family must be either %s or %s.", dQuote("exponential"), dQuote("exponential-gamma")))
	  if (!is.null(tmSpec$acv_multiplier) && tmSpec$acv_multiplier %notin% names(data))
		  stop(gettextf("%s not found in %s", deparse(tmSpec$acv_multiplier), deparse(substitute(data))))
	  cnms <- intersect(names(tmSpec), names(tms))
	  tms[cnms] <- tmSpec[cnms]
  }
  tmSpec <- tms
  tmSpec$numCovariates <- length(attr(terms(mf$formula), "term.labels"))
  # Control parameters for snowfall
  sfc <- sf.control()                    
  if (length(sfControl)) {
	  if (any(names(sfControl) %notin% names(sfc)))
		  warning(gettextf("Ignoring invalid elements %s of sfControl.", paste(names(sfControl), names(sfc), collapse = ", ")))
	  cnms <- intersect(names(sfControl), names(sfc))
	  sfc[cnms] <- sfControl[cnms]
  }
  sfControl <- sfc
  # Control parameters for optim
  if (length(optimControl)) {          # parameters for optim
	  controlopts <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
	  if (any(names(optimControl) %notin% controlopts))
		  warning(gettextf("Ignoring invalid elements %s of optimControl.", paste(setdiff(names(optimControl), controlopts), collapse = ", ")))
	  optimContol <- optimControl[intersect(names(optimControl), controlopts)]	  
  }
  # Create starting values 
  startvals <- tm.setstart(startvals, tmSpec)

  if (!missing(group) && is.character(group)) {
	  lst <- split(data, data[, group, with = FALSE])
	  ### **** PARALLEL VERSION DOES NOT WORK **** 
	  if (sf) {
		  sfInit(parallel = sfControl$parallel, cpus = sfControl$cpus, type = "SOCK")
		  #sfExport("tm.mle.fit", "tm.likelihood", "cumFun", "%notin%", "invlogit", namespace = "newprodforecast")
		  res <- sfLapply(lst, tm.mle.fit, mf, tmSpec, startvals, method, optimControl)
		  sfStop()
	  } else {
		  res <- lapply(lst, tm.mle.fit, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
	  }
  } else {
	  res <- tm.mle.fit(data, mf, startvals, tmSpec, method, optimControl)
  }
  structure(res, call = cl, startvals = startvals, spec = tmSpec, method = method, optimControl = optimControl)
}


tm.predict.onegroup <- function(dat, params, mfcall, tmSpec) {
	m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
	mfcall <- mfcall[c(1L, m)]
	mfcall[[1L]] <- as.name("model.frame")
	mfcall$data <- quote(dat)
	mf <- as.data.table(eval(mfcall))
	mt <- attr(mf, "terms")
	x <- stats::model.matrix(mt, mf)
	if (!is.null(tmSpec$acv_multiplier)) 
		acv <- dat[, tmSpec$acv_multiplier, with = FALSE]
	else
		acv <- 1
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	cumFun(params, predictor, acv)
}


predict.trialmodel <- function(obj, newdata) {
	if (missing(newdata)) {
		group <- attr(obj, "call")$group
		data <- eval(attr(obj, "call")$data, parent.frame())
		mfcall <- attr(obj, "call")
		tmSpec <- attr(obj, "spec")
		if (is.null(group)) {
			if (length(obj))
				tm.predict.onegroup(data, obj, mfcall, tmSpec)
			else 
				simpleError("Model did not converge.")
		} else {
			lst <- split(data, data[, group, with = FALSE])
			mapply(tm.predict.onegroup, lst, obj, MoreArgs = list(mfcall = enquote(mfcall), tmSpec = tmSpec), SIMPLIFY = FALSE)
		}
	}
}



if (FALSE) {
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
}
