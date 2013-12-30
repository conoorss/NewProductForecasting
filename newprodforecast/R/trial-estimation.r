#' @title Specification for Trial Model
#'
#' @description
#' Used to specify options for the trial model:
#' functional family of the cumulative trial curve
#' presence of asymptote parameter
#' presence of multiplicative term for ACV (distribution)
#'
#' @param lst A list with the specification of the model 
#' @return A list with components "family", "p0", "acvMultiplier", "numCovariates
#' 
#' @export

setTrialModelSpec <- function(lst) {
	spec <- default <- list(family = "exponential-gamma", p0 = FALSE, acvMultiplier = NULL, numCovariates = 0L)
	if (length(lst)) {
		if (is.null(names(lst)) || any(identical(names(lst), "")))
			stop("Need a named list of values for model specification")
		if (any(names(lst) %notin% names(default)))
			warning(gettextf("Ignoring unknown options {%s} of model specification",
							 paste(setdiff(names(lst), names(lst)), collapse = ", ")))
		if (!is.null(lst$family) && lst$family %notin% c("exponential", "exponential-gamma"))
			stop(gettextf("family must be either %s or %s.", dQuote("exponential"), dQuote("exponential-gamma")))
		common <- intersect(names(lst), names(default))
		spec[common] <- lst[common]
	}
	return(spec)
}


#' @title Specification of control parameters for \code{snowfall}
#'
#' @description
#' Used to specify options for snowfall
#'
#' @param lst A list of control parameters 
#' @return A list with named elements "parallel", "cpus"
#' 
#' @export

setSnowfallControls <- function(lst) {
	control <- default <- list(parallel = TRUE, cpus = 2L) 
	maxCores <- detectCores(logical = FALSE)
	if (length(lst)) {
		if (is.null(names(lst)) || any(identical(names(lst), "")))
			stop("Need a named list of values for snowfall control")
		if (any(names(lst) %notin% names(default)))
			warning(gettextf("Ignoring unknown options {%s} of model specification",
							 paste(setdiff(names(lst), names(default)), collapse = ", ")))
		if (!is.null(lst$parallel) && !is.logical(lst$parallel))
			stop("%s must be a logical argument", dQuote("parallel"))
		if (!is.null(lst$cpus) && cpus > maxCores) {
			warning(gettextf("Number of cpus specified is greater than available cores (%i). Setting cpus at %i.", maxCores, maxCores))
			lst$cpus <- maxCores
		}
		common <- intersect(names(lst), names(default))
		control[common] <- lst[common]
	}
	return(control)
}

setOptimControls <- function(lst) {
	controlopts <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")

	if (length(lst)) {          # parameters for optim
		if (any(names(lst) %notin% controlopts))
			warning(gettextf("Ignoring invalid elements %s of optimControl.", 
							 paste(setdiff(names(lst), controlopts), collapse = ", ")))
	lst <- lst[intersect(names(lst), controlopts)]	  
	}
	lst
}


transformTrialParams <- function(rawpar, spec) {
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
	nms <- if (!is.null(spec$acvMultiplier)) c("gamma", nms) else nms
	nms <- if (spec$p0) c("p0", nms) else nms
	nms <- if (spec$numCovariates) c(nms, "beta") else nms
	names(params) <- nms
	params
}


mleTrialModel <- function(dat, mfcall, startvals, tmSpec, method, optimControl) {
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
	if (!is.null(tmSpec$acvMultiplier)) 
		acv <- dat[, tmSpec$acvMultiplier, with = FALSE]
	else
		acv <- 1
	res <- try(optim(startvals, likelihoodTrialModel, y = y, x = x, acv = acv, tmSpec = tmSpec, method = method, control = optimControl), silent = TRUE)
	#structure(res, y = y, x = x, terms = mt)
	if (inherits(res, "try-error"))
		list(coefficients = list(), optres = res)
	else
		list(coefficients = transformTrialParams(res$par, tmSpec), optres = res)
}


olsTrialModel <- function(dat, mfcall, startvals, tmSpec, method, optimControl) {
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
	if (!is.null(tmSpec$acvMultiplier)) 
		acv <- dat[, tmSpec$acvMultiplier, with = FALSE]
	else
		acv <- 1
	res <- try(optim(startvals, sqErrorTrialModel, y = y, x = x, acv = acv, tmSpec = tmSpec, method = method, control = optimControl), silent = TRUE)
	#structure(res, y = y, x = x, terms = mt)
	if (inherits(res, "try-error"))
		list(coefficients = list(), optres = res)
	else
		list(coefficients = transformTrialParams(res$par, tmSpec), optres = res)
}


setTrialModelStart <- function(start, spec) {
	# Internal function to check starting values / set default starting values
	numReqPars <- spec$p0 + (!is.null(spec$acvMultiplier)) + spec$numCovariates + switch(spec$fam, "exponential"=1L, "exponential-gamma"=2L)
	if (length(start) == numReqPars)
		return(start)
	else if (length(start) == 0L)
		return(rep(-1, numReqPars))
	else 
		stop(gettextf("Length of startvals (%i) does not match required length (%i)", length(start), numReqPars))
}

autocorrectModelFormula <- function(formula) {
	# Internal function to modify a formula to drop intercept if covariates are supplied
	tm <- terms(formula)
	numCovariates <- length(attr(tm, "term.labels"))
	if (numCovariates && attr(tm, "intercept"))
		formula[[3]] <- as.call(c(as.name("-"), list(formula[[3]]), 1))
	formula
}

likelihoodTrialModel <- function(params, y, x, acv, tmSpec) {
	# Internal function to compute negative loglikelihood
	params <- transformTrialParams(params, tmSpec)
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

sqErrorTrialModel <- function(params, y, x, acv, tmSpec) {
	# Internal function to compute squared error
	params <- transformTrialParams(params, tmSpec)
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	tpredictor <- cumFun(params, predictor, acv)
	sum((y - tpredictor)^2)
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
trialmodel <- function(formula, 
					   data, 
					   group, 
					   startvals = numeric(), 
					   tmSpec = list(), 
					   sf = FALSE, sfControl = list(), 
					   estimation = c("MLE", "OLS"),  
					   method = "Nelder-Mead", 
					   optimControl = list(maxit = 20000)) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  cl$formula <- mf$formula <- autocorrectModelFormula(formula)

  if (!inherits(formula, "formula") || length(formula) < 3L) # Check for two sided formula
	  stop("Require a formula argument")
  if (!inherits(data, "data.frame"))   # Check if data arg is a data frame
	  stop("Require a data.frame or data.table")
  if (!inherits(data, "data.table"))   # Coerce to data.table if needed
	  data <- as.data.table(data)
  if (!missing(group) && !is.character(group))
	  stop("The group argument must be a string giving the name of a column in the data argument")

  estimation <- match.arg(estimation)
					   
  # Specification of the functional form for trial model
  tmSpec <- setTrialModelSpec(tmSpec)
  tmSpec$numCovariates <- length(attr(terms(mf$formula), "term.labels"))
  if (!is.null(tmSpec$acvMultiplier) && tmSpec$acvMultiplier %notin% names(data))
	  stop("The acv multiplier variable in the model specification was not found in the data")
  
  sfControl <- setSnowfallControls(sfControl)	# Control parameters for snowfall
  optimControl <- setOptimControls(optimControl)	# Control parameters for optim
  if (estimation == "MLE") optimControl$fnscale <- -1 # Need maximization if maximum likelihood estimation
  startvals <- setTrialModelStart(startvals, tmSpec)	# Create starting values

  if (!missing(group) && is.character(group)) {
	  lst <- split(data, data[, group, with = FALSE])
	  ### **** PARALLEL VERSION DOES NOT WORK **** 
	  if (sf) {
		  sfInit(parallel = sfControl$parallel, cpus = sfControl$cpus, type = "SOCK")
		  #sfExport("mleTrialModel", "likelihoodTrialModel", "cumFun", "%notin%", "invlogit", namespace = "newprodforecast")		  
		  if (estimation == "MLE")
			  res <- sfLapply(lst, mleTrialModel, mf, tmSpec, startvals, method, optimControl)
		  else
			  message("OLS not yet implemented")
		  sfStop()
	  } else {
		  if (estimation == "MLE")
			  res <- lapply(lst, mleTrialModel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
		  else if (estimation == "OLS")
			 res <- lapply(lst, olsTrialModel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
	  }
  } else {
	  if (estimation == "MLE")
		  res <- mleTrialModel(data, mf, startvals, tmSpec, method, optimControl)
	  else
		  res <- olsTrialModel(data, mf, startvals, tmSpec, method, optimControl)
  }
  out <- list(estimates = res, call = cl, startvals = startvals, spec = tmSpec, estimation = estimation, method = method, optimControl = optimControl)
  class(out) <- "trialmodel"
  out
}


predictTrialOneGroup <- function(dat, params, mfcall, tmSpec) {
	m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
	mfcall <- mfcall[c(1L, m)]
	mfcall[[1L]] <- as.name("model.frame")
	mfcall$data <- quote(dat)
	mf <- as.data.table(eval(mfcall))
	mt <- attr(mf, "terms")
	x <- stats::model.matrix(mt, mf)
	if (!is.null(tmSpec$acvMultiplier)) 
		acv <- dat[, tmSpec$acvMultiplier, with = FALSE]
	else
		acv <- 1
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	cumFun(params, predictor, acv)
}

predict.trialmodel <- function(obj, newdata) {
	# Handling the call object
	mfcall <- obj$call
	callnms <- names(mfcall)
	if ("group" %in% callnms)       # match by name
		group <- mfcall$group
	else if (sum(callnms == "") > 3L)    # match by position
		group <- mfcall[[4L]]
	else
		group <- NULL
	data <- if(missing(newdata)) eval(mfcall[[3L]], parent.frame()) else newdata
	if (is.null(group)) {
		coefs <- obj$estimates$coefficients
		if (length(obj$estimates$coefficients))
			predictTrialOneGroup(data, coefs, mfcall, obj$spec)
		else 
			simpleError("Model did not converge.")
	} else {
		coefs <- lapply(seq_along(obj$estimates), function(i) obj$estimates[[i]]$coefficients)
		lst <- split(data, data[, group, with = FALSE])
		mapply(predictTrialOneGroup, lst, coefs, MoreArgs = list(mfcall = enquote(mfcall), tmSpec = obj$spec), SIMPLIFY = FALSE)
	}
}

