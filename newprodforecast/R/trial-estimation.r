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
  cl$formula <- mf$formula <- autocorrect_formula(formula)

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
  tmSpec <- set_spec_trialmodel(tmSpec)
  tmSpec$numCovariates <- length(attr(terms(mf$formula), "term.labels"))
  if (!is.null(tmSpec$acvMultiplier) && tmSpec$acvMultiplier %notin% names(data))
	  stop("The acv multiplier variable in the model specification was not found in the data")
  
  sfControl <- set_controls_snowfall(sfControl)	# Control parameters for snowfall
  optimControl <- set_controls_optim(optimControl)	# Control parameters for optim
  if (estimation == "MLE") optimControl$fnscale <- -1 # Need maximization if maximum likelihood estimation
  startvals <- set_start_trialmodel(startvals, tmSpec)	# Create starting values

  if (!missing(group) && is.character(group)) {
	  lst <- split(data, data[, group, with = FALSE])
	  ### **** PARALLEL VERSION DOES NOT WORK **** 
	  if (sf) {
      warning("The snowfall implementation is not working. Sorry! Switching to serial estimation.")
      if (estimation == "MLE")
        res <- lapply(lst, mle_trialmodel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
      else if (estimation == "OLS")
        res <- lapply(lst, ols_trialmodel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
      
		  #sfInit(parallel = sfControl$parallel, cpus = sfControl$cpus, type = "SOCK")
		  #sfExport("mle_trialmodel", "loglik_trialmodel", "cumulative_curve", "%notin%", "invlogit", namespace = "newprodforecast")		  
		  #if (estimation == "MLE")
			#  res <- sfLapply(lst, mle_trialmodel, mf, tmSpec, startvals, method, optimControl)
		  #else
			#  message("OLS not yet implemented")
		  #sfStop()
	  } else {
		  if (estimation == "MLE")
			  res <- lapply(lst, mle_trialmodel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
		  else if (estimation == "OLS")
			 res <- lapply(lst, ols_trialmodel, mfcall = mf, tmSpec = tmSpec, startvals = startvals, method = method, optimControl = optimControl)
	  }
  } else {
	  if (estimation == "MLE")
		  res <- mle_trialmodel(data, mf, startvals, tmSpec, method, optimControl)
	  else
		  res <- ols_trialmodel(data, mf, startvals, tmSpec, method, optimControl)
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
	cumulative_curve(params, predictor, acv)
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

