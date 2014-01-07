#' @title Estimation for the Repeat Model
#' 
#' @description
#' Fits a repeat model using either least squares. Returns an object of class \code{repeatmodel}
#
#' @param formula is a two-sided formula object which specifies the dependent variable and the covariates if any.
#' @param data is either a \code{data.frame} or \code{data.table} with the variables needed for the model.
#' @param group is a string with the name of the group variable.
#' @param startvals is a numeric vector with the starting values for the model. 
#' @param repSpec is a list which specifies the model specification (family, p0, acvMultiplier, number of covariates).
#' @param sf is a logical flag for usage of the \code{snowfall} package for parallelization (currently not implemented)
#' @param sfControl is a list of control parameters for \code{snowfall}
#' @param estimation is a string which is either "MLE" or "OLS"
#' @param method is a string which specifies the optimization method
#' @param optimControl is a list of control parameters for \code{optim}
#'
#' @return A list of lists with each sublist containing either the results from the estimation or an object of class \code{try-error} if maximum likelihood fails
#' 
#' @export
repeatmodel <- function(formula, 
						data, 
						group,
						startvals = numeric(),
						repSpec = list(),
						sf = FALSE,
						sfControl = list(),
						estimation = c("OLS"),
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
	# Specification of the functional form for repeat model
	repSpec <- set_spec_repeatmodel(repSpec)
	repSpec$numCovariates <- length(attr(terms(mf$formula), "term.labels"))
	if (!is.null(repSpec$acvMultiplier) && repSpec$acvMultiplier %notin% names(data))
		stop("The acv multiplier variable in the model specification was not found in the data")

	sfControl <- set_controls_snowfall(sfControl)	# Control parameters for snowfall
	optimControl <- set_controls_optim(optimControl)	# Control parameters for optim
	if (estimation == "MLE") optimControl$fnscale <- -1 # Need maximization if maximum likelihood estimation
	startvals <- set_start_repeatmodel(startvals, repSpec)	# Create starting values

	if (!missing(group) && is.character(group)) {
		lst <- split(data, data[, group, with = FALSE])
		### **** PARALLEL VERSION DOES NOT WORK **** 
		if (sf) {
			warning("The snowfall implementation is not working. Sorry! Switching to serial estimation.")
			if (estimation == "MLE")
				res <- lapply(lst, mle_repeatmodel, mfcall = mf, repSpec = repSpec, startvals = startvals, method = method, optimControl = optimControl)
			else if (estimation == "OLS")
				res <- lapply(lst, ols_repeatmodel, mfcall = mf, repSpec = repSpec, startvals = startvals, method = method, optimControl = optimControl)

			#sfInit(parallel = sfControl$parallel, cpus = sfControl$cpus, type = "SOCK")
			#sfExport("mle_repeatmodel", "loglik_repeatmodel", "cumulative_curve", "%notin%", "invlogit", namespace = "newprodforecast")		  
			#if (estimation == "MLE")
			#  res <- sfLapply(lst, mle_repeatmodel, mf, repSpec, startvals, method, optimControl)
			#else
			#  message("OLS not yet implemented")
			#sfStop()
		} else {
			if (estimation == "MLE")
				res <- lapply(lst, mle_repeatmodel, mfcall = mf, repSpec = repSpec, startvals = startvals, method = method, optimControl = optimControl)
			else if (estimation == "OLS")
				res <- lapply(lst, ols_repeatmodel, mfcall = mf, repSpec = repSpec, startvals = startvals, method = method, optimControl = optimControl)
		}
	} else {
		if (estimation == "MLE")
			res <- mle_repeatmodel(data, mf, startvals, repSpec, method, optimControl)
		else
			res <- ols_repeatmodel(data, mf, startvals, repSpec, method, optimControl)
	}
	out <- list(estimates = res, call = cl, startvals = startvals, spec = repSpec, estimation = estimation, method = method, optimControl = optimControl)
	class(out) <- "repeatmodel"
	out
}


