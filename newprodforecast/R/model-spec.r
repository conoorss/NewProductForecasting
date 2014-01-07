#' @title Functions for checking/setting the specification of the trial/repeat models
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
set_spec_trialmodel <- function(lst) {
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

#' @export
set_spec_repeatmodel <- function(lst) {
	spec <- default <- list(family = "exponential-gamma", acvMultiplier = NULL, numCovariates = 0L)
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


