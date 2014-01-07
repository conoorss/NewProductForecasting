#' @title Specification of control parameters for \code{snowfall}
#'
#' @description
#' Used to specify options for snowfall
#'
#' @param lst A list of control parameters 
#' @return A list with named elements "parallel", "cpus"
#' 
#' @export

set_controls_snowfall <- function(lst) {
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

#' @title Specification of control parameters for \code{optim}
#'
#' @description
#' Used to specify options for optim
#'
#' @param lst a named list 
#' @return a valid list of inputs for optim. See \code{optim} for details
#' 
#' @export
set_controls_optim <- function(lst) {
	controlopts <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")

	if (length(lst)) {          # parameters for optim
		if (any(names(lst) %notin% controlopts))
			warning(gettextf("Ignoring invalid elements %s of optimControl.", 
							 paste(setdiff(names(lst), controlopts), collapse = ", ")))
	lst <- lst[intersect(names(lst), controlopts)]	  
	}
	lst
}
