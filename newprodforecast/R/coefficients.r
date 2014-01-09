#' @title Method for retrieving coefficients from trial model object
#' 
#' @description
#' If the trial model was estimated by groups the method returns a matrix of groups by coefficients if the collapse argument is TRUE and a list of lists otherwise. 
#' 
#' @param object An object of class \code{trialmodel}
#' @collapse logical. If TRUE then the result is simplified to an atomic (vector/matrix) object. Otherwise a list or list of lists is returned.
#'
#' @export
coef.trialmodel <- function(object, collapse = TRUE) {
	groupsPresent <- !is.null(groupvar(object))
	if (groupsPresent) {
		res <- lapply(object$estimates, "[[", i = "coefficients")
		if (collapse) {
			res <- lapply(res, unlist)
			return(do.call("rbind", res))
		} else {
			return(res)
		}
	} else {
		res <- object$estimates$coefficients
		if (collapse) 
			return(unlist(res))
		else 
			return(res)
	}
}

#get_coefficients <- function(estimates) estimates$coefficients
