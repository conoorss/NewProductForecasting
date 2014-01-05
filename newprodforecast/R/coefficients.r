coef.trialmodel <- function(object, collapse = TRUE) {
	groupsPresent <- !is.null(groupvar(object))
	if (groupsPresent) {
		res <- lapply(object$estimates, get_coefficients)
		if (collapse)
			return(do.call("rbind", res))
		else 
			return(res)
	} else {
		res <- get_coefficients(object$estimates)
		if (collapse) 
			return(unlist(res))
		else 
			return(res)
	}
}

get_coefficients <- function(estimates) estimates$coefficients
