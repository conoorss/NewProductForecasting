get_objval <- function(estimates) estimates$optres$value
get_mape <- function(estimates) estimates$mape

likelihood <- function(object) UseMethod("likelihood")
likelihood.trialmodel <- function(object) {
	if (object$estimation != "MLE") 
		stop(gettextf("Model was estimated using %s. No likelihood available", object$estimation))
	groupsPresent <- !is.null(groupvar(object))
	#groupsPresent <- all(sapply(object$estimates, is.list))
	if (groupsPresent)
		return(do.call("c", lapply(object$estimates, get_objval)))
	else
		return(get_objval(object$estimates))
}

sumsqerr <- function(object) UseMethod("sumsqerr")
sumsqerr.trialmodel <- function(object) {
	if (object$estimation != "OLS") 
		stop(gettextf("Model was estimated using %s. No sum of squared error available", object$estimation))
	groupsPresent <- !is.null(groupvar(object))
	#groupsPresent <- all(sapply(object$estimates, is.list))
	if (groupsPresent)
		return(do.call("c", lapply(object$estimates, get_objval)))
	else
		return(get_objval(object$estimates))
}


insamplefit <- function(object) UseMethod("insamplefit") 
insamplefit.trialmodel <- function(object) {
	#	groupsPresent <- all(sapply(object$estimates, is.list))
	groupsPresent <- !is.null(groupvar(object))
	if (groupsPresent)
		return(do.call("c", lapply(object$estimates, get_mape)))
	else
		return(get_mape(object$estimates))
}

groupvar <- function(object) object$call$group

