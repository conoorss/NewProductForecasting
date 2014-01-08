
"%notin%" <- Negate("%in%")
invlogit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x) - log (1 - x)

autocorrect_formula <- function(formula) {
	# Internal function to modify a formula to drop intercept if covariates are supplied
	tm <- terms(formula)
	numCovariates <- length(attr(tm, "term.labels"))
	if (numCovariates && attr(tm, "intercept"))
		formula[[3]] <- as.call(c(as.name("-"), list(formula[[3]]), 1))
	formula
}
	
linesep <- function(n = 60) cat(paste(rep("-", n), collapse=""), fill = TRUE)

printsafe <- function(x) {
	if (NROW(x) <= 20) {
		print(x)
	} else {
		print(head(x))
		cat("...", fill = TRUE)
		print(tail(x))
	}
}
