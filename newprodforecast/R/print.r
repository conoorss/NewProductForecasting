print.trialmodel <- function(x) {
	linesep <- function() cat(paste(rep("-", 60), collapse=""), fill = TRUE)
	printsafe <- function(x) {
		if (NROW(x) <= 20) {
			print(x)
		} else {
			print(head(x))
			cat("...", fill = TRUE)
			print(tail(x))
		}
	}
	options(digits = 3)
	on.exit(options(digits = 7))
	coefs <- coefficients(x)
	obj <- if (x$estimation == "MLE") likelihood(x) else sumsqerr(x)
	mape <- insamplefit(x)
	if (x$estimation == "MLE") 
		df_fit <- data.frame(Likelihood = obj, MAPE = mape)
	else if (x$estimation == "OLS")
		df_fit <- data.frame(`Sum of Squared Error` = obj, MAPE = mape)

	cat("Call:", fill = TRUE)
	print(x$call)
	cat("\nCoefficients:", fill = TRUE)
	linesep()
	printsafe(coefs)
	cat("\nModel Fit:", fill = TRUE)
	linesep()
	printsafe(df_fit)
}

