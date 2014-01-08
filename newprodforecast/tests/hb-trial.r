
data("trial-model-input-120313-01")
load_all()
uupcs <- model.upc.exp[, unique(UPC)]

limit <- length(uupcs)
testdata1 <- model.upc.exp[UPC %in% uupcs[1]]
testdata2 <- model.upc.exp[UPC %in% uupcs[1:limit]]
testdata2 <- na.omit(testdata2)

# Run MLE to get starting points
mle1 <- trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata2, "UPC", estimation = "MLE")
#get_coef <- function(x) if (length(x$coefficients)) x$coefficients else as.relistable(list(r = exp(-1), alpha = exp(-1), beta = exp(-1)))
paramsList <- lapply(mle1$estimates, "[[", i="coefficients")

fill_miss_params <- function(paramsList) {
	whichMiss <- which(sapply(paramsList, length) == 0)
	paramsListNonMiss <- paramsList[-whichMiss]
	r_median <- median(sapply(paramsListNonMiss, "[[", i = "r"))
	alpha_median <- median(sapply(paramsListNonMiss, "[[", i = "alpha"))
	beta_median <- apply(do.call("rbind", lapply(paramsListNonMiss, "[[", i = "beta")), 2, median)
	fill_value <- as.relistable(list(r = r_median, alpha = alpha_median, beta = beta_median))
	for (i in whichMiss) paramsList[[i]] <- fill_value
	paramsList
}

paramsList <- fill_miss_params(paramsList)

revt <- generate_transforms(mle1$spec, reverse = TRUE)
paramsList <- lapply(paramsList, function(x) apply_transforms(x, revt))
#                      mapply(function(y, fun) fun(y), x, revt, SIMPLIFY = FALSE))

tmpmat <- do.call("rbind", lapply(paramsList, unlist))
Delta <- apply(tmpmat, 2, median)
Sigma <- diag(3)



# Get initial estimate
sigsq <- log(var(do.call("c", lapply(mle1$estimates, "[[", i = "residuals"))))


Data <- get_xylist(testdata2, "UPC", "ObsTrialPct", "ACV.MultiOutlet")
Data$z <- matrix(1, length(paramsList), 1)

numCovariates <- NCOL(Data$x[[1]])
numUpperCovariates <- NCOL(Data$z)
numLowerPars <- 2 + numCovariates


Priors <- list(Deltabar = matrix(0, numUpperCovariates, numLowerPars), 
			   A = 0.01 * matrix(1, numUpperCovariates, numUpperCovariates), 
			   nu = numLowerPars + 3, V = (numLowerPars + 3) * diag(numLowerPars), 
			   sigsq_mean = 0, sigsq_sd = 10)

Mcmc <- list(paramsList = paramsList, 
			 sigsq = sigsq, 
			 Delta = Delta, 
			 Sigma = Sigma, 
			 rscale = 1.0, alphascale = 1.0, betascale = 1.0, sigsqscale = 1.5, 
			 burn = 100, samples = 0, thin = 1, printThin = 10)

res <- hb_trialmodel(Data, Priors, Mcmc)
