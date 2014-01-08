
data("trial-model-input-120313-01")
load_all()
uupcs <- model.upc.exp[, unique(UPC)]
testdata1 <- model.upc.exp[UPC %in% uupcs[1]]
testdata2 <- model.upc.exp[UPC %in% uupcs[1:10]]
testdata2 <- na.omit(testdata2)

# Run MLE to get starting points
mle1 <- trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata2, "UPC", estimation = "MLE")
paramsList <- lapply(mle1$estimates, "[[", i = "coefficients")
revt <- generate_transforms(mle1$spec, reverse = TRUE)
paramsList <- lapply(paramsList, function(x) apply_transforms(x, revt))
#                      mapply(function(y, fun) fun(y), x, revt, SIMPLIFY = FALSE))

tmpmat <- do.call("rbind", lapply(paramsList, unlist))
Delta <- apply(tmpmat, 2, median)
Sigma <- diag(3)



# Get initial estimate
sigsq <- log(var(do.call("c", lapply(mle1$estimates, "[[", i = "residuals"))))


Data <- get_xy(testdata2, "UPC", "ObsTrialPct", "ACV.MultiOutlet")
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
			 rscale = 0.01, alphascale = 0.01, betascale = 0.01, sisqscale = 0.01, 
			 burn = 1000, samples = 0, thin = 1, printThin = 10)

res <- hb_trialmodel(Data, Priors, Mcmc)
