runTest <- function(call) {
	cat("\n***********************************************************\n")
	print(call)
	time <- system.time( res <- try(eval(call), silent = TRUE) )
	if (any(sapply(res$estimates, function(x) inherits(x, "try-error"))))
		message("STATUS: FAILED")
	else
		message("STATUS: PASSED")
	message(gettextf("Time: %f secs", time[["elapsed"]]))
	cat("\n***********************************************************\n")
	structure(res, call = call)
}


data("trial-model-input-120313-01")
load_all()
uupcs <- model.upc.exp[, unique(UPC)]
testdata1 <- model.upc.exp[UPC %in% uupcs[1]]
testdata2 <- model.upc.exp[UPC %in% uupcs[1:10]]

test1 <- runTest(quote(trialmodel(ObsTrialPct ~ 1, testdata1, estimation = "MLE")))
test2 <- runTest(quote(trialmodel(ObsTrialPct ~ 1, testdata1, estimation = "OLS")))
test3 <- runTest(quote(trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata1, estimation = "MLE")))
test4 <- runTest(quote(trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata1, estimation = "OLS")))
test5 <- runTest(quote(trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata2, "UPC", estimation = "MLE")))
test6 <- runTest(quote(trialmodel(ObsTrialPct ~ ACV.MultiOutlet, testdata2, "UPC", estimation = "OLS")))
