data("trial-model-input-120313-01")
load_all()
uupcs <- model.upc.exp[, unique(UPC)]
testdata1 <- model.upc.exp[UPC %in% uupcs[1]]
testdata2 <- model.upc.exp[UPC %in% uupcs[1:10]]

#if (FALSE){
cat("\n\n")
message("Test 1: Running test with one UPC, no covariates")
st1 <- system.time( testmodel1 <- try(trialmodel.mle(ObsTrialPct ~ 1, data = testdata1), silent = TRUE) )
message(gettextf("Test 1: %s", if (inherits(testmodel1, "try-error")) "FAILED" else "PASSED"))
message(gettextf("Test 1: %f secs", st1[3]))
cat("\n\n")

message("Test 2: Running test with one UPC, 1 covariate")
st2 <- system.time( testmodel2 <- try(trialmodel.mle(ObsTrialPct ~ ACV.MultiOutlet, data = testdata1), silent = TRUE) )
message(gettextf("Test 2: %s", if (inherits(testmodel2, "try-error")) "FAILED" else "PASSED"))
message(gettextf("Test 2: %f secs", st2[3]))
cat("\n\n")

message("Test 3: Running test with ten UPCs, 1 covariate")
st3 <- system.time( testmodel3 <- trialmodel.mle(ObsTrialPct ~ ACV.MultiOutlet, data = testdata2, group = "UPC") )
message(gettextf("Test 3: %s", if (any(sapply(testmodel3, function(x) inherits(x, "try-error")))) "FAILED" else "PASSED"))
message(gettextf("Test 3: %f secs", st3[3]))
cat("\n\n")
#}

if (FALSE) {
message("Test 4: Running test with ten UPCs, 1 covariate, using snowfall (6 cores)")
sfc <- list(parallel = TRUE, cpus = 2)
st4 <- system.time( testmodel4 <- trialmodel.mle(ObsTrialPct ~ ACV.MultiOutlet, data = testdata2, group = "UPC", sf = TRUE, sfControl = sfc) )
message(gettextf("Test 4: %s", if (inherits(testmodel4, "try-error")) "FAILED" else "PASSED"))
message(gettextf("Test 4: %f secs", st4[3]))
cat("\n\n")
}
