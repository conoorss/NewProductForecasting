set_start_trialmodel <- function(start, spec) {
  # Internal function to check starting values / set default starting values
  numReqPars <- spec$p0 + (!is.null(spec$acvMultiplier)) + spec$numCovariates + switch(spec$fam, "exponential"=1L, "exponential-gamma"=2L)
  if (length(start) == numReqPars)
    return(start)
  else if (length(start) == 0L)
    return(rep(-1, numReqPars))
  else 
    stop(gettextf("Length of startvals (%i) does not match required length (%i)", length(start), numReqPars))
}