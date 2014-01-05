generate_vec_splitter <- function(spec) {
  # Generates a function to be used to split a vector of "raw" parameters into a named list  
  numPars <- switch(spec$family, "exponential" = 1, "exponential-gamma" = 2) + spec$p0 + (!is.null(spec$acvMultiplier)) + spec$numCovariates
  covOffset <- numPars - spec$numCovariates
  splitvec <- c(seq(covOffset), if (spec$numCovariates) rep(covOffset + 1L, spec$numCovariates) else c())  
  labs <- c(if (spec$p0) "p0" else c(), 
            if(!is.null(spec$acvMultiplier)) "gamma" else c(), 
            switch(spec$family, "exponential" = c("lambda"), "exponential-gamma" = c("r", "alpha"), stop("Invalid family in model specification")), 
            if (spec$numCovariates) "beta" else c())
  splitvec <- ordered(splitvec, labels = labs)
  
  #names(splitvec) <- nms
  function(x) if (length(x) != length(splitvec)) stop("Length mismatch in splitter") else split(x, splitvec)
}

generate_transform_funss <- function(spec) {
  tfunc <- if (spec$p0) list(p0 = invlogit) else list()
  tfunc <- if (!is.null(spec$acvMultiplier)) c(tfunc, list(gamma = exp)) else tfunc
  tfunc <- c(tfunc, switch(spec$family, "exponential" = list(lambda = exp), "exponential-gamma" = list(r = exp, alpha = exp), stop("Invalid family in model specification")))
  tfunc <- if (spec$numCovariates > 0L) c(tfunc, list(beta = exp)) else tfunc
  tfunc
}

if(0){

transformTrialParams <- function(rawpar, spec) {
  # Internal function to transform an atomic vector 'rawpar' of unconstrained parameters into an appropriate named list of transformed parameters
  covParOffset <- length(rawpar) - spec$numCovariates # Offset for index of covariate effects
  params <- rawpar
  if (spec$p0) {
    params[1] <- invlogit(params[1])
    params[-1] <- exp(params[-1])
  } else {
    params <- exp(params)
  }
  params <- c(as.list(params[seq(covParOffset)]), if (spec$numCovariates) list(params[-seq(covParOffset)]))
  nms <- switch(spec$family, "exponential" = "lambda", "exponential-gamma" = c("r", "alpha"))
  nms <- if (!is.null(spec$acvMultiplier)) c("gamma", nms) else nms
  nms <- if (spec$p0) c("p0", nms) else nms
  nms <- if (spec$numCovariates) c(nms, "beta") else nms
  names(params) <- nms
  params
}

}
