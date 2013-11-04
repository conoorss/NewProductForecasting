source("../try-2.r")
load("campbells-model-data.rdata")

objfn <- function(params, y, time) {
  params[1] <- invlogit(params[1])
  params[-1] <- exp(params[-1])
  names(params) <- c("p0", "r", "alpha")
  ll <- loglikelihood(params, y, time, 1)
  #print(ll)
  ll
}

mle_oneid <- function(dat, startvals = c(-1, -1, -1), ...) {
  y <- dat[, ObsTrialPct]
  time <- seq_along(y)
  try(optim(startvals, objfn, y = y, time = time, control = list(fnscale = -1, maxit = 20000), method = "BFGS"))  
}

mle_allids <- function(lst, ...) {
  res <- vector("list", length(lst))
  start_time <- proc.time()[3]
  for (i in seq_along(lst)) {
    if (i %% 10 == 0L) {
      curr_time <- proc.time()[3]
      cat(i, "ids done in ", curr_time - start_time, " secs\n")
    }
    res[[i]] <- mle_oneid(lst[[i]], ...)
  }
  res
}

reportpars <- function(mle) {
  if (!inherits(mle, "try-error") && ("par" %in% names(mle))) {
    par <- mle$par
    return(list(p0 = invlogit(par[1]), r = exp(par[2]), alpha = exp(par[3])))
  } else {
    return(list())
  }   
}

curvefit <- function(mle, dat){
  res <- reportpars(mle)
  if (length(res)) {
    fit <- cdf(res, seq(nrow(dat)))
    mape <- mean(abs(fit/(dat[, ObsTrialPct] + 1e-5) - 1))
    return(cbind(dat, PredTrialPct = fit, MAPE = mape, p0 = res$p0, r = res$r, alpha = res$alpha))
  } else {
    return(data.table())
  }
}



# Run MLE
mle_oneid(campdat2lst[[1]])
mle.test <- mle_allids(campdat2lst)

# Check convergence status
convcheck <- sapply(mle.test, function(x) if (inherits(x, "try-error")) NA else x$convergence)
table(convcheck, useNA = "ifany")

# Get model fit
modelfit <- mapply(curvefit, mle.test, campdat2lst, SIMPLIFY = FALSE)
modelfit <- rbindlist(modelfit)
modelfit[,`:=`(id = NULL, acv.max = NULL, acv.range = NULL)]

estpars <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("Product", "NPP_BRAND_NUM", "NPP_BRAND", "p0", "r", "alpha")]
mapes <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("MAPE")]

##### Plot model fit
for (i in unique(modelfit$UPC)) {
  mat <- modelfit[UPC == i, as.matrix(list(ObsTrialPct, PredTrialPct))]
  matplot(mat, type = "l", main = paste("UPC:", i), ylab = "Cum. Trial Rate")
  readline("cont?")
}

write.csv(modelfit, file = "UPC-level-mle-2.csv", row.names = FALSE)

save(modelfit, estpars, mapes, file = "model-res-2.rdata")