source("../try-3.r")
load("campbells-model-data.rdata")

objfn <- function(params, y, X) {
  params[1] <- invlogit(params[1])
  params[-1] <- exp(params[-1])
  params <- c(as.list(params[1:3]), list(params[-(1:3)]))
  names(params) <- c("p0", "r", "alpha", "beta")
  if (!identical(length(params$beta), ncol(X)))
    stop("Dimension mismatch between beta and X")
  ll <- loglikelihood(params, y, X, 1)
  #print(ll)
  ll
}

mle_oneid <- function(dat, xvars, startvals = c(-1, -1, -1, rep(-1, length(xvars))), ...) {
  y <- dat[, ObsTrialPct]
  X <- dat[, xvars, with = FALSE]
  X <- as.matrix(X)
  time <- seq_along(y)
  try(optim(startvals, objfn, y = y, X = X, control = list(fnscale = -1, maxit = 20000), method = "BFGS"))  
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

mle1 <- function(lst){
  res <- vector("list", length(lst))
  for (i in seq_along(lst)) {
    if (i %% 10 == 0) cat(i, " done \n")
    yi <- lst[[i]][, ObsTrialPct]
    Xi <- lst[[i]][, ACV.MultiOutlet]
    Xi <- as.matrix(Xi)
    timei <- lst[[i]][, WEEK]
    res[[i]] <- try(optim(c(-1, -1, -1, -1), objfn, y = yi, X = Xi, control = list(fnscale = -1, maxit = 20000), method = "BFGS"))
  }
  res
}

reportpars <- function(mle) {
   par <- mle$par
   p0 <- invlogit(par[1])
   r <- exp(par[2])
   alpha <- exp(par[3])
   beta <- exp(par[-(1:3)])
   list(p0 = p0, r = r, alpha = alpha, beta = beta)
}

curvefit <- function(mle, y, X){
  res <- reportpars(mle)
  fit <- cdf(res, X)
  data.table(ObsTrialPct = y, PredTrialPct = fit)
}

reportpars <- function(mle) {
  if (!inherits(mle, "try-error") && ("par" %in% names(mle))) {
    par <- mle$par
    return(list(p0 = invlogit(par[1]), r = exp(par[2]), alpha = exp(par[3]), beta = exp(par[-(1:3)])))
  } else {
    return(list())
  }   
}

curvefit <- function(mle, dat, xvars){  
  res <- reportpars(mle)
  if (length(res)) {
    y <- dat[, ObsTrialPct]
    X <- dat[, xvars, with = FALSE]
    X <- as.matrix(X)
    
    fit <- cdf(res, X)
    mape <- mean(100 * abs(fit/(y + 1e-5) - 1))
    betamat <- tcrossprod(rep(1, nrow(dat)), res$beta)
    colnames(betamat) <- xvars
    return(cbind(dat, PredTrialPct = fit, MAPE = mape, p0 = res$p0, r = res$r, alpha = res$alpha, beta = betamat))
  } else {
    return(data.table())
  }
}


mle_oneid(campdat2lst[[1]], c("ACV.MultiOutlet"))
mle.res1 <- mle_allids(campdat2lst, xvars = c("ACV.MultiOutlet"))

# Check convergence status
convcheck <- sapply(mle.res1, function(x) if (inherits(x, "try-error")) NA else x$convergence)
table(convcheck, useNA = "ifany")

# Get model fit
modelfit <- mapply(curvefit, mle.res1, campdat2lst, MoreArgs = list(xvars = "ACV.MultiOutlet"), SIMPLIFY = FALSE)
modelfit <- rbindlist(modelfit)
modelfit[,`:=`(id = NULL, acv.max = NULL, acv.range = NULL)]

# Get estimated parameters and MAPE
estpars <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("Product", "NPP_BRAND_NUM", "NPP_BRAND", "p0", "r", "alpha", "beta")]
mapes <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("MAPE")]
 
##########################################################################################
# Run MLE with two covariates
mle.res2 <- mle_allids(campdat2lst, xvars = c("ACV.MultiOutlet", "SCV_CPN_CIRC"))

# Check convergence status
convcheck2 <- sapply(mle.res2, function(x) if (inherits(x, "try-error")) NA else x$convergence)
table(convcheck2, useNA = "ifany")

# Get model fit
modelfit2 <- mapply(curvefit, mle.res2, campdat2lst, MoreArgs = list(xvars = c("ACV.MultiOutlet", "SCV_CPN_CIRC")), SIMPLIFY = FALSE)
modelfit2 <- rbindlist(modelfit2)
modelfit2[,`:=`(id = NULL, acv.max = NULL, acv.range = NULL)]

# Get estimated parameters and MAPE
estpars2 <- modelfit2[, lapply(.SD, unique), keyby = UPC, .SDcols = c("Product", "NPP_BRAND_NUM", "NPP_BRAND", "p0", "r", "alpha", "beta.ACV.MultiOutlet", "beta.SCV_CPN_CIRC")]
mapes2 <- modelfit2[, lapply(.SD, unique), keyby = UPC, .SDcols = c("MAPE")]

write.csv(modelfit, file = "UPC-level-mle-1.csv", row.names = FALSE)

save(modelfit, modelfit2, estpars, estpars2, mapes, mapes2, file = "model-res-1.rdata")



##### Test likelihood and mle with one UPC
y1 <- campdat2[UPC == "11030006421", ObsTrialPct]
X1 <- campdat2[UPC == "11030006421", ACV.MultiOutlet]
X1 <- matrix(X1, ncol = 1)
loglikelihood(c(p0 = invlogit(-1), r = exp(-1), alpha = exp(-1), beta = exp(-1)), y = y1, X = X1, 1)
mle.test <- optim(c(-1, -1, -1, -1), objfn, y = y1, X = X1, control = list(fnscale = -1, maxit = 20000))

##### Perform MLE on all UPCs
mle.test <- mle1(campdat2lst)


##### Create model fit dataset
modelfit <- vector("list", length(mle.test))
for (i in seq_along(mle.test)) {
  if (!is.na(convcheck[i])) {
    pars <- reportpars(mle.test[[i]])
    fit <- curvefit(mle.test[[i]], campda t2lst[[i]]$ObsTrialPct, matrix(campdat2lst[[i]]$ACV.MultiOutlet, ncol = 1))
    modelfit[[i]] <- cbind(campdat2lst[[i]], PredTrialPct = fit$PredTrialPct, p0 = pars$p0, r = pars$r, alpha = pars$alpha, beta = pars$beta)
  } else {
    modelfit[[i]] <- data.table()
  }
}

modelfit <- rbindlist(modelfit)
modelfit[,`:=`(id = NULL, acv.max = NULL, acv.range = NULL)]
setnames(modelfit, "beta", "beta.ACV.MultiOutlet")

##### Plot model fit
for (i in unique(modelfit$UPC)) {
  mat <- modelfit[UPC == i, as.matrix(list(ObsTrialPct, PredTrialPct))]
  matplot(mat, type = "l", main = paste("UPC:", i), ylab = "Cum. Trial Rate")
  readline("cont?")
}

write.csv(modelfit, file = "UPC-level-mle-1.csv", row.names = FALSE)

##### Get parameter estimates
estpars <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("Product", "NPP_BRAND_NUM", "NPP_BRAND", "p0", "r", "alpha", "beta.ACV.MultiOutlet")]

save(modelfit, estpars, file = "model-res-1.rdata")