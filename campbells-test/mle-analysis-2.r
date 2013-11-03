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

mle1 <- function(lst){
  res <- vector("list", length(lst))
  for (i in seq_along(lst)) {
    if (i %% 10 == 0) cat(i, " done \n")
    yi <- lst[[i]][, ObsTrialPct]
    timei <- lst[[i]][, WEEK]
    res[[i]] <- try(optim(c(-1, -1, -1), objfn, y = yi, time = timei, control = list(fnscale = -1, maxit = 20000), method = "BFGS"))
  }
  res
}

reportpars <- function(mle) {
   par <- mle$par
   p0 <- invlogit(par[1])
   r <- exp(par[2])
   alpha <- exp(par[3])
   list(p0 = p0, r = r, alpha = alpha)
}

curvefit <- function(mle, y, time){
  res <- reportpars(mle)
  fit <- cdf(res, time)
  data.table(ObsTrialPct = y, PredTrialPct = fit)
}

##### Test likelihood and mle with one UPC
y1 <- campdat2lst[[1]][, ObsTrialPct]
time1 <- seq_along(y1)
loglikelihood(c(p0 = invlogit(-1), r = exp(-1), alpha = exp(-1)), y = y1, time = time1, 1)
mle.test <- optim(c(-1, -1, -1), objfn, y = y1, time = time1, control = list(fnscale = -1, maxit = 20000))

#### Perform MLE on all UPCs
mle.test <- mle1(campdat2lst)

##### Store convergence status
convcheck <- sapply(mle.test, function(x) if (inherits(x, "try-error")) NA else x$convergence)
table(convcheck)

##### Create model fit dataset
modelfit <- vector("list", length(mle.test))
for (i in seq_along(mle.test)) {
  if (!is.na(convcheck[i])) {
    pars <- reportpars(mle.test[[i]])
    fit <- curvefit(mle.test[[i]], campdat2lst[[i]]$ObsTrialPct, seq_along(campdat2lst[[i]]$ObsTrialPct))
    modelfit[[i]] <- cbind(campdat2lst[[i]], PredTrialPct = fit$PredTrialPct, p0 = pars$p0, r = pars$r, alpha = pars$alpha)
  } else {
    modelfit[[i]] <- data.table()
  }
}

modelfit <- rbindlist(modelfit)
modelfit[,`:=`(id = NULL, acv.max = NULL, acv.range = NULL)]

##### Plot model fit
for (i in unique(modelfit$UPC)) {
  mat <- modelfit[UPC == i, as.matrix(list(ObsTrialPct, PredTrialPct))]
  matplot(mat, type = "l", main = paste("UPC:", i), ylab = "Cum. Trial Rate")
  readline("cont?")
}

write.csv(modelfit, file = "UPC-level-mle-2.csv", row.names = FALSE)

##### Get parameter estimates
estpars <- modelfit[, lapply(.SD, unique), keyby = UPC, .SDcols = c("Product", "NPP_BRAND_NUM", "NPP_BRAND", "p0", "r", "alpha")]

save(modelfit, estpars, file = "model-res-2.rdata")