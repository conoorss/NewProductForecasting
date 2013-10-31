source("try-3.r")
set.seed(10000)

if (!file.exists("simres-3.rdata")) {
  simres <- simExpGammaTrial(1000, 200, 
                             list(p0bar = -2, sigmasq_p0 = 0.5, 
                                  rbar = -3.5, sigmasq_r = 0.5, 
                                  alphabar = 3.5, sigmasq_alpha = 0.5, 
                                  betabar = matrix(0.5, 1, 1), Sigma_beta = 0.5 * diag(1)))
                                  #betabar = matrix(c(0.5, 0.3), 2, 1), Sigma_beta = 0.5 * diag(2)))
  save(simres, file = "simres-3.rdata")
} else {
  load("simres-3.rdata")
}

# Test the likelihood function

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


mle1 <- function(lst){
  res <- vector("list", lst$sample_size)
  for (i in seq(lst$sample_size)) {
    if (i %% 10 == 0) cat(i, " done \n")
    yi <- lst$data[id == i, y]
    sel <- grep("^X", names(lst$data), value = TRUE)
    Xi <- lst$data[id == i, sel, with = FALSE]
    Xi <- as.matrix(Xi)
    timei <- lst$data[id == i, time]
    res[[i]] <- try(optim(c(-1, -1, -1, -1), objfn, y = yi, X = Xi, control = list(fnscale = -1)))
  }
  res
}


# Tests for first id
p01 <- simres$indPars$p0[1]
r1 <- simres$indPars$r[1]
alpha1 <- simres$indPars$alpha[1]
beta1 <- simres$indPars$beta[1,]

y1 <- simres$data[id == 1, y]
time1 <- simres$data[id == 1, time]
X1 <- simres$data[id == 1, grep("^X", names(simres$data), value = TRUE), with = FALSE]
X1 <- as.matrix(X1)

loglikelihood(list(p0 = p01, r = r1, alpha = alpha1, beta = beta1), y = y1, X = X1, 1)
mle.test <- optim(c(-1, -1, -1, -1), objfn, y = y1, X = X1, control = list(fnscale = -1))

#### NOT YET DONE ####

if (FALSE) {
mle.test <- mle1(simres)
conv <- do.call("c", lapply(mle.test, function(x) x$convergence))
estpars <- do.call("rbind", lapply(mle.test, function(x) x$par))
estpars[,-1] <- exp(estpars[,-1])
estpars[,1] <- invlogit(estpars[,1])
colnames(estpars) <- c("p0.hat", "r.hat", "alpha.hat")

fitplot <- function(x, y, param) { 
  plot(x, y, main = paste("Fit of", param), xlab = "Obs", ylab = "Est", pch = ".", cex = 4)
  abline(a = 0, b = 1, col = 2)
}

opar <- par(no.readonly = TRUE)
par(mfrow = c(3,1), mar = c(4.1, 4.1, 4.1, 2.1)) 
fitplot(simres$indPars$p0, estpars[,"p0.hat"], "p0")
fitplot(simres$indPars$r, estpars[,"r.hat"], "r")
fitplot(simres$indPars$alpha, estpars[,"alpha.hat"], "alpha")
par(opar)
}