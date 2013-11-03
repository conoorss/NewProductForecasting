#
# Model Spec:
#  y_{jt} = p0_{j} * (1 - (alpha_{j}/(alpha_{j} + A(t)))^r_{j})
#  A(t) = \sum_{l=1}^{t} exp(x_l * beta)
#
#	logit(p0) = N(p0bar, sigmasq_p0)
#	log(alpha) ~ N(alphabar, sigmasq_alpha)
#	log(r) ~ N(rbar, sigmasq_r)
# log(beta) ~ N(betabar, Sigma_beta)
#

invlogit <- function(x) exp(x) / (1 + exp(x))

simExpGammaTrial <- function(sample_size, horizon, simPars) {
  require(data.table)
  require(MASS)
  p0 <- with(simPars, rnorm(sample_size, p0bar, sqrt(sigmasq_p0)))
  r <- with(simPars, rnorm(sample_size, rbar, sqrt(sigmasq_r)))
  alpha <- with(simPars, rnorm(sample_size, alphabar, sqrt(sigmasq_alpha)))
  beta <- with(simPars, mvrnorm(sample_size, betabar, Sigma_beta))
  
  p0 <- invlogit(p0)
  r <- exp(r)
  alpha <- exp(alpha)
  beta <- exp(beta)
  
  data <- vector("list", sample_size)
  time <- seq(horizon)
  for (i in seq(sample_size)) {
    X <- matrix(runif(horizon * ncol(beta)), horizon, ncol(beta))
    A <- cumsum(exp(X %*% beta[i,]))
    y <- p0[i] * (1 - (alpha[i] / (alpha[i] + A))^r[i])
    data[[i]] <- data.table(id = i, time = time, y = y, X = X, A = A)
  }
  data <- rbindlist(data)
  list(sample_size = sample_size,
       numCovs = ncol(beta),
       potential = rep(1, sample_size),
       horizon = horizon, 
       simPars = simPars,
       indPars = list(p0 = p0, r = r, alpha = alpha, beta = beta), 
       data = data)
}

plotSimData <- function(lst) {
  mat <- matrix(lst$data$y, nrow = lst$horizon)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  layout(matrix(c(1,1,1, 2, 3, 4), 2, 3, byrow = TRUE))
  ix <- seq(ncol(mat)) %% 10 == 0
  matplot(mat[,ix], type = "l", col = 1, main = "Trial curves")
  plot(density(lst$indPars$p0), xlab = "", ylab = "", main = "Density of p0")
  plot(density(lst$indPars$r), xlab = "", ylab = "", main = "Density of r")
  plot(density(lst$indPars$alpha), xlab = "", ylab = "", main = "Density of alpha")
  #invisible()
}

cdf <- function(pars, X, mod = NULL) {
  p0 <- pars[["p0"]]
  r <- pars[["r"]]
  alpha <- pars[["alpha"]]
  beta <- pars[["beta"]]
  
  A <- exp(X %*% beta)
  
  if (is.null(mod))
    A <- cumsum(A)
  else if (mod == "lag")
    A <- cumsum(c(0, head(A, -1)))
  else 
    A <- tail(cumsum(A), 1)
  
  p0 * (1 - (alpha / (alpha + A))^r)
}

loglikelihood <- function(pars, y, X, potential) {	
  incy <- diff(c(0, y))
  cdf1 <- cdf(pars, X)
  cdf2 <- cdf(pars, X, "lag")
  cdf3 <- cdf(pars, X, "max")
  
  ll <- sum(incy * log(cdf1 - cdf2)) + (potential - sum(incy)) * log(1 - cdf3)
  ll
}

