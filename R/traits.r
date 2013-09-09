#' Traits for defining new product forecasting models
#
TrialModel <- 
	proto(expr = 
		  {
			  loglikVector <- function(., params = .$params, time = .$data$time, count = .$data$count, sampleSize = sum(.$data$count)) {
				  time <- na.omit(time)
				  count <- count[-na.action(time)]
				  tmax <- max(time)
				  tmin <- min(time) 
				  if (tmin < 1)
					  stop("Minimum time must be at least 1")
				  cdf1 <- .$trialTimeCdf(time, params)
				  cdf2 <- .$trialTimeCdf(time - 1, params)
				  cdf3 <- .$trialTimeCdf(tmax, params)

				  ll <- count * log(cdf1 - cdf2) + (sampleSize - sum(count)) * log(1 - cdf3)
				  ll
			  }
			  
			  loglik <- function(., params, ..., trace = 0L) {
				  llv <- .$loglikVector(params, ...)
				  ll <- sum(llv)
				  if (trace == 1L) {
					  if (.$iter_count %% 10 == 0)
						  cat("Likelihood: ", ll, "\n")
				  } else if (trace == 2L) {
					  .$iter_count <- .$iter_count + 1L
					  elapsed_time <- Sys.time() - .$est_start_time
					  if (.$iter_count %% 10 == 0)
						  cat("Iterations: ", .$iter_count, " Time Elapsed: ",
							  format(unclass(elapsed_time), digits = 2), " ", 
							  attr(elapsed_time, "units"), " Likelihood: ", ll, "\n")
				  }
			  }

			  transformLogPars <- function(pars) {
				  m1 <- na.omit(match(c("r", "alpha", "lambda"), pars))
				  if (length(m1))
					  pars[m1] <- exp(pars[m1])
				  else
					  stop("Need one of the parameters \"r\", \"alpha\", \"lambda\"")
					  
				  if ("p" %in% names(pars))
					  pars[["p"]] <- exp(pars[["p"]]) / (1 + exp(pars[["p"]]))
					
			  }

			  checkData <- function(data) {
				  if (any(c("times","counts") %notin% names(data)))
					  stop("The required columns \"times\" and \"counts\" are not present in the data")
			  }

			  checkLambda <- function(params) {
				  if ("lambda" %notin% names(params))
					  stop("The parameter vector must contain an element named \"lambda\"") 
				  if (params[["lambda"]] <= 0)
					  stop("\"lambda\" must be positive")
			  }

			  checkGamPars <- function(params) {
				  if (any(c("r", "alpha") %notin% names(params)))
					  stop("The parameter vector must contain elements named \"r\" and \"alpha\"")
				  if (params[["r"]] <= 0 || params[["alpha"]] <= 0)
					  stop("\"r\" and \"alpha\" must be positive")
			  }

			  checkP <- function(params) {
				  if ("p" %notin% names(params))
					  stop("The parameter vector must contain elements named \"p\"")
				  if (params[["p"]] < 0 || params[["p"]] > 1)
					  stop("\"p\" must be between 0 and 1")
			  }


			  mle <- function(., method = "BFGS") {
				  with(., { est_start_time <- Sys.time(); iter_count <- 0L } )

				  eloglik <- function(pars, ...) {
					  pars <- transformLogPars(pars)
					  .$loglik(pars, ...)
				  }
				  pars <- log(.$startParams)
				  .$optRes <- optim(pars0, eloglik, trace = 2L, method = method, control = list(fnscale = -1, maxit = 10000))
				  .$params <- transformLogPars(.$optRes$par)
				  return(.$params) 
			  }

})


ExpTrial <- TrialModel$proto()
ExpTrial$simulate <- function(., sampleSize, numPeriods, params, seed = 1234) {
	if (sampleSize < 0 || numPeriods < 0)
		stop("\"sampleSize\" and \"numPeriods\" must be positive")
	checkLambda(params)
	set.seed(seed)
	trialTimes <- rexp(sampleSize, params[["lambda"]])
	icTrialTimes <- cut(trialTimes, breaks = 0:numPeriods)
	counts <- tapply(seq(sampleSize), icTrialTimes, length)
	.$data <- data.table(times = as.integer(levels(icTrialTimes)), counts = counts)
	.$simPars <- list(sampleSize = sampleSize, numPeriods = numPeriods, params = params)
}

ExpTrial$new <- function(., data = data.table(times = integer(0), counts = integer(0)), startParams = c(lambda = 1), sim = FALSE) {
	if (sim)
		message("Data is initialized to an empty data.table.\nUse \"simulate\" method to generate data.")
	checkData(data)
	checkLambda(startParams)
	.$proto(data = data, startParams = startParams)
}
ExpTrial$trialTimeCdf <- function(., time, params) {
	pexp(time, params[["lambda"]])
}

ExpNeverTrial <- ExpTrial$proto()
ExpNeverTrial$simulate <- function(., sampleSize, numPeriods, params, seed = 1234) {
	if (sampleSize < 0 || numPeriods < 0)
		stop("\"sampleSize\" and \"numPeriods\" must be positive")
	checkLambda(params)
	checkP(params)
	set.seed(seed)
	everTried <- rbinom(sampleSize, 1, params[["p"]])
	trialTimes <- rep(Inf, sampleSize)
	trialTimes[which(everTried == 1L)] <- rexp(sum(everTried), params[["lambda"]])
	icTrialTimes <- cut(trialTimes, breaks = 0:numPeriods)
	counts <- tapply(seq(sampleSize), icTrialTimes, length)
	.$data <- data.table(times = as.integer(levels(icTrialTimes)), counts = counts)
	.$simPars <- list(sampleSize = sampleSize, numPeriods = numPeriods, params = params)
}

ExpNeverTrial$new <- function(., data, startParams = c(lambda = 1, p = 1), sim = FALSE) {
	checkP(startParams)
	.super$new(., data, startParams)
}

ExpNeverTrial$trialTimeCdf <- function(., time, params) {
	params[["p"]] * .$.super$trialTimeCdf(time, params)
}

ExpGammaTrial <- TrialModel$proto()
ExpGammaTrial$new <- function(., data, startParams = c(r = 1, alpha = 1), sim = FALSE) {
	if (sim)
		message("Data is initialized to an empty data.table.\nUse \"simulate\" method to generate data.")
	checkData(data)
	checkGamPars(startParams)
	.$proto(data = data, startParams = startParams)
}
ExpGammaTrial$trialTimeCdf <- function(., time, params) {
	1 - (params[["alpha"]] / (params[["alpha"]] + time))^params[["r"]]
}

ExpGammaNeverTrial <- ExpGammaTrial$proto()
ExpGammaNeverTrial$new <- function(., data, startParams = c(r = 1, alpha = 1, p = 1), sim = FALSE) {
	checkP(startParams)
	.super$new(., data, startParams)
}
ExpGammaNeverTrial$trialTimeCdf <- function(., time, params) {
	params[["p"]] * .$super$trialTimeCdf(time, params)
}


