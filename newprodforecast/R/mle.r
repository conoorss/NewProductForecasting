
mle_trialmodel <- function(dat, mfcall, startvals, tmSpec, method, optimControl) {
	# Internal function for setting up model frame and calling optim
	#cat(unique(dat$UPC), fill = TRUE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
	mfcall <- mfcall[c(1L, m)]
	mfcall[[1L]] <- as.name("model.frame")
	mfcall$data <- quote(dat)
	mf <- as.data.table(eval(mfcall))
	mt <- attr(mf, "terms")
	y <- model.response(mf, "numeric")
	incy <- diff(c(0, y))
	x <- model.matrix(mt, mf)
	if (!is.null(tmSpec$acvMultiplier)) 
		acv <- dat[, tmSpec$acvMultiplier, with = FALSE]
	else
		acv <- 1
	splitFun <- generate_vec_splitter(tmSpec)
	transFuns <- generate_transforms(tmSpec)
	call_loglik <- function(params) {	  
		#print(params)
		params <- splitFun(params)
		params <- apply_transforms(params, transFuns)
		loglik_trialmodel(params, incy, x, acv)
	}
	#res <- try(optim(startvals, loglik_trialmodel, y = y, x = x, acv = acv, tmSpec = tmSpec, method = method, control = optimControl), silent = TRUE)
	#structure(res, y = y, x = x, terms = mt)
	res <- try(optim(startvals, call_loglik, method = method, control = optimControl), silent = TRUE)
	if (inherits(res, "try-error")) {
		coef <- list()
		err <- pred <- mape <- numeric()
	} else {
		coef <- splitFun(res$par)
		coef <- apply_transforms(coef, transFuns)
		pred <- fit_trialmodel(coef, y, x, acv)
		err <- pred - y
		mape <- median(abs(err) / y)
	}
	list(coefficients = coef, optres = res, residuals = err, fit = pred, mape = mape)
}

fit_trialmodel <- function(params, y, x, acv) {
	# Assumes that params is a named list
	predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
	predictor <- cumsum(predictor)
	tpredictor <- cumulative_curve(params, predictor, acv)
	tpredictor
}

if (1) {
	loglik_trialmodel <- function(params, incy, x, acv) {
		# Internal function to compute negative loglikelihood
		# This version assumes params is a named list of constrained parameters
		predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
		predictor <- cumsum(c(0, predictor))
		tpredictor <- cumulative_curve(params, predictor, acv)
		cdf1 <- tail(tpredictor, -1)
		cdf2 <- head(tpredictor, -1)
		cdf3 <- tail(tpredictor, 1)
		ll <- sum(incy * log(cdf1 - cdf2)) + (1 - sum(incy)) * log(1 - cdf3)
		ll
	}



} else {

	loglik_trialmodel <- function(params, y, x, acv, tmSpec) {
		# Internal function to compute negative loglikelihood
		# Older version: assumes params is a vector 

		params <- transformTrialParams(params, tmSpec)
		incy <- diff(c(0, y))
		predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
		predictor <- cumsum(c(0, predictor))
		#  tpredictor <- cumulative_curve(params, tmSpec, predictor, acv)
		tpredictor <- cumulative_curve(params, predictor, acv)
		cdf1 <- tail(tpredictor, -1)
		cdf2 <- head(tpredictor, -1)
		cdf3 <- tail(tpredictor, 1)
		ll <- sum(incy * log(cdf1 - cdf2)) + (1 - sum(incy)) * log(1 - cdf3)
		ll
	}
}

