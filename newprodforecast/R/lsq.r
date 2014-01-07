ols_trialmodel <- function(dat, mfcall, startvals, tmSpec, method, optimControl) {
	# Internal function for setting up model frame and calling optim
	#cat(unique(dat$UPC), fill = TRUE)
	m <- match(c("formula", "data", "subset", "na.action"), names(mfcall), 0L)
	mfcall <- mfcall[c(1L, m)]
	mfcall[[1L]] <- as.name("model.frame")
	mfcall$data <- quote(dat)
	mf <- as.data.table(eval(mfcall))
	mt <- attr(mf, "terms")
	y <- stats::model.response(mf, "numeric")
	x <- stats::model.matrix(mt, mf)
	if (!is.null(tmSpec$acvMultiplier)) 
		acv <- dat[, tmSpec$acvMultiplier, with = FALSE]
	else
		acv <- 1
	splitFun <- generate_vec_splitter(tmSpec)
	transFuns <- generate_transforms(tmSpec)
	call_sumsqerr <- function(params) {
		#print(params)
		params <- splitFun(params)
		params <- mapply(function(x, fun) fun(x), params, transFuns, SIMPLIFY = FALSE)
		sumsqerr_trialmodel(params, y, x, acv)
	}
	#res <- try(optim(startvals, sumsqerr_trialmodel, y = y, x = x, acv = acv, tmSpec = tmSpec, method = method, control = optimControl), silent = TRUE)
	#structure(res, y = y, x = x, terms = mt)
	res <- try(optim(startvals, call_sumsqerr, method = method, control = optimControl), silent = TRUE)
	if (inherits(res, "try-error")) {
		coef <- list()
	} else {
		coef <- splitFun(res$par)
		coef <- mapply(function(x, fun) fun(x), coef, transFuns, SIMPLIFY = FALSE)
		pred <- fit_trialmodel(coef, y, x, acv)
		err <- pred - y
		mape <- median(abs(err) / y)
	}
	list(coefficients = coef, optres = res, residuals = err, fit = pred, mape = mape)
	#    list(coefficients = transformTrialParams(res$par, tmSpec), optres = res)
}

if (1) {
	sumsqerr_trialmodel <- function(params, y, x, acv) {
		# Internal function to compute squared error
		predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
		predictor <- cumsum(predictor)
		tpredictor <- cumulative_curve(params, predictor, acv)
		sum((y - tpredictor)^2)
	}

} else {

	sumsqerr_trialmodel <- function(params, y, x, acv, tmSpec) {
		# Internal function to compute squared error
		params <- transformTrialParams(params, tmSpec)
		predictor <- if (!is.null(params$beta)) drop(exp(x %*% params$beta)) else rep(1, NROW(x))
		predictor <- cumsum(predictor)
		tpredictor <- cumulative_curve(params, predictor, acv)
		sum((y - tpredictor)^2)
	}
}
