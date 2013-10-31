tmp <-
function (Data, Prior, Mcmc) 
{
    pandterm = function(message) {
        stop(message, call. = FALSE)
    }
    if (missing(Data)) {
        pandterm("Requires Data argument -- list of y and X")
    }
    if (is.null(Data$X)) {
        pandterm("Requires Data element X")
    }
    X = Data$X
    if (is.null(Data$y)) {
        pandterm("Requires Data element y")
    }
    y = Data$y
    nvar = ncol(X)
    nobs = length(y)
    if (nobs != nrow(X)) {
        pandterm("length(y) ne nrow(X)")
    }
    if (missing(Prior)) {
        betabar = c(rep(0, nvar))
        A = 0.01 * diag(nvar)
        nu = 3
        ssq = var(y)
    }
    else {
        if (is.null(Prior$betabar)) {
            betabar = c(rep(0, nvar))
        }
        else {
            betabar = Prior$betabar
        }
        if (is.null(Prior$A)) {
            A = 0.01 * diag(nvar)
        }
        else {
            A = Prior$A
        }
        if (is.null(Prior$nu)) {
            nu = 3
        }
        else {
            nu = Prior$nu
        }
        if (is.null(Prior$ssq)) {
            ssq = var(y)
        }
        else {
            ssq = Prior$ssq
        }
    }
    if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
        pandterm(paste("bad dimensions for A", dim(A)))
    }
    if (length(betabar) != nvar) {
        pandterm(paste("betabar wrong length, length= ", length(betabar)))
    }
    if (missing(Mcmc)) {
        pandterm("requires Mcmc argument")
    }
    else {
        if (is.null(Mcmc$R)) {
            pandterm("requires Mcmc element R")
        }
        else {
            R = Mcmc$R
        }
        if (is.null(Mcmc$keep)) {
            keep = 1
        }
        else {
            keep = Mcmc$keep
        }
    }
    cat(" ", fill = TRUE)
    cat("Starting IID Sampler for Univariate Regression Model", 
        fill = TRUE)
    cat("  with ", nobs, " observations", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("Prior Parms: ", fill = TRUE)
    cat("betabar", fill = TRUE)
    print(betabar)
    cat("A", fill = TRUE)
    print(A)
    cat("nu = ", nu, " ssq= ", ssq, fill = TRUE)
    cat(" ", fill = TRUE)
    cat("MCMC parms: ", fill = TRUE)
    cat("R= ", R, " keep= ", keep, fill = TRUE)
    cat(" ", fill = TRUE)
    sigmasqdraw = double(floor(Mcmc$R/keep))
    betadraw = matrix(double(floor(Mcmc$R * nvar/keep)), ncol = nvar)
    itime = proc.time()[3]
    cat("IID Iteration (est time to end - min) ", fill = TRUE)
    fsh()
    for (rep in 1:Mcmc$R) {
        RA = chol(A)
        W = rbind(X, RA)
        z = c(y, as.vector(RA %*% betabar))
        IR = backsolve(chol(crossprod(W)), diag(nvar))
        btilde = crossprod(t(IR)) %*% crossprod(W, z)
        res = z - W %*% btilde
        s = t(res) %*% res
        sigmasq = (nu * ssq + s)/rchisq(1, nu + nobs)
        beta = btilde + as.vector(sqrt(sigmasq)) * IR %*% rnorm(nvar)
        if (rep%%100 == 0) {
            ctime = proc.time()[3]
            timetoend = ((ctime - itime)/rep) * (R - rep)
            cat(" ", rep, " (", round(timetoend/60, 1), ")", 
                fill = TRUE)
            fsh()
        }
        if (rep%%keep == 0) {
            mkeep = rep/keep
            betadraw[mkeep, ] = beta
            sigmasqdraw[mkeep] = sigmasq
        }
    }
    ctime = proc.time()[3]
    cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
        "\n")
    attributes(betadraw)$class = c("bayesm.mat", "mcmc")
    attributes(betadraw)$mcpar = c(1, R, keep)
    attributes(sigmasqdraw)$class = c("bayesm.mat", "mcmc")
    attributes(sigmasqdraw)$mcpar = c(1, R, keep)
    return(list(betadraw = betadraw, sigmasqdraw = sigmasqdraw))
}
