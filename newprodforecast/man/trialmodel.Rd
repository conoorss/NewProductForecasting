\name{trialmodel}
\alias{trialmodel}
\title{Estimation for Trial Model}
\usage{
  trialmodel(formula, data, group, startvals = numeric(),
    tmSpec = list(), sf = FALSE, sfControl = list(),
    estimation = c("MLE", "OLS"), method = "Nelder-Mead",
    optimControl = list(maxit = 20000))
}
\arguments{
  \item{formula}{is a two-sided formula object which
  specifies the dependent variable and the covariates if
  any.}

  \item{data}{is either a \code{data.frame} or
  \code{data.table} with the variables needed for the
  model.}

  \item{group}{is a string with the name of the group
  variable.}

  \item{startvals}{is a numeric vector with the starting
  values for the model.}

  \item{tmSpec}{is a list which specifies the model
  specification (family, p0, acvMultiplier, number of
  covariates).}

  \item{sf}{is a logical flag for usage of the
  \code{snowfall} package for parallelization (currently
  not implemented)}

  \item{sfControl}{is a list of control parameters for
  \code{snowfall}}

  \item{estimation}{is a string which is either "MLE" or
  "OLS"}

  \item{method}{is a string which specifies the
  optimization method}

  \item{optimControl}{is a list of control parameters for
  \code{optim}}
}
\value{
  A list of lists with each sublist containing either the
  results from the estimation or an object of class
  \code{try-error} if maximum likelihood fails
}
\description{
  Fits a trial model using either maximum likelihood or
  least squares. Returns an object of class
  \code{trialmodel}
}

