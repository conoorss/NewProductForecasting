\name{trialmodel.mle}
\alias{trialmodel.mle}
\title{Maximum Likelihood Estimation for Trial Model}
\usage{
  trialmodel.mle(formula, data, group, startvals = list(),
    tmControl = list(), sf = FALSE, sfControl = list(),
    method = "BFGS", optimControl = list())
}
\arguments{
  \item{model}{is a list with two or three components. The
  first component is the dependent variable, the second
  component is the independent variable}

  \item{useSF}{logical flag which turns on usage of the
  \code{snowfall} package if \code{TRUE}. Default is
  \code{FALSE}}

  \item{parallel}{argument to be passed onto
  \code{snowfall} if useSF is TRUE. The default is TRUE.}

  \item{cpus}{argument to be passed onto \code{snowfall} if
  useSF is TRUE. The default is to use 4 processors.}

  \item{...}{are other arguments to be passed onto maximum
  likelihood estimation for each unit.}
}
\value{
  A list of lists with each sublist containing either the
  results from the estimation or an object of class
  \code{try-error} if maximum likelihood fails
}
\description{
  runs maximum likelihood estimation for each unit in the
  data and returns a list of lists
}

