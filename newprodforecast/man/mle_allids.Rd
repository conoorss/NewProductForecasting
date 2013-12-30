\name{mle_allids}
\alias{mle_allids}
\title{Maximum Likelihood Estimation for Trial Model}
\usage{
  mle_allids(lst, useSF = FALSE, parallel = TRUE, cpus = 4,
    ...)
}
\arguments{
  \item{lst}{a list of data frames}

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

