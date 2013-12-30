\name{tm.control}
\alias{tm.control}
\title{Specification for Trial Model}
\usage{
  tm.control(family = "exponential-gamma", p0 = FALSE,
    acv_multiplier = NULL)
}
\arguments{
  \item{family}{string with value "exponential" or
  "exponential-gamma". Default value is
  "exponential-gamma"}

  \item{p0}{logical. If the value is TRUE then the model
  includes an asymptote parameter p0. Default value is
  FALSE.}

  \item{acv_multiplier}{string giving the column name of
  the variable which is to be used in the multiplicative
  term for distribution.}
}
\value{
  A list with named elements "family", "p0",
  "acv_multiplier"
}
\description{
  Used to specify options for the trial model: functional
  family of the cumulative trial curve presence of
  asymptote parameter presence of multiplicative term for
  ACV (distribution)
}

