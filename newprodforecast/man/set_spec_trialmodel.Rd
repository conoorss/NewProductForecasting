\name{set_spec_trialmodel}
\alias{set_spec_trialmodel}
\title{Functions for checking/setting the specification of the trial/repeat models}
\usage{
  set_spec_trialmodel(lst)
}
\arguments{
  \item{lst}{A list with the specification of the model}
}
\value{
  A list with components "family", "p0", "acvMultiplier",
  "numCovariates
}
\description{
  Used to specify options for the trial model: functional
  family of the cumulative trial curve presence of
  asymptote parameter presence of multiplicative term for
  ACV (distribution)
}

