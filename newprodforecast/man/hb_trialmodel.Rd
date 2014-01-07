\name{hb_trialmodel}
\alias{hb_trialmodel}
\title{Hierarchical Bayesian Estimation of the Trial Model}
\usage{
  hb_trialmodel(Data, Priors, Mcmc)
}
\arguments{
  \item{Data}{is a list with components yList, xList, and
  z. yList and xList are lists of vectors and matrices
  respectively, representing the dependent variable and the
  covariates. Each element of the list represents a group
  in the lower level of the hierarchical model. z is a
  matrix of covariates for the upper level of the
  hierarchy.}

  \item{Priors}{is a list with components Deltabar, A, nu,
  V, sigsq_mean, and sigsq_sd. The components Deltabar, A,
  nu, and V are parameters of the hyperprior in the
  multivariate regression for the upper level model (LINK
  TO bayesm::rmultireg). sigsq_mean and sigsq_sd represent
  the mean and standard deviation of the univariate normal
  prior on the log of the variance of the observation error
  term (sigma^2).}

  \item{Mcmc}{is a list with starting values of model
  parameters, scaling parameters for the Random Walk
  Metropolis algorithm, and options controlling various
  aspects of the MCMC sampler. The following are components
  of the Mcmc list:}
}
\value{
  A list of lists with components param_samples,
  prior_samples, rejections, loglikelihood, logprior.
  param_samples is a list of arrays/matrices storing
  samples of parameters from the lower model. prior_samples
  is a list of arrays/matrices storing parameters from the
  upper model. loglikelihood is a matrix of log likelihood
  values. logprior is a list of matrices storing the log
  priors for lower model parameters.
}
\description{
  Fits a hierarchical bayesian trial model
}

