% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_data}
\alias{sim_data}
\title{Simulated Dataset}
\format{
A list containing four datasets, and other elements of interest.
}
\usage{
data(sim_data)
}
\description{
The dataset consists of \eqn{N = 240} individuals and \eqn{R =
  4} data views with different data types. Two of the data views are
  continuous, one is count, and one is binary. The \emph{true} number of
  clusters was set to \eqn{K = 4}, and the cluster proportions were set at \eqn{\pi_1
  = 0.25, \pi_2 = 0.25, \pi_3 = 0.25, \pi_4 = 0.25}, such that we have
  balanced cluster proportions. Each of the data views had \eqn{p_r = 500}
  features, \eqn{r = 1, \dots, 4}, but only 50, or 10\%, were relevant
  features that contributed to the clustering, and the rest were noise
  features that did not contribute to the clustering. In total, there were
  \eqn{p = \sum_{r=1}^4 = 2000} features.

  For data view 1 (continuous), relevant features were generated from the
  following normal distributions: \eqn{\text{N}(10, 1)} for Cluster 1,
  \eqn{\text{N}(5, 1)} for Cluster 2, \eqn{\text{N}(-5, 1)} for Cluster 3,
  and \eqn{\text{N}(-10, 1)} for Cluster 4, while noise features were
  generated from \eqn{\text{N}(0, 1)}. For data view 2 (continuous), relevant
  features were generated from the following normal distributions:
  \eqn{\text{N}(-10, 1)} for Cluster 1,    \eqn{\text{N}(-5, 1)} for Cluster
  2, \eqn{\text{N}(5, 1)} for Cluster 3, and \eqn{\text{N}(10, 1)} for
  Cluster 4, while noise features were generated from \eqn{\text{N}(0, 1)}.
  For data view 3 (binary), relevant features were generated from the
  following Bernoulli distributions: \eqn{\text{Bernoulli}(0.05)} for Cluster
  1,    \eqn{\text{Bernoulli}(0.2)} for Cluster 2,
  \eqn{\text{Bernoulli}(0.4)} for Cluster 3, and \eqn{\text{Bernoulli}(0.6)}
  for Cluster 4, while noise features were generated from
  \eqn{\text{Bernoulli}(0.1)}. For data view 4 (count), relevant features
  were generated from the following Poisson distributions:
  \eqn{\text{Poisson}(50)} for Cluster 1, \eqn{\text{Poisson}(35)} for
  Cluster 2, \eqn{\text{Poisson}(20)} for Cluster 3, and
  \eqn{\text{Poisson}(10)} for Cluster 4, while noise features were generated
  from \eqn{\text{Poisson}(2)}.
}
\keyword{datasets}
