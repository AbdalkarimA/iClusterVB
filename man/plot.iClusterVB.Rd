% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classMethods.R
\name{plot.iClusterVB}
\alias{plot.iClusterVB}
\title{Generic plot method for `iClusterVB` objects}
\usage{
\method{plot}{iClusterVB}(x, ...)
}
\arguments{
\item{x}{A fitted iClusterVB object.}

\item{...}{Potential further arguments (unused)}
}
\value{
Returns an evidence lower bound (ELBO) plot and a barplot of cluster
  percentages.
}
\description{
Generic plot method for `iClusterVB` objects
}
\examples{
# Setting up the data
dat1 <- list(
  gauss_1 = sim_data$continuous1_data[c(1:20, 61:80, 121:140, 181:200), 1:75],
  gauss_2 = sim_data$continuous2_data[c(1:20, 61:80, 121:140, 181:200), 1:75],
  poisson_1 = sim_data$count_data[c(1:20, 61:80, 121:140, 181:200), 1:75])


dist <- c(
  "gaussian", "gaussian",
  "poisson")

fit_iClusterVB <- iClusterVB(
  mydata = dat1,
  dist = dist,
  K = 4,
  initial_method = "VarSelLCM",
  VS_method = 1,
  max_iter = 25
)

plot(fit_iClusterVB)

}
