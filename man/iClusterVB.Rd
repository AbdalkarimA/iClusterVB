\name{iClusterVB}
\alias{iClusterVB}
\title{A Next-Generation R Package for Bayesian Clustering}
\usage{
iClusterVB(mydata, dist, K = 10, initial_method = "VarSelLCM", VS_method = 0,
initial_cluster = NULL, initial_vs_prob = NULL, initial_fit = NULL,
initial_omega = NULL, input_hyper_parameters = NULL, max_iter = 200,
early_stop = 1, per = 10, convergence_threshold = 1e-4)
}
\description{
iClusterVB—an integrative clustering approach grounded in variational Bayesian inference. iClusterVB allows researchers to navigate high-dimensional datasets comprising mixed-type data, including continuous, categorical, and count variables, as well as multi-view data.
}
\examples{
}

\arguments{
\item{mydata}{input data - list of length R data set}
\item{dist}{type of data, or distribution: dist = 'continuous', 'multinomial', 'poisson' (vector of length R)}
\item{K}{maximum number of clusters; Default is 10}
\item{initial_method}{initialization method: VarSelLCM, random, k-prototype, k-means, mclust, or lca; Default is VarSelLCM}
\item{VS_method}{variable selection method: 0 = clustering without variable selection, 1 = clustering with variable selection; Default is 0}
\item{initial_cluster}{initial cluster membership, if it is not NULL, it will overwrite the previous initial values setting for this parameter; Default is NULL}
\item{initial_vs_prob}{initial variable selection probability, a scalar; Default is NULL}
\item{intitial_fit}{initial values based on previously fitted iClusterVB model (an iClusterVB object); Default is NULL}
\item{initial_omega}{customized initial values for variable inclusion probabilities, if it is not NULL, it will overwrite the previous initial values setting for this parameter, if VS_method = 1, initial_omega is a list with length R, each element of the list is an array with dim=c(N,p[[r]])), r = 1,...,R; Default is NULL}
\item{initial_hyper_paramters}{input hyper-parameters of the prior distributions for the model; Default is NULL}
\item{max_iter}{maximum iteration of the VB algorithm; Default is 200}
\item{early_stop}{whether to stop the algorithm when it converges or continue until it reaches max_iter: 1 - algorithm stops when converges, 0 - algorithm stops when it reaches the maximum iteration (regardless of whether it converges or not); Default is 1}
\item{per}{print information every "per" iteration; Default is 10}
\item{convergence_threshold}{define a convergence threshold for the change in ELBO; Default is 0.0001}

}
