#' Outputs a summary based on an iClusterVB object
#'
#' @param object A fitted iClusterVB object.
#' @param vs_prob Input to count the number of variables that meet he minimum probability of inclusion vs_prob. Default is 0.5. Only works for VS_method = 1.
#' @return Returns a summary showing the initialization method, the cluster memberships, and the number of variables for which the probabilty of inclusion is above vs_prob
#' @examples
#' summary(fit_iClusterVB)


summary.iClusterVB <- function(object, vs_prob = NULL, ...) {

  if(is.null(vs_prob)) {
    vs_prob <- 0.5
  }

  mybiglist <- list()

  for (i in 1:length(object$mydata)) {
    name <- paste(object$dist[i])
    mybiglist[[name]] <- sum(object$model_parameters$rho[[i]] >= vs_prob)
  }


  summ <- list("Initialization Method" = object$initial_values$initial_method,
               "Cluster Memberships" = table(object$cluster),
               "# Selected variables above vs_prob" = mybiglist )
  summ
}
