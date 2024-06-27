#' Generic summary method for `iClusterVB` objects
#'
#' @param object A fitted iClusterVB object.
#' @param rho The minimum posterior inclusion probability of interest to count
#'   the number of variables that >= it. Default is 0.5. Only works for
#'   VS_method = 1.
#' @param ... Potential further arguments
#' @return Returns a summary list for an `agnes` object.
#' @examples
#' summary(fit_iClusterVB)
#'
#' @export summary.iClusterVB
#' @method summary iClusterVB
#' @useDynLib iClusterVB, .registration=TRUE


summary.iClusterVB <- function(object, rho = 0.5, ...) {


  cat("Total number of individuals:\n")
  print(object$data_summary$N)

  cat("\n")

  if(!is.null(object$model_parameters$rho)) {
    for (i in 1:length(object$mydata)) {
      name <- paste("View", i, "-", object$dist[i], sep = " ")
      cat(paste("# of variables above the posterior inclusion probability of", rho, "for View", i, "-", object$dist[i], sep = " "))
      cat("\n")
      print(paste(sum(object$model_parameters$rho[[i]] >= rho), "out of a total of", object$data_summary$p[[i]], sep = " "))
      cat("\n")
    }
  }


  cat(paste("User-inputted maximum number of clusters:", object$K, sep = " "))
  cat("\n")
  cat(paste("Number of clusters determined by algorithm:", length(unique(object$cluster))))

  cat("\n")
  cat("\n")

  cat("Cluster Membership:")

  print(table(object$cluster))

}

