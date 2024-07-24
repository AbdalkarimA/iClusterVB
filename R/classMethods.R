#' Generic summary method for `iClusterVB` objects
#'
#' @param object A fitted iClusterVB object.
#' @param rho The minimum posterior inclusion probability of interest to count
#'   the number of features that are >= \code{rho}. Default is 0.5. Only works
#'   for VS_method = 1.
#' @param ... Potential further arguments
#' @return Returns a summary list for an `agnes` object.
#' @examples
#' ## S3 method for class 'iClusterVB'
#' \dontrun{summary(fit_iClusterVB, rho = 0.5)}
#'
#' @export
#' @method summary iClusterVB
#' @useDynLib iClusterVB, .registration=TRUE


summary.iClusterVB <- function(object, rho = 0.5, ...) {

  # args <- list(...)
  #
  # if (!("rho" %in% names(args))) {
  #   rho <- 0.5
  # }


  cat("Total number of individuals:\n")
  print(object$data_summary$N)

  cat("\n")


  cat(paste("User-inputted maximum number of clusters:", object$K, sep = " "))
  cat("\n")
  cat(paste("Number of clusters determined by algorithm:", length(unique(object$cluster))))

  cat("\n")
  cat("\n")

  cat("Cluster Membership:")

  print(table(object$cluster))

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

}








#' Generic plot method for `iClusterVB` objects
#'
#' @param x A fitted iClusterVB object.
#' @param ... Potential further arguments (unused)
#'
#' @return Returns an evidence lower bound (ELBO) plot and a barplot of cluster
#'   percentages.
#'
#'
#' @export
#' @method plot iClusterVB
#' @examples
#' \dontrun{plot(fit_iClusterVB)}


plot.iClusterVB <- function(x, ...) {

  fit <- x

  rm(x)

  plot(x = 1:fit$iter,
       y = fit$elbo[1:1:fit$iter],
       type = "o",
       lwd = 2,
       xlab = "Iteration",
       ylab = "ELBO")

  par(ask = TRUE)

  bar_plot_y <- round(fit[["model_parameters"]][["ppi"]][round(fit[["model_parameters"]][["ppi"]], digits = 4) > 0] * 100, digits = 2)

  bar_plot <-  barplot(bar_plot_y,
                       ylim = c(-10, 100),
                       yaxt = "n",
                       ylab = "Cluster Percentage (%)",
                       xlab = "Cluster"
  )
  abline(h = 0)
  axis(side = 2, at = seq(0, 100, 5),
       labels = paste(seq(0, 100, 5), "%", sep = ""))

  text(x = bar_plot, y = bar_plot_y + 5, labels = paste(bar_plot_y, "%", sep = ""))
  text(x = bar_plot, y = -5, labels = paste("Cluster", sort(unique(fit$cluster)), sep = " "))

  par(ask = FALSE)
}



























