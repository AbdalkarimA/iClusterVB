#' Generic summary method for `iClusterVB` objects
#'
#' @param object A fitted iClusterVB object.
#' @param rho The minimum posterior inclusion probability of interest to count
#'   the number of variables that are >= \code{rho}. Default is 0.5. Only works for
#'   VS_method = 1.
#' @param ... Potential further arguments
#' @return Returns a summary list for an `agnes` object.
#' @examples
#' ## S3 method for class 'iClusterVB'
#' summary(fit_iClusterVB, rho = 0.5, ...)
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
#' @param fit A fitted iClusterVB object.
#'
#' @return Returns an ELBO plot, a barplot of cluster proportions, and heatmaps
#'   for each data view.
#'
#'
#' @export
#' @method plot iClusterVB
#' @examples
#' plot(fit_iClusterVB)


plot.iClusterVB <- function(fit) {

  plot(x = 1:fit$iter,
       y = fit$elbo[1:1:fit$iter],
       type = "o",
       lwd = 2,
       xlab = "Iteration",
       ylab = "ELBO")

  devAskNewPage(ask = TRUE)

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
  text(x = bar_plot, y = -5, labels = paste("Cluster", unique(fit$cluster), sep = " "))


  formals(pheatmap)[c(
    "cluster_rows", "cluster_cols", "color", "treeheight_row", "treeheight_col", "scale",
    "show_colnames", "show_rownames", "annotation_names_row", "annotation_names_col"
  )] <- list(
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50), treeheight_row = 0,
    treeheight_col = 0, scale = "row", show_colnames = FALSE,
    show_rownames = FALSE, annotation_names_row = FALSE,
    annotation_names_col = FALSE
  )



  rho <- 0

  title <- paste(
    "View", 1:length(fit$mydata), "-",
    tools::toTitleCase(fit$dist[1:length(fit$mydata)]))


  if (is.null(fit$model_parameters$rho)) {
    for (i in 1:length(fit$mydata)) {
      df <- as.data.frame(t(data.matrix(fit$mydata[[i]])))
      mat_col <- data.frame(Clusters = as.numeric(fit$cluster))
      rownames(mat_col) <- colnames(df)
      df <- df[, order(as.numeric(fit$cluster))]
      pheatmap(df,
               main = as.character(title[i]),
               annotation_col = mat_col
      )
    }
  } else if (!is.null(fit$model_parameters$rho)) {
    for (i in 1:length(fit$mydata)) {
      df <- as.data.frame(t(data.matrix(fit$mydata[[i]][, fit$model_parameters$rho[[i]] >= rho])))
      mat_col <- data.frame(Clusters = as.numeric(fit$cluster))
      rownames(mat_col) <- colnames(df)
      df <- df[, order(as.numeric(fit$cluster))]
      pheatmap(df,
               main = as.character(title[i]),
               annotation_col = mat_col
      )
    }
  }

  devAskNewPage(ask = FALSE)

}



























