#' Generates a heat map based on an iClusterVB object
#'
#' @param fit A fitted iClusterVB object.
#' @param rho The minimum probability of inclusion for variables shown on the
#'   heatmap. Default is 0. Only useful for VS_method = 1.
#' @param title A character vector or a single value. Title of the heat map. The
#'   default is "View 1 - Distribution 1", ..., "View R - Distribution R".
#' @param cols A vector of colors to use for the clusters. The default is a random selection of colors.
#' @param ... Additional arguments to be passed down to
#'   \code{\link[pheatmap]{pheatmap}}
#'
#' @return Returns a heat map for each data view.
#' @examples
#' chmap(fit_iClusterVB, rho = 0.75, title = c("View 1", "View 2"), cols = c("green", "blue", "purple", "red"))
#' chmap(fit_iClusterVB)
#'
#' @export
#' @import pheatmap
#' @useDynLib iClusterVB, .registration=TRUE


chmap <- function(fit, rho = 0, cols = NULL, title = NULL, ...) {

  if (is.null(cols)) {
    cols <- colors()[sample(1:600, size = length(unique(fit$cluster)))]
  }


  formals(pheatmap)[c(
    "cluster_rows", "cluster_cols", "color", "treeheight_row", "treeheight_col", "scale",
    "show_colnames", "show_rownames", "annotation_names_row", "annotation_names_col", "annotation_colors"
  )] <- list(
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50), treeheight_row = 0,
    treeheight_col = 0, scale = "row", show_colnames = FALSE,
    show_rownames = FALSE, annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    annotation_colors = setNames(as.list(cols), paste("Cluster", sort(unique(fit$cluster))))
  )



  ifelse(is.null(title), title <- paste(
    "View", 1:length(fit$mydata), "-",
    tools::toTitleCase(fit$dist[1:length(fit$mydata)])
  ),
  title
  )


  plot_list <- list()

  if (is.null(fit$model_parameters$rho)) {
    for (i in 1:length(fit$mydata)) {
      df <- as.data.frame(t(data.matrix(fit$mydata[[i]])))
      mat_col <- data.frame(Clusters = paste("Cluster", as.numeric(fit$cluster)))
      rownames(mat_col) <- colnames(df)
      df <- df[, order(as.numeric(fit$cluster))]
      mat_colors <- list(Clusters = cols)
      plot_list[[i]] <- pheatmap(df,
                                 main = as.character(title[i]),
                                 annotation_col = mat_col,
                                 ...
      )[[4]]
    }
  } else if (!is.null(fit$model_parameters$rho)) {
    for (i in 1:length(fit$mydata)) {
      df <- as.data.frame(t(data.matrix(fit$mydata[[i]][, fit$model_parameters$rho[[i]] >= rho])))
      mat_col <- data.frame(Clusters = paste("Cluster", as.numeric(fit$cluster)))
      rownames(mat_col) <- colnames(df)
      df <- df[, order(as.numeric(fit$cluster))]
      mat_colors <- list(Clusters = cols)
      plot_list[[i]] <- pheatmap(df,
                                 main = as.character(title[i]),
                                 annotation_col = mat_col,
                                 ...
      )[[4]]
    }
  }

  return(invisible(plot_list))
}
