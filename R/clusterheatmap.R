#' Generates a heat map based on an iClusterVB object
#'
#' @param fit A fitted iClusterVB object.
#' @param nvars The number of variables to show; defaults to all.
#' @param title Title of the heat map.
#' @param cols A vector of colors to use for the clusters; defaults to a random selection of colors.
#' @return Returns a heat map.
#' @examples
#' clusterheatmap(fit_iClusterVB, nvars = 50, title = "View 1", cols = c("green","blue","purple", "red"))


clusterheatmap <- function(fit, nvars = NULL, title = NULL, cols = NULL) {
  if(!is.null(nvars)) {
    nvars
  } else{
    nvars <- sum(sapply(fit$mydata, ncol))
  }
  if(!is.null(cols)) {
    cols
  } else{
    cols <- colors()[sample(1:600, size = fit$K)]
  }
  if(!is.null(title)) {
    title
  } else{
    title <- paste("Method = ",fit_iClusterVB$initial_values[[1]],
                   ", K = ", fit_iClusterVB$K, sep = "")
  }

  if(nvars > sum(sapply(fit$mydata, ncol))){warning("Warning: nvars exceeds the number of variables in the dataset")}

  df <- as.data.frame(t(matrix(as.numeric(unlist(fit$mydata)), nrow = nrow(fit$mydata[[1]]))[,1:nvars]) )
  mat_col <- data.frame(Clusters = as.numeric(fit$cluster))
  rownames(mat_col) <- colnames(df)
  df <- df[, order(as.numeric(fit$cluster))]
  mat_colors <- list(Clusters = cols)
  pheatmap::pheatmap(df,
                     main = as.character(title),
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                     treeheight_row=0,
                     treeheight_col=0,
                     scale = "row",
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE,
                     annotation_col   = mat_col,
                     annotation_colors= mat_colors)
}
