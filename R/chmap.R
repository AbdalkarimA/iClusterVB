#' Generates a heat map based on an iClusterVB object
#'
#' @param fit A fitted iClusterVB object.
#' @param nvars A numeric vector or a single value. The number of variables to show; defaults to all.
#' @param title A character vector or a single value. Title of the heat map; defaults to "View 1,..R - Distribution".
#' @param cols A vector of colors to use for the clusters; defaults to a random selection of colors.
#' @return Returns a heat map.
#' @examples
#' chmap(fit_iClusterVB, nvars = 50, title = c("View 1", "View 2"), cols = c("green","blue","purple", "red"))
#' chmap(fit_iClusterVB)

chmap <- function(fit, nvars = NULL, title = NULL, cols = NULL) {

  if(!is.null(nvars) & length(nvars) == 1) {
    nvars <- rep(nvars, length(fit$mydata))
  }

  if(!is.null(title) & length(title) == 1) {
    title <- rep(title, length(fit$mydata))
  }

  if(is.null(nvars)) {
    nvars <- unlist(fit$data_summary[[3]])
  }

  if(is.null(cols)) {
    cols <- colors()[sample(1:600, size = length(unique(fit$cluster)))]
  }

  ifelse(is.null(title), title <-paste("View", 1:length(fit$mydata), "-",
                                       tools::toTitleCase(fit$dist[1:length(fit$mydata)])),
         title)

  for (i in 1:length(fit$mydata)) {
    df <- as.data.frame(t(matrix(as.numeric(unlist(fit$mydata[[i]])), nrow = nrow(fit$mydata[[i]]))[,1:nvars[i]]) )
    mat_col <- data.frame(Clusters = as.numeric(fit$cluster))
    rownames(mat_col) <- colnames(df)
    df <- df[, order(as.numeric(fit$cluster))]
    mat_colors <- list(Clusters = cols)
    pheatmap::pheatmap(df,
                       main = as.character(title[i]),
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
}
