#' Generates a probability inclusion plot based on an iClusterVB object
#'
#' @param fit A fitted iClusterVB object.
#' @param plot_grid LOGICAL. Whether to use the \code{\link[cowplot]{plot_grid}}
#'   function from the \bold{cowplot} package. The default is TRUE.
#' @param ylab The y-axis label. The default is "Probability of Inclusion".
#' @param title The title of the plots. It can be a character vector or a single
#'   value. The default output is "View 1 - Distribution 1", ..., "View R -
#'   Distribution R".
#' @param ... Additional arguments to add to the
#'   \code{\link[cowplot]{plot_grid}} function.
#'
#' @return Returns a probability inclusion plot or plots.
#' @examples
#' \dontrun{
#' piplot(fit_iClusterVB,
#'   ylab = "Probability", title = NULL, align = "hv",
#'   nrow = 2, ncol = 2
#' )
#' piplot(fit_iClusterVB)
#' }
#' @import ggplot2 cowplot
#' @export
#' @useDynLib iClusterVB, .registration=TRUE


piplot <- function(fit,
                   plot_grid = TRUE,
                   ylab = "Probability of Inclusion",
                   title = NULL, ...) {
  for (i in 1:length(fit$model_parameters$rho)) {
    assign(paste("dat", i, sep = ""), data.frame(
      varid = 1:length(fit$model_parameters$rho[[i]]),
      rho = t(fit$model_parameters$rho[[i]])
    ))
  }

  gp_rho <- list()

  if (is.null(title)) {
    title <- c()
    title <- paste(
      "View", 1:length(fit$model_parameters$rho), "-",
      tools::toTitleCase(fit$dist[1:length(fit$model_parameters$rho)])
    )
  }

  for (i in 1:length(fit$model_parameters$rho)) {
    gp_rho[[i]] <- ggplot2::ggplot(
      data = get(paste("dat", i, sep = "")),
      mapping = ggplot2::aes(x = reorder(varid, rho), y = rho, fill = rho)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::ggtitle(ifelse(length(title) > 1, title[i], title)) +
      ggplot2::ylim(c(0, 1)) +
      ggplot2::geom_hline(yintercept = 0.5, color = "red", linewidth = 1, linetype = 2) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 15, face = "bold"),
        axis.text = ggplot2::element_text(size = 1),
        axis.title = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text.x = ggplot2::element_text(size = 1, angle = 0),
        strip.text.y = ggplot2::element_text(size = 15, face = "bold")
      ) +
      ggplot2::xlab(" ") +
      ggplot2::ylab(if (ylab != "Probability of Inclusion") {
        ylab
      } else {
        ylab
      }) +
      ggplot2::scale_fill_gradientn(colors = topo.colors(2), limits = c(0, 1)) +
      ggplot2::coord_polar()
  }


  args <- list(...)

  if (!("labels" %in% names(args))) {
    labels <- paste("(", LETTERS[1:length(fit$model_parameters$rho)], ")", sep = "")
  }

  if (plot_grid == TRUE) {
    cowplot::plot_grid(
      plotlist = gp_rho,
      labels = labels,
      ...
    )
  } else {
    gp_rho[1:length(fit$model_parameters$rho)]
  }
}

utils::globalVariables(c("rho", "varid"))
