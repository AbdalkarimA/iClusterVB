#' Generates a probability inclusion plot based on an iClusterVB object
#'
#' @param fit A fitted iClusterVB object.
#' @param combined_plot LOGICAL. Whether to plot all the plots together or separately; defaults to TRUE (together).
#' @param ylab The y-axis label; defaults to "Probability of Inclusion".
#' @param default_title LOGICAL; defaults is TRUE, which outputs "View 1,..R - Distribution".
#' @param title The title of the plots. It can be a character vector or a single value; default outputs the default title.
#' @return Returns a probability inclusion plot or plots.
#' @examples
#' piplot(fit_iClusterVB, ylab = "Probability" , default_title = TRUE, title = NULL)
#' piplot(fit_iClusterVB)

piplot <- function(fit, combined_plot = TRUE, ylab = "Probability of Inclusion", default_title = TRUE, title = NULL) {
  for (i in 1:length(fit$model_parameters$rho)) {
    assign(paste("dat", i, sep = ""), data.frame(varid = 1:length(fit$model_parameters$rho[[i]]),
                                                 rho = t(fit$model_parameters$rho[[i]])))
  }

  gp_rho <- list()

  if (default_title == TRUE & is.null(title)) {
    title <- c()
    title <- paste("View", 1:length(fit$model_parameters$rho), "-",
                   tools::toTitleCase(fit$dist[1:length(fit$model_parameters$rho)]))
  }

  for(i in 1:length(fit$model_parameters$rho)) {

    gp_rho[[i]] <-  ggplot2::ggplot(data = get(paste("dat", i, sep = "")), mapping = ggplot2::aes(x=reorder(varid, rho), y=rho, fill=rho )) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::ggtitle(ifelse(length(title) >1, title[i], title)) +
      ggplot2::ylim (c(0,1))+
      ggplot2::geom_hline( yintercept = 0.5, color = "red",  linewidth = 1, linetype = 2) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 15, face = "bold"),
                     axis.text= ggplot2::element_text(size=1),
                     axis.title= ggplot2::element_text(size=15),
                     axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5, hjust=1),
                     strip.text.x = ggplot2::element_text(size = 1, angle = 0),
                     strip.text.y = ggplot2::element_text(size = 15,face="bold")) +
      ggplot2::xlab(" ") + ggplot2::ylab(if(ylab != "Probability of Inclusion") {ylab}
                                         else{ylab})  +
      ggplot2::scale_fill_gradientn( colors = topo.colors(2), limits=c(0,1)) +
      ggplot2::coord_polar()
  }

  if(combined_plot == TRUE) {
    cowplot::plot_grid(plotlist = gp_rho,
                       labels= paste("(", LETTERS[1:length(fit$model_parameters$rho)], ")", sep = ""),
                       nrow = 2,
                       align = "hv")
  } else {
    gp_rho[1:length(fit$model_parameters$rho)]
  }
}
