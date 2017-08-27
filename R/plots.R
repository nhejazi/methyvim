#' Common plots for methytmle objects
#'
#' Several types of plots providing...
#'
#'
#' @param x object of class \code{biotmle} as produced by an appropriate call to
#'        \code{biomarkertmle}
#' @param type character describing whether to provide a plot of unadjusted or
#'        adjusted p-values (adjustment performed via Benjamini-Hochberg)
#' @param ... additional arguments passed \code{plot} as necessary
#'
#'
#' @return object of class \code{ggplot} containing one of several types of
#'         plots: a heatmap of the top 50 CpG sites, a volcano plot displaying
#'         the parameter estimate against the log-transformed raw p-values, or
#'         side-by-side histograms of the raw and corrected p-values.
#'
#' @export
#'
#' @method plot methytmle
#'

plot.methytmle <- function(x, ..., type = c("heat", "volcano", "pvals")) {
  message("Plot methods have not yet been implemented.")
}

################################################################################

#' a heatmap for methytmle objects with more features goes here...

################################################################################

#' fancy annotation plot using GViz goes here...
