utils::globalVariables(c("..count..", "color", "log_pval", "param pval",
                         "param", "pval"))

#' Plot p-values of methytmle objects
#'
#' @param x Object of class \code{methytmle} as produced by an appropriate call
#'        to \code{methyvim}.
#' @param ... Additional arguments passed \code{plot} as necessary.
#' @param type The type of plot to build: one of side-by-side histograms (type
#'        "both") comparing raw p-values to FDR-adjusted p-values (using the
#'        FDR-MSA correction) or either of these two histogram separately. Set
#'        this argument to "raw_pvals" for a histogram of the raw p-values, and
#'        to "fdr_pvals" for a histogram of the FDR-corrected p-values.
#'
#' @return Object of class \code{ggplot} containing a histogram or side-by-side
#'         histograms of the raw (marginal) and corrected p-values, with the
#'         latter computed automatically using the method of Tuglus and van der
#'         Laan.
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select slice arrange transmute
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram xlab ylab ggtitle
#'             scale_colour_manual scale_fill_gradientn guides guide_legend
#'             theme_minimal
#' @importFrom wesanderson wes_palette
#' @importFrom gridExtra grid.arrange
#' @importFrom superheat superheat
#'
#' @export
#'
#' @method plot methytmle
#'
#' @examples
#' suppressMessages(library(SummarizedExperiment))
#' library(methyvimData)
#' data(grsExample)
#' var_int <- as.numeric(colData(grsExample)[, 1])
# TMLE procedure for the ATE parameter over M-values with Limma filtering
#' methyvim_out_ate <- suppressWarnings(
#'  methyvim(data_grs = grsExample, sites_comp = 25, var_int = var_int,
#'           vim = "ate", type = "Mval", filter = "limma", filter_cutoff = 0.1,
#'           parallel = FALSE, tmle_type = "glm"
#'          )
#' )
#' plot(methyvim_out_ate)
#
plot.methytmle <- function(x, ..., type = "both") {
  # use an older color palette to ensure compatibility with CRAN version
  pal <- wesanderson::wes_palette("Royal2", 100, type = "continuous")

  # get corrected p-values and add them to output object
  pval_fdr <- fdr_msa(pvals = x@vim$pval, total_obs = nrow(x))
  vim_table <- as.data.frame(cbind(x@vim, pval_fdr))

  # plot of raw p-values
  p1 <- ggplot2::ggplot(vim_table, ggplot2::aes(pval, xmin = 0, xmax = 1))
  p1 <- p1 + ggplot2::geom_histogram(ggplot2::aes(y = ..count..,
                                                  fill = ..count..),
                                     colour = "white", na.rm = TRUE,
                                     binwidth = 0.025)
  p1 <- p1 + ggplot2::ggtitle("Histogram of raw p-values")
  p1 <- p1 + ggplot2::xlab("magnitude of raw p-values")
  p1 <- p1 + ggplot2::scale_fill_gradientn("Count", colors = pal)
  p1 <- p1 + ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))
  p1 <- p1 + ggplot2::theme_minimal()

  # plot of FDR-corrected p-values
  p2 <- ggplot2::ggplot(vim_table, ggplot2::aes(pval_fdr, xmin = 0, xmax = 1))
  p2 <- p2 + ggplot2::geom_histogram(ggplot2::aes(y = ..count..,
                                                  fill = ..count..),
                                     colour = "white", na.rm = TRUE,
                                     binwidth = 0.025)
  p2 <- p2 + ggplot2::ggtitle("Histogram of FDR-corrected p-values")
  p2 <- p2 + ggplot2::xlab("magnitude of FDR-corrected p-values")
  p2 <- p2 + ggplot2::scale_fill_gradientn("Count", colors = pal)
  p2 <- p2 + ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))
  p2 <- p2 + ggplot2::theme_minimal()

  if (type == "both") {
    # return a grob with the two plots side-by-side
    gridExtra::grid.arrange(p1, p2, nrow = 1)
  } else if (type == "raw_pvals") {
    return(p1)
  } else if (type == "fdr_pvals") {
    return(p2)
  }
}

################################################################################

#' Heatmap for methytmle objects
#'
#' @param x Object of class \code{methytmle} as produced by an appropriate call
#'        to \code{methyvim}.
#' @param ... Additional arguments passed to \code{superheat}. Consult the
#'        documentation of the \code{superheat} package for a list of options.
#' @param n_sites Numeric indicating the number of CpG sites to be shown in the
#'        plot. If the number of sites analyzed is greater than this cutoff,
#'        sites to be displayed are chosen by ranking sites based on their raw
#'        (marginal) p-values.
#' @param type Whether to plot the original data (M-values or Beta-values) for
#'        the set of top CpG sites or to plot the measurements after applying a
#'        transformation into influence curve space (with respect to the target
#'        parameter of interest). The latter uses the fact that the parameters
#'        have asymptotically linear representations to obtain a rotation of the
#'        raw data into an alternative space; moreover, in this setting, the
#'        heatmap reduces to visualizing a supervised clustering procedure.
#'
#' @return Nothing. This function is called for its side-effect of outputting a
#'         heatmap to the graphics device. The heatmap is constructed using the
#'         \code{superheat} package.
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select slice arrange transmute
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram xlab ylab ggtitle
#'             scale_colour_manual scale_fill_gradientn guides guide_legend
#'             theme_minimal
#' @importFrom wesanderson wes_palette
#' @importFrom superheat superheat
#'
#' @export
#'
#' @examples
#' suppressMessages(library(SummarizedExperiment))
#' library(methyvimData)
#' data(grsExample)
#' var_int <- as.numeric(colData(grsExample)[, 1])
# TMLE procedure for the ATE parameter over M-values with Limma filtering
#' methyvim_out_ate <- suppressWarnings(
#'  methyvim(data_grs = grsExample, sites_comp = 25, var_int = var_int,
#'           vim = "ate", type = "Mval", filter = "limma", filter_cutoff = 0.1,
#'           parallel = FALSE, tmle_type = "glm"
#'          )
#' )
#' methyheat(methyvim_out_ate, type = "raw")
#
methyheat <- function(x, ..., n_sites = 25, type = "raw") {
  # elementary type checking
  #type <- check.args(type)

  # need observations in influence curve space to plot on heatmap
  if(type == "ic" & sum(dim(x@ic)) == 0) {
    stop("Please re-run 'methyvim' and set argument 'return_ic' to 'TRUE'.")
  }

  # set up annotations
  tx_annot <- ifelse(x@var_int == 0, "Control", "Treated")

  # rank sites based on raw p-value
  sites_ranked <- x@vim %>%
    data.frame() %>%
    dplyr::select(pval) %>%
    unlist() %>%
    as.numeric() %>%
    order()
  sites_mask <- x@screen_ind[sites_ranked]

  # subset matrix of measures/estimates to those for the top ranked sites
  if (type == "ic") {
    sites_mat <- x@ic %>%
      data.frame() %>%
      dplyr::slice(sites_mask) %>%
      as.matrix()
  } else if (type == "raw") {
    sites_mat <- SummarizedExperiment::assay(x) %>%
      data.frame() %>%
      dplyr::slice(sites_mask) %>%
      as.matrix()
  }

  # limit to the specified maximum number of sites
  if (!is.null(n_sites) & (nrow(sites_mat) > n_sites)) {
    sites_mat <- sites_mat %>%
      data.frame() %>%
      dplyr::slice(seq_len(n_sites)) %>%
      as.matrix()
  }

  # plot the (super) heatmap
  superheat::superheat(sites_mat, row.dendrogram = TRUE,
                       grid.hline.col = "white", force.grid.hline = TRUE,
                       grid.vline.col = "white", force.grid.vline = TRUE,
                       membership.cols = tx_annot,
                       title = paste("Heatmap of Top", nrow(sites_mat), "CpGs"),
                       ...
                      )
}

################################################################################

#' Volcano plot for methytmle objects
#'
#' @param x Object of class \code{methytmle} as produced by an appropriate call
#'        to \code{methyvim}.
#' @param param_bound Numeric for a threshold indicating the magnitude of the
#'        size of the effect considered to be interesting. This is used to
#'        assign groupings and colors to individual CpG sites.
#' @param pval_bound Numeric for a threshold indicating the magnitude of
#'        p-values deemed to be interesting. This is used to assign groupings
#'        and colors to individual CpG sites.
#'
#' @return Object of class \code{ggplot} containing a volcano plot of the
#'         estimated effect size on the x-axis and the -log10(p-value) on the
#'         y-axis. The volcano plot is used to detect possibly false positive
#'         cases, where a test statistic is significant due to low variance.
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select slice arrange transmute
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram xlab ylab ggtitle
#'             scale_colour_manual scale_fill_gradientn guides guide_legend
#'             theme_minimal xlim
#' @importFrom wesanderson wes_palette
#' @importFrom superheat superheat
#'
#' @export
#'
#' @examples
#' suppressMessages(library(SummarizedExperiment))
#' library(methyvimData)
#' data(grsExample)
#' var_int <- as.numeric(colData(grsExample)[, 1])
# TMLE procedure for the ATE parameter over M-values with Limma filtering
#' methyvim_out_ate <- suppressWarnings(
#'  methyvim(data_grs = grsExample, sites_comp = 25, var_int = var_int,
#'           vim = "ate", type = "Mval", filter = "limma", filter_cutoff = 0.1,
#'           parallel = FALSE, tmle_type = "glm"
#'          )
#' )
#' methyvolc(methyvim_out_ate)
#
methyvolc <- function(x, param_bound = 2.0, pval_bound = 0.2) {
  # use an older color palette to ensure compatibility with CRAN version
  pal <- wesanderson::wes_palette(name = "Darjeeling2", type = "continuous")

  # get corrected p-values
  pval_fdr <- fdr_msa(pvals = x@vim$pval, total_obs = nrow(x))
  vim_table <- as.data.frame(cbind(x@vim, pval_fdr))

  # set up object for plotting
  into_volcano <- vim_table %>%
    data.frame() %>%
    dplyr::arrange(pval) %>%
    dplyr::transmute(
      param = if(x@param == "Average Treatment Effect") {
          as.numeric(x@vim$est_ATE)
        } else if (x@param == "Relative Risk") {
          as.numeric(x@vim$est_logRR)
        },
      log_pval = -log10(pval),
      pval_fdr = I(pval_fdr),
      color = ifelse((param > param_bound) & (pval_fdr < pval_bound), "1",
                      ifelse((param < -param_bound) & (pval_fdr < pval_bound),
                              "-1", "0")
                    )
    )

  # create and return plot object
  p <- ggplot2::ggplot(into_volcano, ggplot2::aes(x = param, y = log_pval))
  p <- p + ggplot2::geom_point(ggplot2::aes(colour = color))
  p <- p + ggplot2::xlab(ifelse(x@param == "Relative Risk",
                                paste("Estimated log-Change in", x@param),
                                paste("Estimated Change in", x@param)
                               )
                         )
  p <- p + ggplot2::ylab("-log10(raw p-value)")
  p <- p + ggplot2::xlim(max(abs(into_volcano$param)) * c(-1, 1))
  p <- p + ggplot2::scale_colour_manual(values = pal[seq_len(3)],
                                        guide = FALSE)
  p <- p + ggplot2::theme_minimal()
  return(p)
}

