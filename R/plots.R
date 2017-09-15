utils::globalVariables(c("..count..", "color", "log_pval", "param pval",
                         "param", "pval"))

#' Plot utility for methytmle objects
#'
#' Several types of plots providing...
#'
#'
#' @param x Object of class \code{methytmle} as produced by an appropriate call
#'        to \code{methyvim}.
#' @param ... additional arguments passed \code{superheat}. This argument is
#'        used only if \code{type} below is set to the heatmap option. Please
#'        consult the documentation of the \code{superheat} package for options
#'        that may be appropriately set.
#' @param type Character describing the type of plot to be generated: one among
#'        the following collection: heatmap, volcano plot, or histograms of
#'        raw and FDR-corrected p-values.
#' @param num_sites Numeric indcating the number of CpG sites to be shown in the
#'        plot. This argument is only used if \code{type} is set to the heatmap
#'        option. If the number of sites analyzed is greater than this cutoff,
#'        sites to be displayed are chosen by ranking sites based on the raw
#'        (marginal) p-values.
#' @param param_bound Numeric for a threshold indicating the magnitude of the
#'        size of the effect considered to be interesting. This is used to
#'        assign groupings and colors to individual CpG sites. This argument is
#'        only used when \code{type} is set to the volcano plot option.
#' @param pval_bound Numeric for a threshold indicating the magnitude of
#'        p-values deemed to be interesting. This is used to assign groupings
#'        and colors to individual CpG sites. This argument is only used when
#'        \code{type} is set to the volcano plot option.
#'
#' @importFrom SummarizedExperiment colData
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

plot.methytmle <- function(x,
                           ...,
                           type = c("heatmap", "volcano", "pvals"),
                           num_sites = 25,
                           param_bound = 2.0,
                           pval_bound = 0.2) {
  if (type == "heatmap") {
    # need observations in influence curve space to plot on heatmap
    if(sum(dim(x@ic)) == 0) {
      stop("Please re-run 'methyvim' and set argument 'return_ic' to 'TRUE'.")
    }

    # use an older color palette to ensure compatibility with CRAN version
    pal <- wesanderson::wes_palette(name = "Moonrise3", n = 100,
                                    type = "continuous")

    # set up annotations
    tx_annot <- ifelse(SummarizedExperiment::colData(x)[, x@var_int] == 0,
                       "Control", "Treated")

    if (nrow(x@ic) > num_sites) {
      # rank sites based on raw p-value
      sites_ranked <- x@vim %>%
        data.frame() %>%
        dplyr::select(pval) %>%
        unlist() %>%
        as.numeric() %>%
        order()
      # subset matrix of IC estimates to those for the top ranked sites
      sites_mat <- x@ic %>%
        data.frame() %>%
        dplyr::slice(sites_ranked) %>%
        as.matrix()
    } else {
      sites_mat <- x@ic %>%
        as.matrix()
    }

    # plot the super heatmap
    superheat::superheat(sites_mat, row.dendrogram = TRUE,
                         grid.hline.col = "white", force.grid.hline = TRUE,
                         grid.vline.col = "white", force.grid.vline = TRUE,
                         membership.cols = tx_annot, heat.pal = pal,
                         title = paste("Heatmap of Top", num_sites, "CpGs"),
                         ...
                        )

  } else if (type == "volcano") {
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
    p <- p + ggplot2::scale_colour_manual(values = pal[seq_len(3)],
                                          guide = FALSE)
    p <- p + ggplot2::theme_minimal()
    return(p)

  } else {
    # use an older color palette to ensure compatibility with CRAN version
    pal <- wesanderson::wes_palette("Royal2", 100, type = "continuous")

    # get corrected p-values and add them to output object
    pval_fdr <- fdr_msa(pvals = x@vim$pval, total_obs = nrow(x))
    vim_table <- as.data.frame(cbind(x@vim, pval_fdr))

    # plot of raw p-values
    p1 <- ggplot2::ggplot(vim_table, ggplot2::aes(pval))
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
    p2 <- ggplot2::ggplot(vim_table, ggplot2::aes(pval_fdr))
    p2 <- p2 + ggplot2::geom_histogram(ggplot2::aes(y = ..count..,
                                                    fill = ..count..),
                                       colour = "white", na.rm = TRUE,
                                       binwidth = 0.025)
    p2 <- p2 + ggplot2::ggtitle("Histogram of BH-corrected FDR p-values")
    p2 <- p2 + ggplot2::xlab("magnitude of BH-corrected p-values")
    p2 <- p2 + ggplot2::scale_fill_gradientn("Count", colors = pal)
    p2 <- p2 + ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))
    p2 <- p2 + ggplot2::theme_minimal()

    # return a grob with the two plots side-by-side
    gridExtra::grid.arrange(p1, p2, nrow = 1)
  }
}
