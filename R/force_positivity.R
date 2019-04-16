#' Enforce the Assumption of Positivity
#'
#' Discretize continuous variables in the adjustment set (W) of a TMLE procedure
#' in order to avoid practical violations of the assumption of positivity.
#' Discretizes the numeric columns of an input matrix such that the newly
#' created levels of each variable individually contain at least a specified
#' mass when considering each level against levels of the treatment variable.
#' INTERNAL USE ONLY.
#'
#' @param A Numeric giving the levels of the (discretized) treatment variable.
#' @param W Data.Frame or Matrix containing the covariates in the adjustment set
#'  to be discretized against the levels of the treatment variable.
#' @param pos_min Numeric indicating the minimum mass (as a proportion) of the
#'  observations to be included in any cell of the table composed of the levels
#'  of the treatment against levels of an adjustment covariate.
#' @param q_init Numeric indicating the initial number of levels to discretize a
#'  given adjustment variable into. This defaults to quantiles.
#'
#' @return A numeric vector with the adjustment variables re-coded into discrete
#'  levels respecting the minimum mass requested in each table comparing levels
#'  of the treatment against levels of an adjustment covariate.
#'
#' @keywords internal
#'
#' @importFrom gtools quantcut
force_positivity <- function(A, W, pos_min = 0.1, q_init = 10) {
  stopifnot(length(A) == nrow(W))

  if (class(W) != "data.frame") W <- as.data.frame(W) # cover use of "ncol"
  out_w <- NULL # concatenate W columnwise as we discretize each covar below

  for (obs_w in seq_len(ncol(W))) {
    in_w <- as.numeric(W[, obs_w])
    discr_w <- as.numeric(as.factor(gtools::quantcut(x = in_w, q = q_init)))
    check <- sum((table(A, discr_w) / length(A)) < pos_min)
    next_guess_q <- q_init
    while (check > 0) {
      next_guess_q <- (next_guess_q - 1)
      discr_w <- as.numeric(as.factor(gtools::quantcut(
        x = in_w,
        q = next_guess_q
      )))
      check <- sum((table(A, discr_w) / length(A)) < pos_min)
    }
    out_w <- cbind(out_w, discr_w)
  }
  out <- as.data.frame(out_w)
  colnames(out) <- colnames(W)
  rownames(out) <- rownames(W)
  if (length(which(colSums(out) == length(A))) > 0) {
    out <- out[, -which(colSums(out) == length(A)), drop = FALSE]
  }
  return(out)
}
