#' Parallelization with Futures and BiocParallel
#'
#' Easily set up a suitable parallelization scheme using the various options
#' provided in \code{BiocParallel} and packages of the \code{future} ecosystem.
#' INTERNAL USE ONLY.
#'
#' @param parallel Logical indicating whether parallelization ought to be used.
#'        Parallelization is invoked via a combination of \code{BiocParallel}
#'        and \code{future}. If \code{TRUE} the default method uses multiprocess
#'        evaluation, though other \code{future::plan}s may be specified using
#'        an optional argument. If \code{FALSE}, sequential computation is used
#'        and a warning message is issued.
#' @param future_param Character (if not \code{NULL}) specifying a particular
#'        parallelization approach to be used. For a list of the options, see
#'        the documentation for \code{future::plan}. If the previous argument
#'        (\code{parallel}) is set to \code{FALSE}, this argument is ignored and
#'        sequential computation is invoked via \code{future::sequential}.
#' @param bppar_type Character specifying the type of backend to be used with
#'        the parallelization invoked by \code{BiocParallel}. Consult the manual
#'        page for \code{BiocParallel::BiocParallelParam} for possible types and
#'        descriptions on their appropriate uses. The default for this argument
#'        is \code{NULL}, which silently uses \code{BiocParallel::DoparParam}.
#'
#' @return Nothing. This function is designed to be called for its side-effect
#'         of registering a parallel backend (for \code{BiocParallel}) and/or
#'         \code{future::plan}, making parallel computation a trivial process.
#'
#' @keywords internal
#'
#' @importFrom BiocParallel register bpprogressbar DoparParam
#' @importFrom future plan multiprocess sequential
#' @importFrom doFuture registerDoFuture
#
set_parallel <- function(parallel = c(TRUE, FALSE),
                         future_param = NULL,
                         bppar_type = NULL) {
  # invoke a future-based backend
  doFuture::registerDoFuture()

  if (parallel == TRUE) {
    if (!is.null(future_param)) {
      set_future_param <- parse(text = paste0("future", "::", future_param))
      future::plan(eval(set_future_param))
    } else {
      future::plan(future::multiprocess)
    }
  } else if (parallel == FALSE) {
    warning(paste(
      "Sequential evaluation is strongly discouraged.",
      "\n Proceed with caution."
    ))
    future::plan(future::sequential)
  }
  if (!is.null(bppar_type)) {
    bp_type <- eval(parse(text = paste0(
      "BiocParallel", "::",
      bppar_type, "()"
    )))
  } else {
    bp_type <- BiocParallel::DoparParam()
  }
  # try to use a progress bar is supported in the parallelization plan
  BiocParallel::bpprogressbar(bp_type) <- TRUE
  # register the chosen parallelization plan
  BiocParallel::register(bp_type, default = TRUE)
}

################################################################################

#' Wrap a Function in a Try Statement
#'
#' Function factory that generates versions of functions wrapped in \code{try}.
#' Originally fround in and borrowed from package \code{origami}.
#'
#' @param fun A \code{function} to be wrapped in a \code{try} statement.
#' @param ... Additional arguments passed to the previous argument \code{fun}.
#
wrap_in_try <- function(fun, ...) {
  wrapped <- function(...)
    try({
      fun(...)
    }, silent = TRUE)
  return(wrapped)
}

