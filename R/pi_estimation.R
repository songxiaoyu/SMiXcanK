#' Estimate cell-type proportions (SMiXcan wrapper)
#'
#' This is a thin wrapper that calls MiXcan's internal `pi_estimate`
#' (if available). You must have MiXcan installed.
#'
#' @param ... Arguments passed to MiXcan's internal function.
#' @return The value returned by MiXcan's `pi_estimate`.
#' @importFrom utils getFromNamespace
#' @export
#' @seealso \code{MiXcan} package.
pi_estimation <- function(...) {
  if (!requireNamespace("MiXcan", quietly = TRUE)) {
    stop("MiXcan is required but not installed.", call. = FALSE)
  }
  fn <- getFromNamespace("pi_estimation", "MiXcan")  # internal access
  fn(...)
}

