#' PRIMO MAP Classification
#'
#' Run PRIMO separately on cell-type-specific and non-specific rows,
#' then combine posterior probabilities and MAP classifications
#' back into the original data frame.
#'
#' @param merged A data.frame containing p-values and a type column.
#' @param pvals_names Character vector of column names used as PRIMO inputs.
#' @param type_col Character. Column name indicating specificity
#'   (e.g. \code{"type"}, \code{"type_ct2"}).
#' @param alt_props Numeric scalar or numeric vector giving the prior
#'   proportion of non-null tests for each p-value column.
#'   If scalar, it is recycled to \code{length(pvals_names)}.
#' @param specific_label Character label identifying cell-type-specific rows.
#' @param unspecific_label Character label identifying non-specific rows.
#' @param ... Additional arguments passed to \code{Primo::Primo_pval()}
#'
#' @return A list with:
#' \describe{
#'   \item{out}{Original data.frame with posterior probabilities
#'   (prefixed by \code{"post__"}) and a \code{MAP} column appended.}
#'   \item{primo_specific}{PRIMO result object for specific rows (or NULL).}
#'   \item{primo_unspecific}{PRIMO result object for non-specific rows (or NULL).}
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Splits rows by \code{type_col}
#'   \item Runs PRIMO separately on specific and non-specific rows
#'   \item Combines posterior probabilities in original row order
#'   \item Assigns MAP class via maximum posterior probability
#' }
#'
#' The \code{alt_props} argument represents the assumed prior proportion
#' of non-null tests per input p-value (see PRIMO documentation).
#'
#' @seealso \code{\link[Primo]{Primo_pval}}
#'
#' @examples
#' \dontrun{
#' res <- primo_map_classify(
#'   merged,
#'   pvals_names = c("p_1_ct2", "p_2_ct2", "p_3"),
#'   type_col = "type_ct2"
#' )
#' head(res$out)
#' }
#'
#' @importFrom Primo Primo_pval
#' @export
primo_map_classify <- function(
    merged,
    pvals_names,
    type_col = "type",
    alt_props = 0.05,
    specific_label = "CellTypeSpecific",
    unspecific_label = "UnSpecific",
    ...
) {
  stopifnot(is.data.frame(merged))
  if (!type_col %in% names(merged)) {
    stop("Missing column: ", type_col)
  }

  missing <- setdiff(pvals_names, names(merged))
  if (length(missing) > 0) {
    stop("Missing p-value columns: ", paste(missing, collapse = ", "))
  }

  idx_spec <- which(merged[[type_col]] == specific_label)
  idx_uns  <- which(merged[[type_col]] == unspecific_label)

  # normalize alt_props
  if (length(alt_props) == 1) {
    alt_props <- rep(alt_props, length(pvals_names))
  }
  if (length(alt_props) != length(pvals_names)) {
    stop("alt_props must be length 1 or length(pvals_names)")
  }

  post_cols <- NULL
  res1 <- res2 <- NULL

  # --- specific ---
  if (length(idx_spec) > 0) {
    pmat1 <- as.matrix(merged[idx_spec, pvals_names, drop = FALSE])
    res1 <- Primo::Primo_pval(pvals = pmat1, alt_props = alt_props, ...)
    post_cols <- colnames(res1$post_prob)
  }

  # --- unspecific ---
  if (length(idx_uns) > 0) {
    pmat2 <- as.matrix(merged[idx_uns, pvals_names, drop = FALSE])
    res2 <- Primo::Primo_pval(pvals = pmat2, alt_props = alt_props, ...)
    if (is.null(post_cols)) {
      post_cols <- colnames(res2$post_prob)
    }
  }

  # allocate posterior matrix
  post_mat <- matrix(NA_real_, nrow = nrow(merged), ncol = length(post_cols))
  colnames(post_mat) <- paste0("post__", post_cols)

  if (!is.null(res1)) post_mat[idx_spec, ] <- as.matrix(res1$post_prob)
  if (!is.null(res2)) post_mat[idx_uns, ]  <- as.matrix(res2$post_prob)

  # MAP classification
  MAP <- apply(post_mat, 1, function(x) {
    if (all(is.na(x))) NA_character_
    else colnames(post_mat)[which.max(x)]
  })

  out <- cbind(merged, post_mat)
  out$MAP <- MAP

  list(
    out = out,
    primo_specific = res1,
    primo_unspecific = res2
  )
}
