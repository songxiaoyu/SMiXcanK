#' PRIMO-based cell-type pattern classification with FDR filtering
#'
#' This function implements the complete PRIMO-based classification pipeline
#' used in S-MiXcan analyses. Marginal and joint p-values are first adjusted
#' using the Benjamini--Hochberg procedure. Genes labeled as cell-type-specific
#' are evaluated using marginal FDR (OR rule across cell types) and subsequently
#' classified into cell-type association patterns using PRIMO posterior
#' probabilities. Genes labeled as non-specific are evaluated using the fdr
#' test only and are not passed to PRIMO.
#'
#' Importantly, PRIMO posterior probabilities are not modified. Maximum a
#' posteriori (MAP) classification is performed only among significant
#' cell-type-specific genes after conditioning on being non-null (i.e., excluding
#' the all-zero association pattern).
#'
#' The function supports an arbitrary number of cell types (\eqn{K}) and
#' preserves PRIMO's native association-pattern ordering using
#' \code{Primo::make_qmat()}.
#'
#' @param merged A data.frame containing marginal p-values, joint p-values,
#'   and a column indicating whether each gene is cell-type-specific or
#'   non-specific.
#' @param pvals_names Character vector of column names for marginal p-values,
#'   one per cell type (length \eqn{K}).
#' @param p_join_name Character. Column name for the joint (non-specific)
#'   p-value.
#' @param type_col Character. Column indicating specificity labels.
#' @param fdr_cutoff Numeric. FDR threshold for declaring significance
#'   (default 0.1).
#' @param specific_label Character. Label identifying cell-type-specific rows.
#' @param unspecific_label Character. Label identifying non-specific rows.
#' @param ... Additional arguments passed to \code{Primo::Primo_pval()}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{out}{Original data.frame augmented with BH-adjusted p-values,
#'     PRIMO posterior probabilities (\code{post_*}), and MAP pattern labels.}
#'   \item{alt_props_used}{Prior non-null proportion used in PRIMO.}
#'   \item{sig_gene_percentage}{Estimated fraction of genes with marginal
#'     association signal.}
#'   \item{sig_spec_idx}{Indices of significant cell-type-specific genes.}
#'   \item{sig_uns_idx}{Indices of significant non-specific genes.}
#'   \item{tab_specific_patterns}{Table of MAP pattern counts among significant
#'     cell-type-specific genes.}
#'   \item{n_trait}{Named vector giving counts of genes uniquely associated
#'     with each cell type.}
#'   \item{n_shared_specific}{Number of significant cell-type-specific genes
#'     associated with two or more cell types.}
#'   \item{n_shared_nonspecific}{Number of significant non-specific genes.}
#'   \item{n_shared_total}{Total number of shared genes.}
#'   \item{post_colnames}{Names of PRIMO posterior probability columns.}
#' }
#'
#' @details
#' The prior non-null proportion supplied to PRIMO is estimated as
#' \deqn{\Pr(\text{any marginal FDR} < \alpha) / K,}
#' where \eqn{K} is the number of cell types and \eqn{\alpha} is the FDR cutoff.
#'
#' For MAP classification, posterior probabilities are renormalized over all
#' non-null patterns (excluding the all-zero pattern) and the most likely
#' pattern is selected.
#'
#' @seealso \code{\link[Primo]{Primo_pval}}, \code{\link[Primo]{make_qmat}}
#'
#' @examples
#' \dontrun{
#' res <- primo_pipeline_wrap(
#'   merged,
#'   pvals_names = c("p_1_ct2", "p_2_ct2"),
#'   p_join_name = "p_join_ct2",
#'   type_col = "type_ct2"
#' )
#'
#' res$n_shared_total
#' table(res$out$MAP_pattern_nonnull)
#' }
#'
#' @importFrom stats p.adjust
#' @importFrom Primo Primo_pval make_qmat
#' @export
primo_pipeline_wrap <- function(
    merged,
    pvals_names,
    p_join_name,
    type_col = "type_ct2",
    fdr_cutoff = 0.1,
    specific_label = "CellTypeSpecific",
    unspecific_label = "NonSpecific",
    ...
) {
  stopifnot(is.data.frame(merged))
  if (!type_col %in% names(merged)) stop("Missing column: ", type_col)

  missing <- setdiff(c(pvals_names, p_join_name), names(merged))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  K <- length(pvals_names)

  ## 1) BH FDR
  fdr_marg_names <- paste0("fdr_", pvals_names)
  for (j in seq_along(pvals_names)) {
    merged[[fdr_marg_names[j]]] <- p.adjust(merged[[pvals_names[j]]], method = "BH")
  }
  fdr_join_name <- paste0("fdr_", p_join_name)
  merged[[fdr_join_name]] <- p.adjust(merged[[p_join_name]], method = "BH")

  ## 2) alt_props
  sig_gene_percentage <- mean(apply(
    merged[, fdr_marg_names, drop = FALSE],
    1,
    function(r) any(r < fdr_cutoff)
  ))
  alt_props <- sig_gene_percentage / K

  ## 3) significance rules
  sig_spec_idx <- which(
    merged[[type_col]] == specific_label &
      apply(merged[, fdr_marg_names, drop = FALSE], 1, function(r) any(r < fdr_cutoff))
  )

  sig_uns_idx <- which(
    merged[[type_col]] == unspecific_label &
      merged[[fdr_join_name]] < fdr_cutoff
  )

  ## 4) PRIMO patterns
  Q <- suppressWarnings(Primo::make_qmat(1:K))
  patterns <- apply(Q, 1, paste0, collapse = "")
  post_colnames <- paste0("post_", patterns)

  post_mat <- as.data.frame(matrix(NA_real_, nrow = nrow(merged), ncol = length(post_colnames)))
  colnames(post_mat) <- post_colnames

  idx_spec <- which(merged[[type_col]] == specific_label)
  if (length(idx_spec) > 0) {
    res <- Primo::Primo_pval(
      pvals = as.matrix(merged[idx_spec, pvals_names, drop = FALSE]),
      alt_props = rep(alt_props, K),
      ...
    )
    pp <- as.data.frame(res$post_prob)
    colnames(pp) <- post_colnames
    post_mat[idx_spec, ] <- pp
  }

  ## 5) MAP excluding 00...0
  post00 <- paste0("post_", paste(rep("0", K), collapse = ""))
  non00_cols <- setdiff(post_colnames, post00)

  MAP_pattern_nonnull <- rep(NA_character_, nrow(merged))
  if (length(sig_spec_idx) > 0) {
    post_sig <- as.matrix(post_mat[sig_spec_idx, non00_cols, drop = FALSE])
    den <- rowSums(post_sig, na.rm = TRUE)
    ok <- den > 0
    post_sig[ok, ] <- post_sig[ok, , drop = FALSE] / den[ok]
    cls <- max.col(post_sig[ok, , drop = FALSE])
    MAP_pattern_nonnull[sig_spec_idx][ok] <- sub("^post_", "", non00_cols[cls])
  }

  ## 6) summary counts
  count_ones <- function(x) sum(strsplit(x, "")[[1]] == "1")
  n_shared_specific <- sum(vapply(MAP_pattern_nonnull[sig_spec_idx], count_ones, integer(1)) >= 2, na.rm = TRUE)

  n_trait <- sapply(seq_len(K), function(k) {
    pat <- paste0(rep("0", k - 1), "1", rep("0", K - k))
    sum(MAP_pattern_nonnull == pat, na.rm = TRUE)
  })
  names(n_trait) <- paste0("cell", seq_len(K))

  out <- cbind(merged, post_mat)
  out$MAP_pattern_nonnull <- MAP_pattern_nonnull

  list(
    out = out,
    alt_props_used = alt_props,
    sig_gene_percentage = sig_gene_percentage,
    sig_spec_idx = sig_spec_idx,
    sig_uns_idx = sig_uns_idx,
    tab_specific_patterns = table(MAP_pattern_nonnull[sig_spec_idx], useNA = "no"),
    n_trait = n_trait,
    n_shared_specific = n_shared_specific,
    n_shared_nonspecific = length(sig_uns_idx),
    n_shared_total = n_shared_specific + length(sig_uns_idx),
    post_colnames = post_colnames
  )
}

