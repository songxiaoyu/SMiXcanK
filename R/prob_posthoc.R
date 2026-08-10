#' Infer cell-type association patterns from SMiXcan results
#'
#' This function post-processes genome-wide SMiXcan results using
#' Benjamini-Hochberg FDR adjustment and PRIMO posterior probabilities.
#' It is designed for tables that contain:
#' marginal p-values for each cell type, a joint p-value, and a label
#' indicating whether a gene is cell-type-specific or non-specific.
#'
#' @param merged A data frame containing gene-level SMiXcan results.
#' @param pvals_names Character vector of marginal p-value column names
#'   for the K cell types.
#' @param p_join_name Character scalar giving the joint p-value column name.
#' @param type_col Character scalar giving the column that stores
#'   specificity labels. Default: \code{"type_ct2"}.
#' @param gene_id_col Character scalar giving the preferred gene identifier
#'   column name in the output. If this column is missing, the function
#'   falls back to \code{"gene_id"}, then \code{"gene_name"},
#'   then \code{"gene"}, and finally \code{"ENSG"}.
#' @param fdr_cutoff Numeric FDR threshold used to define significance.
#'   Default: \code{0.1}.
#' @param specific_label Character label used for cell-type-specific genes.
#' @param unspecific_label Character label used for non-specific genes.
#' @param ... Additional arguments passed to \code{Primo::Primo_pval()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{out}{Input table augmented with BH-adjusted p-values,
#'     PRIMO posterior probabilities, and MAP pattern assignments.}
#'   \item{gene_id_col_used}{Gene identifier column used in the output.}
#'   \item{fdr_join_name}{Name of the BH-adjusted joint p-value column.}
#'   \item{sig_join_idx}{Row indices of genes passing the joint FDR cutoff.}
#'   \item{sig_spec_idx}{Row indices of significant cell-type-specific genes.}
#'   \item{sig_uns_idx}{Row indices of significant non-specific genes.}
#'   \item{alt_props_used}{Prior non-null proportion supplied to PRIMO.}
#'   \item{sig_gene_percentage}{Observed fraction of genes passing the joint
#'     FDR cutoff.}
#'   \item{post_colnames}{Names of the posterior-probability columns added
#'     to \code{out}.}
#'   \item{tab_specific_patterns}{Table of MAP pattern counts among
#'     significant cell-type-specific genes.}
#'   \item{genes_sig_specific}{Gene IDs of significant cell-type-specific genes.}
#'   \item{genes_sig_nonspecific}{Gene IDs of significant non-specific genes.}
#'   \item{genes_by_pattern}{Named list of significant cell-type-specific genes
#'     grouped by inferred MAP pattern.}
#'   \item{unique_by_cell}{Named list of genes uniquely assigned to each cell
#'     type pattern.}
#'   \item{shared_specific_genes}{Gene IDs of significant cell-type-specific
#'     genes assigned to patterns involving two or more cell types.}
#' }
#'
#' @details
#' The joint p-value is adjusted across all genes using the BH procedure.
#' PRIMO is then applied only to genes labeled as cell-type-specific.
#' MAP pattern assignment is performed only among joint-significant
#' cell-type-specific genes, after excluding the all-null pattern.
#'
#' @importFrom stats p.adjust
#' @importFrom Primo Primo_pval make_qmat
#' @export
infer_celltype_patterns <- function(
    merged,
    pvals_names,
    p_join_name,
    type_col = "type_ct2",
    gene_id_col = "gene_id_clean",
    fdr_cutoff = 0.1,
    specific_label = "CellTypeSpecific",
    unspecific_label = "NonSpecific",
    alt_props=NULL,
    ...
) {
  stopifnot(is.data.frame(merged))
  if (!type_col %in% names(merged)) stop("Missing column: ", type_col)

  missing <- setdiff(c(pvals_names, p_join_name), names(merged))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  K <- length(pvals_names)
  if (K < 1) stop("pvals_names must have length >= 1")

  # Decide which column to use as the gene identifier in downstream outputs.
  if (!gene_id_col %in% names(merged)) {
    fallback <- c("gene_id", "gene_name", "gene", "ENSG")
    hit <- fallback[fallback %in% names(merged)]
    if (length(hit) == 0) {
      stop("gene_id_col not found, and no fallback gene ID columns found. ",
           "Please provide gene_id_col that exists in merged.")
    }
    gene_id_col <- hit[1]
  }

  ## 1) BH adjustment on the joint p-value
  fdr_join_name <- paste0("fdr_", p_join_name)
  merged[[fdr_join_name]] <- stats::p.adjust(merged[[p_join_name]], method = "BH")

  ## 2) Define significance using the joint FDR only
  sig_join_idx <- which(merged[[fdr_join_name]] < fdr_cutoff)

  # Among joint-significant genes, split by specificity label.
  sig_spec_idx <- intersect(sig_join_idx, which(merged[[type_col]] == specific_label))
  sig_uns_idx  <- intersect(sig_join_idx, which(merged[[type_col]] == unspecific_label))

  ## 3) Choose alt_props for Primo
  # Since marginal p-values are not thresholded directly here, use a simple
  # heuristic: the fraction of joint-significant genes, distributed across K
  # cell types.
  if (is.null(alt_props)) {
    sig_gene_percentage <- mean(merged[[fdr_join_name]] < fdr_cutoff, na.rm = TRUE)
    alt_props <- sig_gene_percentage / K
    alt_props <- min(max(alt_props, 1e-6), 1 - 1e-6)
  }
  ## 4) Fit PRIMO for all CellTypeSpecific genes
  Q <- suppressWarnings(Primo::make_qmat(1:K))
  patterns <- apply(Q, 1, paste0, collapse = "")
  post_colnames <- paste0("post_", patterns)

  post_mat <- as.data.frame(matrix(NA_real_, nrow = nrow(merged), ncol = length(post_colnames)))
  colnames(post_mat) <- post_colnames

  idx_spec_all <- which(merged[[type_col]] == specific_label)
  if (length(idx_spec_all) > 0) {
    res <- Primo::Primo_pval(
      pvals = as.matrix(merged[idx_spec_all, pvals_names, drop = FALSE]),
      alt_props = rep(alt_props, K)
    )
    pp <- as.data.frame(res$post_prob)
    colnames(pp) <- post_colnames
    post_mat[idx_spec_all, ] <- pp
  }

  ## 5) Assign MAP patterns excluding the all-null configuration
  post00 <- paste0("post_", paste(rep("0", K), collapse = ""))
  non00_cols <- setdiff(post_colnames, post00)

  MAP_pattern_nonnull <- rep(NA_character_, nrow(merged))

  if (length(sig_spec_idx) > 0) {
    post_sig <- as.matrix(post_mat[sig_spec_idx, non00_cols, drop = FALSE])
    den <- rowSums(post_sig, na.rm = TRUE)
    ok <- den > 0

    if (any(ok)) {
      # Normalize over non-null patterns only.
      post_sig[ok, ] <- sweep(post_sig[ok, , drop = FALSE], 1, den[ok], "/")
      cls <- max.col(post_sig[ok, , drop = FALSE])

      idx_ok <- sig_spec_idx[ok]
      MAP_pattern_nonnull[idx_ok] <- sub("^post_", "", non00_cols[cls])
    }
  }

  ## 6) Assemble the output table
  out <- cbind(merged, post_mat)
  out$MAP_pattern_nonnull <- MAP_pattern_nonnull

  ## 7) Create convenient summaries and gene lists
  count_ones <- function(x) {
    if (is.na(x)) return(NA_integer_)
    sum(strsplit(x, "")[[1]] == "1")
  }

  # Joint-significant CellTypeSpecific genes with a valid MAP assignment.
  spec_sig <- out[sig_spec_idx, , drop = FALSE]
  spec_sig <- spec_sig[!is.na(spec_sig$MAP_pattern_nonnull), , drop = FALSE]

  genes_sig_specific <- spec_sig[[gene_id_col]]
  genes_sig_nonspecific <- out[sig_uns_idx, gene_id_col, drop = TRUE]

  genes_by_pattern <- split(spec_sig[[gene_id_col]], spec_sig$MAP_pattern_nonnull)

  unique_by_cell <- lapply(seq_len(K), function(k) {
    pat <- paste0(rep("0", k - 1), "1", rep("0", K - k))
    spec_sig[[gene_id_col]][spec_sig$MAP_pattern_nonnull == pat]
  })
  names(unique_by_cell) <- paste0("cell", seq_len(K))

  shared_specific_genes <- spec_sig[[gene_id_col]][
    vapply(spec_sig$MAP_pattern_nonnull, count_ones, integer(1)) >= 2
  ]

  tab_specific_patterns <- table(spec_sig$MAP_pattern_nonnull, useNA = "no")

  list(
    out = out,
    gene_id_col_used = gene_id_col,

    # joint significance
    fdr_join_name = fdr_join_name,
    sig_join_idx = sig_join_idx,
    sig_spec_idx = sig_spec_idx,
    sig_uns_idx = sig_uns_idx,

    # PRIMO bookkeeping
    alt_props_used = alt_props,
    sig_gene_percentage = sig_gene_percentage,
    post_colnames = post_colnames,

    # pattern summaries (among joint-significant specific genes)
    tab_specific_patterns = tab_specific_patterns,

    # gene lists
    genes_sig_specific = genes_sig_specific,
    genes_sig_nonspecific = genes_sig_nonspecific,
    genes_by_pattern = genes_by_pattern,
    unique_by_cell = unique_by_cell,
    shared_specific_genes = shared_specific_genes
  )
}
