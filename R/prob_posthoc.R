#' PRIMO MAP classification with FDR filtering
#'
#' @export
primo_map_classify <- function(
    pvals,
    merged_ctspec,
    fdr_cols = c("fdr_p_1_ct2", "fdr_p_2_ct2"),
    fdr_thresh = 0.1,
    alt_props = rep(0.05, ncol(pvals)),
    ...
) {
  stopifnot(
    nrow(pvals) == nrow(merged_ctspec),
    length(fdr_cols) == 2,
    all(fdr_cols %in% names(merged_ctspec)),
    length(alt_props) == ncol(pvals)
  )

  res <- Primo::Primo_pval(pvals = as.matrix(pvals),
                           alt_props = alt_props, ...)

  keep <- merged_ctspec[[fdr_cols[1]]] < fdr_thresh |
    merged_ctspec[[fdr_cols[2]]] < fdr_thresh

  post_prob_sel <- res$post_prob[keep, , drop = FALSE]
  MAP_class <- max.col(post_prob_sel)

  list(
    res = res,
    keep = keep,
    post_prob_selected = post_prob_sel,
    MAP_class = MAP_class
  )
}
