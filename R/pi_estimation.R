#' Estimate cell-type fractions for K cell types using BayesDeBulk
#'
#' This function runs \code{BayesDeBulk} on bulk expression data and returns
#' estimated cell-type fractions (π) for K cell types, along with the
#' estimated cell-type–specific expression.
#'
#' @param exprB A matrix or data frame of bulk expression values (genes × samples
#'   or samples × genes, formatted as expected by \code{BayesDeBulk}).
#' @param markers A list of length K, where \code{markers[[c]]} is a character
#'   vector of marker genes for cell type \code{cell.type[c]}.
#' @param seed Random seed.
#' @param n.iter Total number of MCMC iterations passed to \code{BayesDeBulk}.
#' @param burn.in Number of burn-in iterations passed to \code{BayesDeBulk}.
#' @param ... Additional arguments passed to \code{BayesDeBulk}.
#'
#' @return A list with components:
#' \describe{
#'   \item{cell_fraction}{A data frame of estimated cell-type fractions with
#'     one row per sample and one column per cell type, plus a \code{SampleID}
#'     column from row names.}
#'   \item{cell_expression}{Estimated cell-type–specific expression as returned
#'     by \code{BayesDeBulk}.}
#' }
#'
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom BayesDeBulk BayesDeBulk
#' @export
pi_estimation_K <- function(exprB,
                            markers,
                            seed   = 1,
                            n.iter = 10000,
                            burn.in = 1000,
                            ...) {
  cell.type<-names(markers)
  index.matrix<-NULL
  for (s in 1:length(cell.type)){
    for (k in 1:length(cell.type)){
      if (s!=k){
        mg<-match(markers[[s]],markers[[k]])
        index.matrix<-rbind(index.matrix,cbind(rep(cell.type[s],sum(is.na(mg))),
                                               rep(cell.type[k],sum(is.na(mg))),markers[[s]][is.na(mg)]))
      }
    }
  }
  set.seed(seed)

  fit <- BayesDeBulk(
    n.iter  = n.iter,
    burn.in = burn.in,
    Y       = list(exprB),
    markers = index.matrix,
    ...
  )

  cell_fraction <- fit$cell.fraction %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID")

  list(
    cell_fraction  = cell_fraction,
    cell_expression = fit$cell.expression
  )
}
