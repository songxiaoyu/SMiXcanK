#' Train Cell-Type-Specific Prediction Model
#'
#' Trains a cell-type-aware prediction model using MiXcan, estimating weights for SNPs across cell types.
#'
#' @name train_prediction_model
#' @title Train Cell-Type-Specific Prediction Model Using MiXcan
#'
#' @param y.train Gene expression levels in training data.
#' @param x.train Genotype matrix in training data.
#' @param pi.train Estimated cell-type proportions in training data.
#' @param cov A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param xNameMatrix Default is NULL. A matrix to save theX matrix information,
#' such as variable ID, position, rsid, ref_allele, eff_allele.
#' @param yName Default is NULL. A row vector to save the expression information, such as gene ID, gene name.
#' @param foldid Default is NULL. 10-fold cross-validation (CV) is used in our pipeline. A random split
#'
#' @return A list containing:
#' \describe{
#'   \item{W1}{Weight vector for cell type 1}
#'   \item{W2}{Weight vector for cell type 2}
#'   \item{selected_snp}{Indices of selected SNPs}
#'   \item{MiXcan_weight_result}{Total result matrix}
#' }
#'
#' @importFrom MiXcan MiXcan MiXcan_extract_weight
#' @export

train_prediction_model <- function(y.train, x.train, pi.train, cov=NULL, xNameMatrix = NULL, yName = NULL, foldid=NULL){
  n_train = nrow(x.train)
  MiXcan_result <- MiXcan(y= y.train, x = x.train,
                          pi = pi.train,
                          cov = cov,xNameMatrix = xNameMatrix, yName = yName,
                          foldid = foldid)

  # To get training weights


  MiXcan_weight_result <- MiXcan_extract_weight(model = MiXcan_result)
  W1 = matrix(MiXcan_weight_result[, 'weight_cell_1'])
  W2 = matrix(MiXcan_weight_result[, 'weight_cell_2'])
  if(!is.null(xNameMatrix)){
    selected_snp = MiXcan_weight_result$varID
    results <- lst(W1, W2, selected_snp, MiXcan_weight_result)
    return(results)
  }

  selected_snp = as.numeric(substr(MiXcan_weight_result$xNameMatrix, 4, nchar(MiXcan_weight_result$xNameMatrix)))
  results <- lst(W1, W2, selected_snp, MiXcan_weight_result)
  return(results)
}
