#' Estimate the cell-type level prediction weights of a gene. ****** Sinan to update *******
#'
#' The core function of MiXcan package for estimating the
#' cell-type level prediction weights of a gene.
#' @param y The pre-cleaned expression level data for a single gene in N samples.
#' @param x A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi_k An estimation of cell-type fraction for the cell type of interest (e.g.
#' epithelial). It can be obtained using existing methods
#' in the literature or from the output of pi_estimation function.
#' @param xNameMatrix Default is NULL. A matrix to save theX matrix information,
#' such as variable ID, position, rsid, ref_allele, eff_allele.
#' @param yName Default is NULL. A row vector to save the expression information, such as gene ID, gene name.
#' @param foldid Default is NULL. 10-fold cross-validation (CV) is used in our pipeline. A random split
#' is considered if foldid is NULL. Otherwise foldid is used to split the data for CV.
#' @return list with 9 elements. It contains
#'
#' \item{type:}{Whether the prediction model is "CellTypeSpecific" or "NonSpecific".}
#' \item{beta.SNP.cell1:}{The prediction weights of the genetic predictors in cell type 1 (the cell type of interest).}
#' \item{beta.SNP.cell2:}{The prediction weights of the genetic predictors in cell type 2 (other cell types).}
#' \item{beta.all.models:}{All regression coefficients are saved in beta.all.models, including intercepts,
#' coefficients of genetic and non-genetic predictors in cell-type specific and non-specific models.}

#' \item{glmnet.cell:}{The cell-type-level prediction model selected using elastic-net.
#' This model may not be the final model as elastic-net selected parameters in the two
#' cell types may not be robustly different.  }
#' \item{glmnet.tissue:}{The tissue-level prediction model, which does not consider
#' cell type composition (same as PrediXcan).}
#' @export
#'
#'
#'
#'
# set.seed(123)
# n=200
# x=matrix(rnorm(n*10), ncol=10)
# y=rnorm(n)
# library(MCMCpack)
# pi_k=rdirichlet(n, alpha=c(1,1,1))
# cov=NULL # Sinan to check if it works with cov.
# xNameMatrix=NULL
# yName=NULL
# foldid=NULL

MiXcan_train_K=function(y, x, cov=NULL, pi_k, xNameMatrix=NULL, yName=NULL,
                foldid=NULL) {
  # harmonize input
  x=as.matrix(x); y=as.matrix(y);  pi_k=as.matrix(pi_k)
  n=nrow(x); p=ncol(x); K=ncol(pi_k)

  if(is.null(xNameMatrix)) {xNameMatrix=paste0("SNP", 1:p)}
  if (is.null(foldid)) {foldid= sample(1:10, n, replace=T)}

  # prepare variables for the model
  y=scale(y, center=T, scale=F)
  x=scale(x, center=T, scale=F)
  pi_k=scale(pi_k, center=T, scale=F)
  z=t(sapply(1:n, function(f) outer(x[f, ], pi_k[f, ])))
  z=scale(z, center=T, scale=F)
  # z_names=c(outer(xNameMatrix, 1:K, paste, sep = "_"))

  if (is.null(cov)) { xcov=x; xx=as.matrix(cbind(pi_k, z))}

  if (is.null(cov)==F) {
    cov=as.matrix(cov);cov=scale(cov, center=T, scale=F)
    xcov=as.matrix(cbind(x, cov))
    xx=as.matrix(cbind(pi_k, z, cov))
  }

  # tissue-level model
  ft00=glmnet::cv.glmnet(x=xcov, y=y,family="gaussian",  intercept=F,
                         foldid=foldid, alpha=0.5)
  ft0=glmnet::glmnet(x=xcov, y=y,  family="gaussian", intercept=F,
                     lambda = ft00$lambda.1se, alpha=0.5)
  est.tissue=as.numeric(ft0$beta)

  # cell-type-level model
  ft11=glmnet::cv.glmnet(x=xx, y=y, intercept=F,
                         family="gaussian", foldid=foldid, alpha=0.5)
  ft=glmnet::glmnet(x=xx, y=y, intercept=F,  family="gaussian",
                    lambda = ft11$lambda.1se, alpha=0.5)
  est=as.numeric(ft$beta)
  pi_coef=est[1:K]
  x_pi_coef=matrix(est[-c(1:K)],ncol=K, byrow = F)

  if (suppressWarnings(all(x_pi_coef==0))) {
    Type ="NoPredictor";
  } else {
    if (suppressWarnings(all( apply(x_pi_coef, 1, var)==0 ))) {
      Type ="NonSpecific"
    } else {Type ="CellTypeSpecific"}
  }


  ## add inference to ensure CellTypeSpecific are robust
  if (Type =="CellTypeSpecific") {
    idx.nonzero=which(apply(x_pi_coef, 1, mean)!=0)
    idx.nonzero.diff=which(apply(x_pi_coef, 1, mean)!=0 & apply(x_pi_coef, 1, var)!=0)
    x.diff=x_pi_coef[idx.nonzero.diff, ]
    x.diff=matrix(x.diff, ncol=K)
    idx.cell.max=apply(x.diff, 1, which.max )
    idx.cell.min=apply(x.diff, 1, which.min )
    # create xx.select

    z.select=t(sapply(1:n, function(f)
      outer(x[f, idx.nonzero], pi_k[f, ])))
    z.select=scale(z.select, center=T, scale=F)
    # z.select_names=c(outer( idx.nonzero, 1:K, paste, sep = "_"))

    if (is.null(cov)) {xx.select=as.matrix(cbind(pi_k, z.select))} else{
      xx.select=as.matrix(cbind(pi_k, z.select, cov))
    }

    beta.ols.boot=NULL; # beta.en.boot=NULL
    for(boot in 1:200) {
      id=sample(n, n, replace =T)
      gfit = glmnet::glmnet(x=xx.select[id,], y=y[id,],
                            alpha=0, lambda=0.00001, intercept=F)
      refit=matrix(gfit$beta[-c(1:K)], ncol=K, byrow=F)

      # For each SNP use the two cell types with the largest diff in the CTS model

      # Sinan, make this code works without matrix (No. of SNP =1 )
      # if (length(idx.nonzero.diff)== 1) {} else { XX }
      diff=sapply(idx.nonzero.diff, function(f) refit[f, idx.cell.max[f]])-
        sapply(idx.nonzero.diff, function(f) refit[f, idx.cell.min[f]])


      beta.ols.boot =rbind(beta.ols.boot, diff)
    }
    beta.diff.range=apply(beta.ols.boot, 2, function(f) quantile(f, prob=c(0.025, 0.975), na.rm=T))

    if (is.null(dim(beta.diff.range))) {
      any.nonzero= beta.diff.range[1] * beta.diff.range[2]>0
      } else {
        print(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0))
        any.nonzero= any(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0), na.rm=T)
      }

    if (any.nonzero==F) {
      Type ="NonSpecific";
      x_pi_coef =do.call(cbind, replicate(K, est.tissue, simplify = FALSE))

    }
  }
  colnames(x_pi_coef)=paste0("Cell", 1:K)
  rownames(x_pi_coef)=paste0("SNP", 1:p)
  x_pi_coef_all=cbind(x_pi_coef, tissue=est.tissue)


  return(list(type=Type,
              weight.matrix= x_pi_coef,
              weight.all.models=x_pi_coef_all,
              yName=yName,
              xNameMatrix=xNameMatrix))

}



































