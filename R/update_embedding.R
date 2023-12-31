#' @useDynLib fpta
#' @import Rcpp
#' @import RcppArmadillo
#'
NULL

#' This calculates the disc game embedding of the coefficients associated with the basis function provided by the user
#'
#' @param projection k x k projection matrix generated by update_projection
#' @param GramMatrix the matrix generated by Gram-Schmit process
#'
#' @return a list containing a vector lambda, which calculates the relative importance of the embedding, an 2p x k coefficient matrix, Mat_Embed_Coef
#' where 2i-1 and 2i row is the coefficients of the basis function associated with the ith embedding. The 2i-1 row is the embedding for the x-coordinate, the 2i th row is the embedding for the y-coordinate.
#'
update_embedding = function(projection, GramMatrix){
  SchurDecomp = Schur_wrapperCpp(projection)
  U = SchurDecomp$U
  S = SchurDecomp$S
  mcoef = matrix(0, nrow = dim(U)[1], ncol = dim(U)[2])
  lambda = rep(0, dim(U)[1]/2)
  for(i in 1:(dim(U)[1]/2)){
    lambda_temp = S[2*i - 1,2*i]
    lambda[i] = sqrt(abs(lambda_temp))
    if(lambda_temp > 0){
      mcoef[2*i - 1, ] = lambda[i] * U[, 2*i - 1]
      mcoef[2*i,] = lambda[i] * U[, 2*i]
    }else{
      mcoef[2*i - 1, ] = lambda[i] * U[, 2*i]
      mcoef[2*i, ] = lambda[i] * U[, 2*i - 1]
    }
  }
  k = dim(GramMatrix)[1]
  Mat_Embed_Coef = tcrossprod(mcoef[,1:k], GramMatrix)
  return(list(lambda = lambda, Mat_Embed_Coef = Mat_Embed_Coef))
}
