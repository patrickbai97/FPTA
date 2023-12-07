#' @useDynLib fpta
#' @import Rcpp
#' @import RcppArmadillo
#'
NULL

#' Predict the payout function f(x,y) with first person's trait = x and second person's trait as y
#'
#' @param list_basis a list of k vectorized function provided by the user, the basis function should be the same as used in solve_embedding()
#' @param Mat_Embed_Coef 2p x k coefficient matrix, output from solve_embedding()
#' @param X1 n x p matrix, containing the traits of first n players
#' @param X2 n x p matrix, containing the traits of second n players
#'
#' @return y, vector of length n, which is the estimated f(x,y)
#' @export
#'
#' @examples a = solve_embedding(f_val, list_f, c(0,1,1.5,2))
#' fpta_predict(list_f, a$Mat_Embed_Coef, rnorm(100), rnorm(100))
fpta_predict = function(list_basis = NULL, Mat_Embed_Coef = NULL, X1 = NULL, X2 = NULL){
  #sanity check
  if(is.null(list_basis) || is.null(Mat_Embed_Coef) || is.null(X1) || is.null(X2)){
    stop("All inputs should not be null.")
  }

  if (length(list_basis) != dim(Mat_Embed_Coef)[2]){
    stop("The number of basis in list_basis not compatible with the embedding coefficient matrix Mat_Embed_Coef")
  }

  if(length(X1) != length(X2)){
    stop("X1 and X2 should have equal length")
  }

  n = length(X1)
  k = length(list_basis)
  B1 = matrix(0, nrow = n, ncol = k)
  B2 = matrix(0, nrow = n, ncol = k)
  for(i in 1:length(list_basis)){
    B1[,i] = list_basis[[i]](X1)
    B2[,i] = list_basis[[i]](X2)
  }
  y = predictCpp(B1, B2, Mat_Embed_Coef)

  return(y)
}
