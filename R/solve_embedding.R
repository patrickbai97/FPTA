#' @useDynLib fpta
#' @import Rcpp
#' @import RcppArmadillo
#'
NULL

#'Solve Embedding
#'
#' This functions is the function which calculates the functional form of the disc game embedding via 3 steps.
#'  1. Gram-Schmidt 2. Calculate projection matrix 3. Calculate the sequence of disc games by Schur decomposition of the projection matrix.
#'  It gives the coefficients of every embedding constructed using user-supplied basis functions, as well as the relative importance of each embedding.
#'
#' @param f_val n x n skew symmetric matrix, with i,j th entry representing f(x_i, x_j)
#' @param list_basis a list of k vectorized function provided by the user
#' @param Xsample n samples, possible n x p matrix associated with f_val
#'
#' @return a list containing a vector lambda, which calculates the relative importance of the embedding, an 2p x k coefficient matrix, Mat_Embed_Coef
#' where 2i-1 and 2i row is the coefficients of the basis function associated with the ith embedding. The 2i-1 row is the embedding for the x-coordinate, the 2i th row is the embedding for the y-coordinate.
#' @export
#'
#' @examples M = matrix(c(1,1,1,1,0,1,1.5,2,0,1,2.25,4), nrow = 4, ncol = 3)
#' list_f = list(x0,x1,x2), where xn is the function x^n
#' embedding_list = solve_embedding(f_val, list_f, c(0,1,1.5,2))
#'
solve_embedding = function(f_val = NULL, list_basis = NULL, Xsample = NULL){
  #sanity check
  if(is.null(f_val) || is.null(list_basis) || is.null(Xsample)){
    stop("All inputs cannot be null.")
  }

  if(dim(f_val)[1] != dim(f_val)[2]){
    stop("f_val must be a skew-symmetric matrix.")
  }

  if(sum(abs(f_val + t(f_val))) > 0.001){
    stop("f_val must be a skew-symmetric matrix.")
  }

  if(dim(f_val)[1] != length(Xsample)){
    stop("The dimension of f_val and Xsample is not compatible.")
  }

  n = length(Xsample)
  k = length(list_basis)
  B = matrix(0, nrow = n, ncol = k)
  for(i in 1:length(list_basis)){
    B[,i] = list_basis[[i]](Xsample)
  }
  GramMatrix = gram_schimdtCpp(B)
  projection = update_projection(B, f_val, GramMatrix)
  embedding = update_embedding(projection, GramMatrix)

  return(embedding)
}
