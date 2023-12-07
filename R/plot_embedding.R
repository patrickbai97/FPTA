#' Plot Embedding
#'
#' This function plots the embedding for the chosen samples for better visualization.
#' The payout of two agents can be understand as the area of the triangle formed by the two points and the origin.
#'
#' @param list_basis a list of k vectorized function provided by the user, the basis function should be the same as used in solve_embedding()
#' @param Mat_Embed_Coef 2p x k coefficient matrix, output from solve_embedding()
#' @param X n x p matrix, containing the traits of n players
#' @param i integer, 1<= i <= p. Generate the plot for the i th embedding
#'
#' @return a plot
#' @export
#' @examples plot_embedding(list_basis, Mat_Embed_Coef, X, 1)
#'
plot_embedding = function(list_basis = NULL, Mat_Embed_Coef = NULL, X = NULL, i = 1){
  if(is.null(list_basis) || is.null(Mat_Embed_Coef) || is.null(X) || is.null(i)){
    stop("All inputs should not be null.")
  }

  if(i < 1){
   stop("i must be at least 1")
  }

  if(i > dim(Mat_Embed_Coef)[1]/2){
    stop("The inputted i is larger than the total number of embedding. Please enter a smaller i.")
  }

  if (length(list_basis) != dim(Mat_Embed_Coef)[2]){
    stop("The number of basis in list_basis not compatible with the embedding coefficient matrix Mat_Embed_Coef")
  }
  n = length(X)
  k = length(list_basis)
  B1 = matrix(0, nrow = n, ncol = k)
  for(j in 1:length(list_basis)){
    B1[,j] = list_basis[[j]](X)
  }
  x_coor = B1 %*% Mat_Embed_Coef[2*i-1,]
  y_coor = B1 %*% Mat_Embed_Coef[2*i,]
  embedding.plot = plot(x_coor, y_coor)
  return(embedding.plot)
}
