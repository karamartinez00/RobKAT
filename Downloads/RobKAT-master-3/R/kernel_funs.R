# # IBS kernel
# IBS.kernel = function(X)
# {
#   ## X is the full nxp SNP matrix, n: sample size, m: number of variables
#   ## output: nxn kernel matrix
#   n = nrow(X)
#   K = matrix(nrow = n, ncol = n)
#
#   for (i in 1:n) {
#     for (j in i:n) {
#       K[i,j] = K[j,i] = sum(  2*(X[j,] == X[i,]) + (abs(X[j,]-X[i,])==1) )
#     }
#   }
#   K = K/(2*ncol(X))
#   return(K)
# }
# attr(IBS.kernel, "name") <- "IBS.kernel"


#' Compute the linear kernel
#'
#' @param \code{Z} A SNP matrix of size \code{n x p}. Here \code{n} is the
#'                   number of subjects and \code{p} is the number of SNPs. Each
#'                   entry of \code{Z} should be 0, 1 or 2 denoting the minor
#'                   allele frequency.
#' @export
linearKernel = function(Z){
  Z %*% t(Z)
}

#' Compute the quadratic kernel
#'
#' @param \code{Z} A SNP matrix of size \code{n x M}. Here \code{n} is the
#'                   number of subjects and \code{M} is the number of SNPs. Each
#'                   entry of \code{Z} should be 0, 1 or 2 denoting the minor
#'                   allele frequency.
#' @export
quadraticKernel = function(Z){
  (1 + Z %*% t(Z))^2
}
