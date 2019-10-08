#' @useDynLib RobKAT
#' @importFrom Rcpp sourceCpp
NULL





#' Huber's psi function
#'
#' @param \code{u} An nx1 vector
#' @param \code{k} Cut-off parameter
#' @details Huber loss behaves like least squares loss until a cut-off point, k,
#'           after which the loss behaves like absolute value loss. 
#'           The default value, k = 1.345, is chosen for 95 \% efficiency 
#'           compared to the least squares estimate when the error distribution is normal. 
#'           Huber's psi function is defined as
#'           psi(u) = u*[abs(u) <= k] + k*sign(u)*[abs(u) > k]. It 
#'           is the derivative of rho(x) = 0.5*x^2*[abs(x)<= k] + k*(abs(x)-0.5*k)*[abs(x) > k]
#'
#' @export
#' @examples
#' u <- seq(-10, 10, len = 101)
#' plot(u, huberPsiF(u), type="b")
huberPsiF = function(u, k=1.345){
  u*(abs(u)<=k) + k*sign(u)*(abs(u)>k)
}

#' Tukey's Bi-square psi function
#'
#' @param \code{u} An nx1 vector
#' @param \code{k} Cut-off parameter
#' @details Tukey's Bi-square psi function is defined as
#'           psi(u) = u * (1-(u/k)^2)^2 * (abs(u)<=k)
#'
#' @export
#' @examples
#' u <- seq(-10, 10, len = 101)
#' plot(u, bisquarePsiF(u), type="b")
bisquarePsiF <- function(u, k=4.685){
  u * (1-(u/k)^2)^2 * (abs(u)<=k)
}


#' Hampel's psi function
#'
#' @param \code{u} An nx1 vector
#' @param \code{k} Cut-off parameter
#' @details Tukey's Bi-square psi function is defined as
#'           psi(u) = u*(abs(u)<=a) +
#'           a*sign(u)*((abs(u) > a) & (abs(u) <=b)) +
#'           a*sign(u)*((r-abs(u))/(r-b))*((abs(u) > b) & (abs(u) <=r))
#'           
#'           The default values of a=1.353, b=3.157, and r=7.216 are 
#'           chosen for 95 \% efficiency compared to the least squares 
#'           estimate when the error distribution is normal.
#'
#' @export
#' @examples
#' u <- seq(-10, 10, len = 101)
#' plot(u, hampelPsiF(u), type="b")
hampelPsiF <- function(u, a=1.353, b=3.157, r=7.216){
  #For 95% efficiency of the regression estimator
  # a = 1.5, b=3.5k, r=8k, k=0.902
  u*(abs(u)<=a) +
    a*sign(u)*((abs(u) > a) & (abs(u) <=b)) +
    a*sign(u)*((r-abs(u))/(r-b))*((abs(u) > b) & (abs(u) <=r))
}

#' LAD Check function
#'
#' @param \code{u} An nx1 vector
#' @param \code{k} Cut-off parameter
#' @details The check function for LAD is defined as
#'         psi(u) = sign(u)/2.
#'         If u  is zero, then psi(u) is assigned either -1/2 or 1/2 with
#'         probability 1/2. LAD loss represents the subgradient of rho(x) = abs(x).
#'
#' @export
#' @examples
#' u <- seq(-10, 10, len = 101)
#' plot(u, LADPsiF(u), type="b")
LADPsiF <- function(u){

  w <- 0.5*sign(u)
  zeroRes = sum(u == 0)
  w[u == 0] <- stats::rbinom(zeroRes, 1, 0.5) - 0.5
  w
}


#' Function for checking whether minor allele frequencies of SNPs are
#'                 correctly supplied
#'
#' @param \code{Z} An nxp matrix of genetic covariates to be tested. Each column
#'               corresponds to one SNP.
#' @param \code{silent} Whether to print a warning if genotypes are flipped
#' @details For each column of \code{Z}, we compute [2*(# of 2's) +
#'              (# of 1's)] and  [2*(# of 0's) + (# of 1's)]. If the former
#'              is smaller than the latter quantity, the MAF is correct. Otherwise,
#'              genotypes of that column are exchanged (0 -> 2 and 2 -> 0)
#' @export
#' @examples
#' n <- 100 #sample size
#' p <- 5 #no. of SNPs variables
#' Z <- matrix(sample(0:2, n*p, replace = T), ncol=p) # pretending to be SNP data
#' Znew <- check.Z.MinorAllels(Z, silent = FALSE)

check.Z.MinorAllels <- function(Z, silent = FALSE){

  n <- nrow(Z)
  p <- ncol(Z)

  a <- 2*colSums(Z == 2) + colSums(Z == 1) # count of a (supposed minor allele)
  A <- 2*colSums(Z == 0) + colSums(Z == 1) # count of A

  if(any(a + A != 2*n)){ stop("Z has values other than 0/1/2.") }

  # a should be always smaller than A. Thus this indices (when TRUE),
  # corresponds to the columns whoose genotypes need to be flipped
  flip.index <- (a > A)

  if(!all(!flip.index)){
    Znew <- Z
    Znew[ ,flip.index] <- 2 - Znew[ ,flip.index]

    if(!silent){
      wstr <- warning( paste("Genotypes of columns ",
                             paste(paste(c(1:p)[flip.index], sep=",")),
                             " are flipped to minor alleles. ")
      )
    }

    return(Znew)
  } else {
    return(Z)
  }

}

