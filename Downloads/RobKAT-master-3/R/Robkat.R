#' Robust Kernel Association Test
#'
#' @param \code{Y} An nx1 response vector
#' @param \code{Z} An nxM matrix of genetic covariates to be tested
#' @param \code{X} An nxp matrix of non-genetic covariates, if present.
#'                  If no other covariates, set X = NULL. Do not
#'                  include an intercept column.
#' @param \code{kernel.fun} A string ("linear", "quadratic", "IBS")
#'                          indicating which kernel function to use
#' @param \code{loss.fun} A string ("huber", "hampel", "bisquare")
#'                          indicating which loss function to use
#' @param \code{silent} Whether to print a warning if genotypes are
#'                        flipped when the MAFs are checked for validity.
#' @keywords Robust
#' @export
#' @examples
#' n <- 100 #sample size
#' p <- 5 #no. of SNPs variables
#' q <- 2 #no. of parametric variables
#' bet = rep(1, q) #params for parametric effects
#' Z <- matrix(sample(0:2, n*p, replace = T), ncol=p) # pretending to be SNP data
#' X <- matrix(rnorm(n*q), n, q) # some other covariates
#' Y <- X%*%bet + 0.5*rowSums(Z) + rcauchy(n)
#' robkat(Y, Z, X, kernel.fun = "IBS", loss.fun = "huber")
#' robkat(Y, Z, X, kernel.fun = "IBS", loss.fun = "LAD")
robkat <- function(Y, Z, X = NULL,
                   kernel.fun = "linear", loss.fun = "huber",
                   silent = FALSE){

  n <- length(Y)

  # Check whether the MAFs are properly provided. If not, adjust accordingly
  Znew <- check.Z.MinorAllels(Z, silent = silent)

  # Kernel
  K <- do.call(base::paste(kernel.fun, "Kernel", sep=""), list(Znew))

  # Check if X iss NULL. If so, only setup the intercept
  Xfull <- cbind(1, X)
  if( is.null(X) ){
    Xfull <- rep(1, n)
  }


  # Psi function for robust fitting
  psi.par <- base::paste("psi.", loss.fun, sep="")
  psi.fun <- base::match.fun(paste(loss.fun, "PsiF", sep=""))

  if(loss.fun == "LAD"){
    # LAD regression
    quantresult <- quantreg::rq(Y ~ Xfull - 1, tau = 0.5) #estimating the median t=0.5

    # Residuals (u)
    ustd <- quantresult$residuals
  } else{
    # Robust regression
    rob = MASS::rlm(Y ~ Xfull - 1, psi = psi.par, maxit = 1000, scale.est = "Huber") ###

    # Standardized residuals (u)
    ustd <- (Y - rob$fitted.values) / rob$s
  }

  # w = Psi(u)
  w <- psi.fun(ustd)

  # P-value
  Kw <- w %*% t(w)
  pv <- KRV(Kw, K)

  c(pv)
}
