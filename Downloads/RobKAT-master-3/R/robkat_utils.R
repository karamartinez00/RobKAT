#' @useDynLib RobKAT
#' @importFrom Rcpp sourceCpp
NULL

# Simple trace
tr <- function(x){
  return( sum(diag(x)) )
  }

# Moments of the RV stat
KRV=function(K,L){
  n=nrow(K)
  A=scale(K ,scale = F) ## that is A=HK
  W=scale(L, scale = F) ## that is W=HL
  Fstar=tr(A%*%W)
  mean.krv=tr(A)*tr(W)/(n-1)  ## mean of KRV

  T=tr(A)
  T2=tr(A%*%A)
  S2=sum(diag(A)^2)

  Ts=tr(W)
  T2s=tr(W%*%W)
  S2s=sum(diag(W)^2)

  temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
  temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
  temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
  temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
  temp2=temp21*temp22/temp23
  variance.krv=temp1+temp2        ## variance of KRV

  T3=tr(A%*%A%*%A)
  S3=sum(diag(A)^3)
  U=sum(A^3)
  R=t(diag(A))%*%diag(A%*%A)
  B=t(diag(A))%*%A%*%diag(A)

  T3s=tr(W%*%W%*%W)
  S3s=sum(diag(W)^3)
  Us=sum(W^3)
  Rs=t(diag(W))%*%diag(W%*%W)
  Bs=t(diag(W))%*%W%*%diag(W)

  t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
  t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
  t3=24*(n^2-n-4)*(U*Bs+B*Us)
  t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
  t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
  t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
  t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
  t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
  t8=24*(t81+t82)
  t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
  t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
  t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
  t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
  t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
  t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
  t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
  t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
  t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
  t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
  t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
  t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
  t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
  t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
  t20=-(n-2)*(t201+t202+t203)
  temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
  temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
  mom3=temp31/temp32
  skewness.krv = (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
  m1=mean.krv
  m2=variance.krv
  m3=skewness.krv
  shape=4/m3^2
  scale=sqrt(m2)*m3/2
  location = m1-2*sqrt(m2)/m3
  PIIIpars = list(shape, location, scale)
  pv = 1 - PearsonDS::ppearsonIII(Fstar, params = PIIIpars)
  return(pv)
}


# Psi Functions --------------------------------
# Huber
huberPsiF = function(u, k=1.345){
  u*(abs(u)<=k) + k*sign(u)*(abs(u)>k)
}

# Tukey's Bi-square
bisquarePsiF <- function(u, k=4.685){
  u * (1-(u/k)^2)^2 * (abs(u)<=k)
}

# Hampel
hampelPsiF <- function(u, a=1.353, b=3.157, r=7.216){
  #For 95% efficiency of the regression estimator
  # a = 1.5, b=3.5k, r=8k, k=0.902
  u*(abs(u)<=a) +
    a*sign(u)*((abs(u) > a) & (abs(u) <=b)) +
    a*sign(u)*((r-abs(u))/(r-b))*((abs(u) > b) & (abs(u) <=r))
}

# LAD
LADPsiF <- function(u){

  w <- 0.5*sign(u)
  zeroRes = sum(u == 0)
  w[u == 0] <- stats::rbinom(zeroRes, 1, 0.5) - 0.5
  w
}


#' Function for checking whether minor allele frequencies of SNPs are
#'                 correctly supplied
#'
#' @param \code{Z} An nxM matrix of genetic covariates to be tested. Each column
#'               corresponds to one SNP.
#' @param \code{silent} Whether to print a warning if genotypes are flipped
#' @details For each column of \code{Z}, we compute [2*(# of 2's) +
#'              (# of 1's)] and  [2*(# of 0's) + (# of 1's)]. If the former
#'              is smaller than the latter quantity, the MAF is correct. Otherwise,
#'              genotypes of that column are exchanged (0 -> 2 and 2 -> 0)
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
                 " are flipped to minor allels. ")
                 )
      }

    return(Znew)
  } else {
    return(Z)
  }

}
