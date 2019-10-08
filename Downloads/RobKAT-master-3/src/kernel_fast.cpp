#include <Rcpp.h>
using namespace Rcpp;

//' Compute the IBS kernel
//'
//' @param \code{z} A SNP matrix of size \code{n x M}. Here \code{n} is the
//'                   number of subjects and \code{M} is the number of SNPs. Each
//'                   entry of \code{z} should be 0, 1 or 2 denoting the minor
//'                   allele frequency.
//' @export
// [[Rcpp::export]]
NumericMatrix IBSKernel(NumericMatrix z) {
  int i, j, k, n, p, df;
  double tmp;

  n = z.nrow();
  p = z.ncol();

  NumericMatrix kernel(n, n);
  for(i = 0; i < n; i++){
    for(j = i+1; j < n; j++){
      tmp = 0;
      for(k = 0; k < p; k++){
        df = abs(z(i, k) - z(j, k));
        tmp += df;
      }
      kernel(i, j) = 1 - tmp / (2*p);
      kernel(j, i) = 1 - tmp / (2*p);
    }
  }

  for(i = 0; i < n; i++){
    kernel(i, i) = 1;
  }

  return kernel;
}
