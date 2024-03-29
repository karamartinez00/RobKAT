% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{IBSKernel}
\alias{IBSKernel}
\title{Compute the IBS kernel}
\usage{
IBSKernel(z)
}
\arguments{
\item{\code{z}}{A SNP matrix of size \code{n x M}. Here \code{n} is the
number of subjects and \code{M} is the number of SNPs. Each
entry of \code{z} should be 0, 1 or 2 denoting the minor
allele frequency.}
}
\description{
Compute the IBS kernel
}
