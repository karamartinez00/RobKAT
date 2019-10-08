context("Is C++ IBS same as R IBS?")


IBS.kernel <- function(X)
{
  ## X is the full nxp SNP matrix, n: sample size, m: number of variables
  ## output: nxn kernel matrix
  n = nrow(X)
  K = matrix(nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in i:n) {
      K[i,j] = K[j,i] = sum(  2*(X[j,] == X[i,]) + (abs(X[j,]-X[i,])==1) )
    }
  }
  K = K/(2*ncol(X))
  return(K)
}

test_that("IBS with M = 9, n = 10",{
  z <- matrix( sample(0:2, size = 90, replace = T), ncol = 9 )
  expect_equal(IBS.kernel(z), IBSKernel(z))
  }
  )


test_that("IBS with M = 90, n = 100",{
  z <- matrix( sample(0:2, size = 9000, replace = T), ncol = 9 )
  expect_equal(IBS.kernel(z), IBSKernel(z))
}
)

test_that("IBS with M = 90, n = 10",{
  z <- matrix( sample(0:2, size = 900, replace = T), ncol = 9 )
  expect_equal(IBS.kernel(z), IBSKernel(z))
}
)

