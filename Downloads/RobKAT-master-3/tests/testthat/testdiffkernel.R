context("Different kernels")


n <- 100 #sample size
p <- 5 #no. of nonpar variables
q <- 2 #no. of parametric variables
bet = rep(1, q) #params for parametric effects
X <- matrix(sample(0:2, n*p, replace = T), ncol=p) # pretending to be SNP data
Z <- matrix(rnorm(n*q),n,q) # some other covariates
Y <- Z%*%bet + 0.5*rowSums(X) + rcauchy(n)

test_that("IBS",{
  expect_is(robkat(Y, X, Z, kernel.fun = "IBS", loss.fun = "huber"), "numeric")
  }
  )


test_that("linear",{
  expect_is(robkat(Y, X, Z, kernel.fun = "linear", loss.fun = "huber"), "numeric")
}
)


test_that("quadratic",{
  expect_is(robkat(Y, X, Z, kernel.fun = "quadratic", loss.fun  = "huber"), "numeric")
}
)


test_that("IBS",{
  expect_is(robkat(Y, X, Z, kernel.fun = "IBS", loss.fun = "LAD"), "numeric")
}
)


test_that("linear",{
  expect_is(robkat(Y, X, Z, kernel.fun = "linear", loss.fun = "LAD"), "numeric")
}
)


test_that("quadratic",{
  expect_is(robkat(Y, X, Z, kernel.fun = "quadratic", loss.fun  = "LAD"), "numeric")
}
)
