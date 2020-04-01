context("Spectral methods")
library(RSpectra)

# get_first_eigs

test_that("get_first_eigs returns correct values on dense matrix", {

  G = matrix(c(7,-4,14,0,-4,19,10,0,14,10,10,0,0,0,0,100), nrow=4)
  vals = c(-9, 18, 27)
  vects = matrix(c(-2,-1,2,0,-2,2,-1,0,-1,-2,-2,0)/3, nrow=4)

  spect = get_first_eigs(G, 3)

  expect_equal(spect$vals, vals)
  expect_equal(spect$vects, vects)
})

test_that("get_first_eigs returns correct values on sparse matrix", {

  G = as(matrix(c(7,-4,14,0,-4,19,10,0,14,10,10,0,0,0,0,100), nrow=4), "sparseMatrix")
  print(G)
  vals = c(-9, 18, 27)
  vects = matrix(c(-2,-1,2,0,-2,2,-1,0,-1,-2,-2,0)/3, nrow=4)

  spect = get_first_eigs(G, 3)

  expect_equal(spect$vals, vals)
  expect_equal(spect$vects, vects)
})
