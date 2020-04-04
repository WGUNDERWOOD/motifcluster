context("Motif adjacency")

# get_largest_component

test_that("get_largest_component returns correct indices", {

  G_dense = matrix(c(
    0,0,0,0,0,
    0,0,1,0,0,
    0,0,0,0,2,
    3,0,0,0,0,
    0,4,0,0,0
  ), nrow=5, byrow=TRUE)
  G_sparse = Matrix::drop0(G_dense)

  ans = c(2,3,5)

  expect_equal(get_largest_component(G_dense), ans)
  expect_equal(get_largest_component(G_sparse), ans)
})

# drop0_killdiag

test_that("drop0_killdiag returns correct matrix", {

  G_dense = matrix(-1:7, nrow=3)
  G_sparse = Matrix::drop0(G_dense)

  ans = Matrix::drop0(matrix(c(0,0,1,2,0,4,5,6,0), nrow=3))

  expect_equal(drop0_killdiag(G_dense), ans)
  expect_equal(drop0_killdiag(G_sparse), ans)
})
