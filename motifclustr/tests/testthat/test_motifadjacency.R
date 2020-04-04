context("Motif adjacency")

# build_motif_adjacency_matrix
# test against igraph too?

test_that("build_motif_adjacency_matrix returns correct unweighted functional matrix", {

  G_dense = matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    2, 0, 3, 0, 6, 8, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 5, 0, 0, 0,14, 0, 0,18,19, 0, 0,
    0, 7, 9, 0,13, 0, 0, 0, 0, 0,21, 0,
    0, 0,11,12, 0,15, 0,17, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,16, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,24, 0, 0,
    0, 0, 0, 0, 0,20, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,22, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,23, 0, 0, 0, 0, 0,
  ), nrow=12, byrow=TRUE)

})

# build_indicator_matrices

test_that("build_indicator_matrices returns correct matrices", {

  G_dense = matrix(c(0,2,0,3,0,4,0,0,0), nrow=3, byrow=TRUE)
  G_sparse = Matrix::drop0(G_dense)

  ans = list()

  ans$G  = unname(Matrix::drop0(matrix(c(0,2,0,3,0,4,0,0,0), nrow=3, byrow=TRUE)))
  ans$J  = unname(Matrix::drop0(matrix(c(0,1,0,1,0,1,0,0,0), nrow=3, byrow=TRUE)))
  ans$J0 = unname(Matrix::drop0(matrix(c(0,0,1,0,0,0,1,0,0), nrow=3, byrow=TRUE)))
  ans$Jn = unname(Matrix::drop0(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3, byrow=TRUE)))
  ans$Gs = unname(Matrix::drop0(matrix(c(0,0,0,0,0,4,0,0,0), nrow=3, byrow=TRUE)))
  ans$Gd = unname(Matrix::drop0(matrix(c(0,5,0,5,0,0,0,0,0), nrow=3, byrow=TRUE)))
  ans$Js = unname(Matrix::drop0(matrix(c(0,0,0,0,0,1,0,0,0), nrow=3, byrow=TRUE)))
  ans$Jd = unname(Matrix::drop0(matrix(c(0,1,0,1,0,0,0,0,0), nrow=3, byrow=TRUE)))

  expect_equal(build_indicator_matrices(G_dense), ans)
  expect_equal(build_indicator_matrices(G_sparse), ans)
})

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
