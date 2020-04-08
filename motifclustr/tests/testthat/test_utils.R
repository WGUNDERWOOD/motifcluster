context("Utilities")

# a_b_one

test_that("a_b_one returns correct matrix",{

  A_dense = matrix(-4:4, nrow=3, byrow=TRUE)
  B_dense = matrix(-1:7, nrow=3, byrow=TRUE)
  A_sparse = drop0(A_dense)
  B_sparse = drop0(B_dense)

  ans = drop0(matrix(c(0,0,0,-9,0,9,36,54,72), nrow=3, byrow=TRUE))

  expect_equal(a_b_one(A_dense, B_dense), ans)
  expect_equal(a_b_one(A_sparse, B_sparse), ans)

})

# a_one_b

test_that("a_one_b returns correct matrix",{

  A_dense = matrix(-4:4, nrow=3, byrow=TRUE)
  B_dense = matrix(-1:7, nrow=3, byrow=TRUE)
  A_sparse = drop0(A_dense)
  B_sparse = drop0(B_dense)

  ans = drop0(matrix(c(-24,-27,-24,-6,0,12,12,27,48), nrow=3, byrow=TRUE))

  expect_equal(a_one_b(A_dense, B_dense), ans)
  expect_equal(a_one_b(A_sparse, B_sparse), ans)

})

# drop0_killdiag

test_that("drop0_killdiag returns correct matrix", {

  adj_mat_dense = matrix(-1:7, nrow=3)
  adj_mat_sparse = drop0(adj_mat_dense)

  ans = drop0(matrix(c(0,0,1,2,0,4,5,6,0), nrow=3))

  expect_equal(drop0_killdiag(adj_mat_dense), ans)
  expect_equal(drop0_killdiag(adj_mat_sparse), ans)
})

# get_largest_component

test_that("get_largest_component returns correct indices", {

  adj_mat_dense = matrix(c(
    0,0,0,0,0,
    0,0,1,0,0,
    0,0,0,0,2,
    3,0,0,0,0,
    0,4,0,0,0
  ), nrow=5, byrow=TRUE)
  adj_mat_sparse = drop0(adj_mat_dense)

  ans = c(2,3,5)

  expect_equal(get_largest_component(adj_mat_dense), ans)
  expect_equal(get_largest_component(adj_mat_sparse), ans)
})
