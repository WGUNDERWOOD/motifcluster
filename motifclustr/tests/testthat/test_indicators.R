# build_indicator_matrices

test_that("build_indicator_matrices returns correct matrices", {

  G_dense = matrix(c(0,2,0,3,0,4,0,0,0), nrow=3, byrow=TRUE)
  G_sparse = drop0(G_dense)

  J  = drop0(matrix(c(0,1,0,1,0,1,0,0,0), nrow=3, byrow=TRUE))
  Gs = drop0(matrix(c(0,0,0,0,0,4,0,0,0), nrow=3, byrow=TRUE))
  Gd = drop0(matrix(c(0,5,0,5,0,0,0,0,0), nrow=3, byrow=TRUE))
  Js = drop0(matrix(c(0,0,0,0,0,1,0,0,0), nrow=3, byrow=TRUE))
  Jd = drop0(matrix(c(0,1,0,1,0,0,0,0,0), nrow=3, byrow=TRUE))
  J0 = drop0(matrix(c(0,0,1,0,0,0,1,0,0), nrow=3, byrow=TRUE))
  Jn = drop0(matrix(c(0,1,1,1,0,1,1,1,0), nrow=3, byrow=TRUE))

  expect_equal(build_J(G_dense), J)
  expect_equal(build_J(G_sparse), J)

  expect_equal(build_Gs(G_dense), Gs)
  expect_equal(build_Gs(G_sparse), Gs)

  expect_equal(build_Gd(G_dense), Gd)
  expect_equal(build_Gd(G_sparse), Gd)

  expect_equal(build_Js(G_dense), Js)
  expect_equal(build_Js(G_sparse), Js)

  expect_equal(build_Jd(G_dense), Jd)
  expect_equal(build_Jd(G_sparse), Jd)

  expect_equal(build_J0(G_dense), J0)
  expect_equal(build_J0(G_sparse), J0)

  expect_equal(build_Jn(G_dense), Jn)
  expect_equal(build_Jn(G_sparse), Jn)

})
