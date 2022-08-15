context("Indicator matrices")

test_that("indicator matrices return correct matrices", {

  G_dense <- matrix(c(0, 2, 0, 3, 0, 4, 0, 0, 0), nrow = 3, byrow = TRUE)
  G_sparse <- drop0(G_dense)

  G  <- drop0(matrix(c(0, 2, 0, 3, 0, 4, 0, 0, 0), nrow = 3, byrow = TRUE))
  J  <- drop0(matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 0), nrow = 3, byrow = TRUE))
  Gs <- drop0(matrix(c(0, 0, 0, 0, 0, 4, 0, 0, 0), nrow = 3, byrow = TRUE))
  Js <- drop0(matrix(c(0, 0, 0, 0, 0, 1, 0, 0, 0), nrow = 3, byrow = TRUE))
  Gd <- drop0(matrix(c(0, 5, 0, 5, 0, 0, 0, 0, 0), nrow = 3, byrow = TRUE))
  Jd <- drop0(matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), nrow = 3, byrow = TRUE))
  J0 <- drop0(matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow = 3, byrow = TRUE))
  Jn <- drop0(matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3, byrow = TRUE))
  Id <- Diagonal(3)
  Je <- drop0(matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), nrow = 3, byrow = TRUE))
  Gp <- drop0(matrix(c(0, 6, 0, 6, 0, 0, 0, 0, 0), nrow = 3, byrow = TRUE))

  expect_true(all(G - build_G(G_dense) == 0))
  expect_true(all(G - build_G(G_sparse) == 0))

  expect_true(all(J - build_J(G_dense) == 0))
  expect_true(all(J - build_J(G_sparse) == 0))

  expect_true(all(Gs - build_Gs(G_dense) == 0))
  expect_true(all(Gs - build_Gs(G_sparse) == 0))

  expect_true(all(Js - build_Js(G_dense) == 0))
  expect_true(all(Js - build_Js(G_sparse) == 0))

  expect_true(all(Gd - build_Gd(G_dense) == 0))
  expect_true(all(Gd - build_Gd(G_sparse) == 0))

  expect_true(all(Jd - build_Jd(G_dense) == 0))
  expect_true(all(Jd - build_Jd(G_sparse) == 0))

  expect_true(all(J0 - build_J0(G_dense) == 0))
  expect_true(all(J0 - build_J0(G_sparse) == 0))

  expect_true(all(Jn - build_Jn(G_dense) == 0))
  expect_true(all(Jn - build_Jn(G_sparse) == 0))

  expect_true(all(Id - build_Id(G_dense) == 0))
  expect_true(all(Id - build_Id(G_sparse) == 0))

  expect_true(all(Je - build_Je(G_dense) == 0))
  expect_true(all(Je - build_Je(G_sparse) == 0))

  expect_true(all(Gp - build_Gp(G_dense) == 0))
  expect_true(all(Gp - build_Gp(G_sparse) == 0))

})
