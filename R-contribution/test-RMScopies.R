library(CAs)
test_that("RMScopies works", {
  ## base construction (Raaphorst-Moura-Stevens 2014)
  D <- RMScopies(3)
  expect_equal(dim(D), c(53, 13))
  expect_true(caverify::ca_verify(D, 3)$covered)
  D <- RMScopies(5)
  expect_equal(dim(D), c(249, 31))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## one shifted copy
  D <- RMScopies(5, 1)
  expect_equal(dim(D), c(324, 32))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## several shifted copies (layers)
  for (m in 2:7){
    D <- RMScopies(5, m)
    expect_equal(nrow(D), N_RMScopiesCA(3, 31 + m, 5))
    expect_equal(ncol(D), 31 + m)
    expect_true(caverify::ca_verify(D, 3)$covered)
  }
  D <- RMScopies(7, 8)
  expect_equal(dim(D), c(980, 65))
  expect_true(caverify::ca_verify(D, 3)$covered)

  ## interface function
  E <- RMScopiesCA(37, 5)
  expect_s3_class(E, "ca")
  expect_equal(dim(E), c(350, 37))
  expect_equal(attr(E, "t"), 3)
  expect_equal(dim(RMScopiesCA(20, 5)), c(249, 20))

  ## run sizes
  expect_equal(N_RMScopiesCA(3, 31, 5), 249)
  expect_equal(N_RMScopiesCA(3, 32, 5), 324)
  expect_equal(N_RMScopiesCA(3, 36, 5), 325)
  expect_equal(N_RMScopiesCA(3, 37, 5), 350)
  expect_equal(N_RMScopiesCA(3, 58, 7), 930)
  expect_true(is.na(N_RMScopiesCA(4, 32, 5)))
  expect_true(is.na(N_RMScopiesCA(3, 32, 6)))

  ## invalid input
  expect_error(RMScopies(6))
  expect_error(RMScopies(5, 32))
  expect_error(RMScopiesCA(32, 6))
})
