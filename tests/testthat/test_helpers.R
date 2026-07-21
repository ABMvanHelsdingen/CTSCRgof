test_that("empirical AC",{

  Cx <- seq(0.3, 0.75, length.out = 4)
  C <- matrix(0, nrow = 16, ncol = 2)
  C[,1] <- rep(Cx, times = 4)
  C[,2] <- rep(Cx, each = 4)

  capthist <- matrix(0, nrow = 5, ncol = 3)
  capthist[,2] <- c(1,1,1,2,2)
  capthist[,3] <- c(15,14,11,10,5)

  out <- empiricalAC(C, capthist)
  answer <- matrix(c(0.55,0.7,0.375,0.525), byrow=TRUE, nrow = 2)

  expect_equal(out, answer)
})

test_that("FT Stat",{

  obs <- c(4,9,16)
  exp <- c(9,16,25)

  expect_equal(FT_Stat(obs, exp), 3)

  obs <- matrix(c(4,16,9,1), nrow = 2)
  exp <- matrix(c(1,9,4,16), nrow = 2)

  expect_equal(FT_Stat(obs, exp), 12)
})

test_that("convertCaptHist",{

  ids <- c(1,3,4)
  N <- 4
  out <- convertCaptHist(ids, N)
  expect_equal(out$ids, c(1,2,3))
  expect_equal(out$new_order, c(1,3,4,2))

  ids <- c(2,5,3)
  N <- 6
  out <- convertCaptHist(ids, N)
  expect_equal(out$ids, c(1,3,2))
  expect_equal(out$new_order, c(2,3,5,1,4,6))

})
