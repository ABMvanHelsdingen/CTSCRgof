
test_that("Single NLL",{

  s <- c(0,0)
  C <- matrix(0, nrow = 4, ncol = 2)
  C[,1] <- c(0,1,0,1)
  C[,2] <- c(0,0,1,1)
  duration <- 4
  pars <- c(1,1)
  times <- c(1,2,3)
  tau <- c(times, duration)
  out <- hazard_SCR(s,C,duration,pars,tau,NULL)

  capthist <- matrix(0, nrow = 3, ncol = 3)
  capthist[,1] <- times
  capthist[,2] <- c(1,1,1)
  capthist[,3] <- c(1,2,1)

  # Test multiple detections
  NLL <- sum(out$H[4,]) - log(out$h[1,1]) - log(out$h[2,2]) - log(out$h[3,1])
  expect_equal(singleNLL(s, C, duration, pars, hazard_SCR, capthist), NLL)

  # Test one detection
  NLL <- sum(out$H[4,]) - log(out$h[1,1])
  expect_equal(singleNLL(s, C, duration, pars, hazard_SCR, matrix(capthist[1,],nrow=1)), NLL)

  # Test no detection
  NLL <- sum(out$H[4,])
  expect_equal(singleNLL(s, C, duration, pars, hazard_SCR, matrix(0,nrow=0,ncol=3)), NLL)


})

test_that("Residuals",{

  s <- c(0,0)
  C <- matrix(0, nrow = 4, ncol = 2)
  C[,1] <- c(0,1,0,1)
  C[,2] <- c(0,0,1,1)
  duration <- 4
  pars <- c(1,1)
  times <- c(1,2,3,4)
  tau <- c(times, duration)
  out <- hazard_SCR(s,C,duration,pars,tau,NULL)
  R <- sum(out$h[1,])


  capthist <- matrix(0, nrow = 4, ncol = 3)
  capthist[,1] <- times
  capthist[,2] <- c(1,1,1,1)
  capthist[,3] <- c(1,2,1,2)

  out <- getResiduals(matrix(s,nrow=1), C, 4, pars, hazard_SCR, 1, NULL, NULL, capthist)

  # Ind
  vec <- (1 - exp(-rep(R,4))) / (1 - exp(-c(4*R,3*R,2*R,R)))
  expect_equal(out$resid_ind, ks.test(vec, "punif", 0, 1)$p.value)

  # Trap 1
  vec <- (1 - exp(-c(1,2))) / (1 - exp(-c(4,3)))
  expect_equal(out$resid_cam[1], ks.test(vec, "punif", 0, 1)$p.value)

  # Trap 2
  vec <- (1 - exp(-2*c(exp(-0.5), exp(-0.5)))) / (1 - exp(-c(4*exp(-0.5),2*exp(-0.5))))
  expect_equal(out$resid_cam[2], ks.test(vec, "punif", 0, 1)$p.value)

})
