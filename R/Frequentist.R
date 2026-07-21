#' Negative log-likelihood of a continuous-time SCR model
#'
#' Applies the generic likelihood formula for any Poisson SCR model where
#' individuals are independent
#' @inheritParams timeRescaling
#' @inheritParams hazard_SCR
#' @export

singleNLL <- function(s, C, duration, pars, hazard, capthist){

  times <- capthist[,1]
  cameras <- capthist[,3]
  I <- length(times)
  tau <- c(times, duration)
  out <- hazard(s, C, duration, pars, tau, capthist)

  # NLL is integral minus sum of log hazards at detections
  if(I > 1){
    hazards <- diag(out$h[1:I,cameras])
  } else if (I == 1){
    hazards <- out$h[1, cameras]
  } else {
    hazards <- 1
  }
  return(sum(out$H[I+1,]) - sum(log(hazards)))
}

#' Find the activity center for a single individual
#'
#' Applies Bayes' Rule to return a point estimate for an observed individual
#'
#' Note there are no constraints to ensure that the estimated ACs are within the survey region
#' @inheritParams timeRescaling
#' @param s starting point for the activity center
#' @export

findAC <- function(s, C, duration, pars, hazard, capthist){

  out = stats::optim(par = s, fn = singleNLL, gr = NULL, C, duration, pars, hazard, capthist)
  return(out$par)
}

#' Find the activity center for all individuals
#'
#' For observed individuals, use maximum likelihood estimation
#' For unobserved individuals, sample from the likelihood.
#' @inheritParams timeRescaling
#' @param S starting points for the activity centers of observed individuals

findACs <- function(S, C, duration, pars, hazard, N, mask, capthist){

  times <- capthist[,1]
  ids <- capthist[,2]
  cameras <- capthist[,3]

  m <- max(ids)

  centers <- matrix(0, nrow = N, ncol = 2)

  # Observed individuals
  for(ani in 1:m){
    Ij = sum(capthist[,2] == ani)
    capthist_ani <- matrix(0, nrow = Ij, ncol = 3)
    capthist_ani[1:Ij,] <- capthist[capthist[,2] == ani,]
    centers[ani,] <- findAC(S[ani,], C, duration, pars, hazard, capthist_ani)
  }

  # Do not proceed further if no unobserved individuals
  if (m == N){
    return(centers)
  }

  # Unobserved individuals
  nrand <- nrow(mask)
  lkhd <- numeric(nrand)
  empty <- matrix(0, nrow = 0, ncol = 3)
  for(pt in 1:nrand){
    lkhd[pt] <- exp(-sum(hazard(mask[pt,], C, duration, pars, duration, empty)$H))
  }

  lkhd <- lkhd/sum(lkhd)
  lkhd <- cumsum(lkhd)

  for(ani in (m+1):N){
    index <- min(which(lkhd > stats::runif(1,0,1)))
    centers[ani,] = mask[index,]
  }

  return(centers)
}

#' Time Rescaling Tests
#'
#' Performs cumulation and superposition time rescaling tests for a fitted frequentist SCR model
#'
#'
#' @inheritParams hazard_SCR
#' @inheritParams simulate_SCR
#' @param S matrix of activity centers. If not provided, `findACs()` will be used to estimate them
#' @param pars vector of parameter values for `hazard`
#' @param N estimated population size
#' @param hazard function that returns hazard rates, e.g. `hazard_SCR()`
#' @param comps (optional) output of `Hazards()` function.
#' @param mask (only needed if `S` not provided) points that lie within survey region and which serve as an approximation of it
#' @returns a list with four elements: \cr
#'  * `pvalues` a 2x2 dataframe containing p-values for the cumulation and superposition tests, with and without subsampling \cr
#'  * `S` a Nx2 matrix of activity centers used for the tests \cr
#'  * `CD_cumulation` the compensator differences from cumulation \cr
#'  * `CD_aggregation` the compensator differences from aggregation \cr
#' @export

timeRescaling <- function(S = NULL, C, duration, pars, hazard, N, comps = NULL,
                          mask = NULL, capthist){

  # Constants
  K <- nrow(C)
  I <- nrow(capthist)
  S1 <- N * K

  # Unpack capthist
  checkCaptHist(capthist, duration, K)
  times <- capthist[,1]
  ids <- capthist[,2]
  cameras <- capthist[,3]


  if (N %% 1 != 0){
    stop("N is the population size and must be integer")
  }


  if (is.null(S)){
    if(is.null(mask)){
      stop("If S is not supplied, a mask must be supplied")
    }
    if(ncol(mask) != 2){
      stop("mask must have two columns")
    }
    obsA <- empiricalAC(C, capthist)
    S <- findACs(obsA, C, duration, pars, hazard, N, mask, capthist)
  }

  # Constants
  K <- nrow(C)
  I <- nrow(capthist)
  S1 <- N * K


  # Calculate comps if necessary
  if (is.null(comps)){
    comps <- Hazards(S, C, duration, pars, times, hazard, N, capthist)
  }

  # Calculate compensator differences in each stream
  CD <- list(); counter = 0
  integrals <- numeric(S1)

  for(ani in 1:N){
    for(cam in 1:K){
      counter = counter + 1
      indices = which(cameras == cam & ids == ani)
      integrals[counter] = comps[[ani]][(I+1), cam]

      if (length(indices) > 0){
        CD[[counter]] = comps[[ani]][indices, cam]
      } else{
        CD[[counter]] = numeric(0)
      }
    }
  }

  # Calculate residuals for each stream
  res <- numeric(S1)
  for(s in 1:S1){

    l <- length(CD[[s]])
    if (l >= 1){
      res[s] <- integrals[s] - CD[[s]][l]
    } else {
      res[s] <- integrals[s]
    }

    if (l >= 2){
      CD[[s]][2:l] = CD[[s]][2:l] - CD[[s]][1:(l-1)]
    }

  }

  # Randomize order
  order = sample.int(S1)
  residual = 0; CD1 = numeric(0)
  for(s in order){
    l <- length(CD[[s]])
    if (l == 0){
      residual = residual + res[s]
    } else if (l == 1){
      CD1 = append(CD1, CD[[s]] + residual)
      residual = res[s]
    } else { # l >=2
      CD1 = append(CD1, CD[[s]][1] + residual)
      CD1 = append(CD1, CD[[s]][2:l])
      residual = res[s]
    }
  }

  # Perform superposition/aggregation
  CD2 <- numeric(I)
  for(ani in 1:N){
    CD2 <- CD2 + rowSums(comps[[ani]][1:I, ])
  }

  CD2 <- CD2 - c(0, CD2[1:(I-1)])


  pvalues <- data.frame(Cumulation = c(0), Aggregation = c(0),
                        row.names = c("Conventional"))

  # Cumulation Test
  pvalues[1, 1] <- stats::ks.test(CD1, "pexp", 1)$p.value

  # Superposition/Aggregation Tests
  pvalues[1, 2] <- stats::ks.test(CD2, "pexp", 1)$p.value


  return(list(pvalues = pvalues, S_hat = S,
              CD_cumulation = CD1, CD_aggregation = CD2))
}

#' Residuals for Frequentist SCR Models
#'
#' Calculates both trap and individual level residuals for frequentist SCR models
#'
#' @inheritParams timeRescaling
#' @returns a list with three elements: \cr
#'  * `S` a Nx2 matrix of activity centers used for the tests \cr
#'  * `resid_trap` the trap residuals \cr
#'  * `resid_ind` the individual residuals \cr
#' @export

getResiduals <- function(S = NULL, C, duration, pars, hazard, N, comps = NULL,
                      mask, capthist){

  # Constants
  K <- nrow(C)
  I <- nrow(capthist)
  S1 <- N * K

  # Unpack capthist
  checkCaptHist(capthist, duration, K)
  times <- capthist[,1]
  ids <- capthist[,2]
  cameras <- capthist[,3]
  m <- max(ids)


  if (N %% 1 != 0){
    stop("If S is not supplied, N must be the population size and integer")
  }


  if (is.null(S)){
    if(is.null(mask)){
      stop("If S is not supplied, a mask must be supplied")
    }
    if(ncol(mask) != 2){
      stop("mask must have two columns")
    }

    obsA <- empiricalAC(C, capthist)
    S <- findACs(obsA, C, duration, pars, hazard, N, mask, capthist)
  }

  # Calculate comps if necessary
  if (is.null(comps)){
    comps <- Hazards(S, C, duration, pars, times, hazard, N, capthist)
  }


  # Store information
  CD <- list()
  integrals <- numeric(K)

  # RESIDUALS BY TRAP
  for(k in 1:K){
    CD[[k]] = numeric(0)
    indices = which(cameras == k)
    if(length(indices) == 0){
      next
    }
    # For each event at camera, calculate compensator
    for(i in 1:length(indices)){
      integral <- 0
      for(ani in 1:N){
        integral <- integral + comps[[ani]][[indices[i], k]]
      }
      CD[[k]] <- append(CD[[k]], integral)
    }

    # Calculate compensator at end of survey
    for(ani in 1:N){
      integrals[k] = integrals[k] + comps[[ani]][(I+1), k]
    }

  }

  # For cameras with observations, apply corrected Time Rescaling
  p1 <- rep(NA, K)
  for(k in 1:K){
    l <- length(CD[[k]])
    if (l >= 2){
      maxs <- integrals[k] - c(0, CD[[k]][1:(l-1)])
      CD[[k]][2:l] = CD[[k]][2:l] - CD[[k]][1:(l-1)]
      CD[[k]] = (1 - exp(-CD[[k]])) / (1 - exp(-maxs))
      p1[k] = as.numeric(stats::ks.test(CD[[k]], "punif", 0, 1)$p.value)
    } else if (l == 1){
      CD[[k]] = (1 - exp(-CD[[k]])) / (1 - exp(-integrals[k]))
      p1[k] = as.numeric(stats::ks.test(CD[[k]], "punif", 0, 1)$p.value)
    }
  }

  # RESIDUALS BY INDIVIDUAL
  CD = list()
  for(ani in 1:m){
    CD[[ani]] = numeric(0)
    indices = which(ids == ani)
    for(i in 1:length(indices)){
      CD[[ani]] <- append(CD[[ani]], sum(comps[[ani]][indices[i], ]))
    }
  }

  # Calculate total integrals
  integrals = numeric(m)
  for(ani in 1:m){
    integrals[ani] <- sum(comps[[ani]][(I+1), ])
  }

  # For individuals with observations, apply corrected Time Rescaling
  p2 <- numeric(m)
  for(ani in 1:m){
    l <- length(CD[[ani]])
    if (l >= 2){
      maxs <- integrals[ani] - c(0, CD[[ani]][1:(l-1)])
      CD[[ani]][2:l] = CD[[ani]][2:l] - CD[[ani]][1:(l-1)]
      CD[[ani]] = (1 - exp(-CD[[ani]])) / (1 - exp(-maxs))
      p2[ani] = as.numeric(stats::ks.test(CD[[ani]], "punif", 0, 1)$p.value)
    } else {
      CD[[ani]] = (1 - exp(-CD[[ani]])) / (1 - exp(-integrals[ani]))
      p2[ani] = as.numeric(stats::ks.test(CD[[ani]], "punif", 0, 1)$p.value)
    }
  }

  return(list(S = S, resid_cam = p1, resid_ind = p2))
}
