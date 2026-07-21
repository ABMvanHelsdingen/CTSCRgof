#' Chain metrics
#'
#' As discussed in Appendix A. Can serve as the `metric` input to `runPPC()`.
#'
#' `data`, `N`, `K` and `duration` are not required for this metric
#'
#' @inheritParams runPPC
#' @inheritParams Hazards
#' @inheritParams checkCaptHist
#' @param data pre-calculated information needed to compute metric
#' @export

chainMetrics <- function(data = NULL, capthist,
                         N = NULL, K = NULL, duration = NULL){

  # Calculates and returns each of the 3 chain metrics

  output <- numeric(3)

  ids <- capthist[,2]
  cameras <- capthist[,3]

  if (length(ids) == 0){
    return(output)
  }

  m <- length(unique(ids))
  chains <- numeric(m)

  for(ani in unique(ids)){
    cameras_ani <- cameras[ids == ani]
    chains[ani] <- max(rle(cameras_ani)$lengths)
  }

  output[1] <- max(chains)
  output[2] <- mean(chains)
  output[3] <- stats::median(chains)
  return(output)
}

#' FT Metric Preparation
#'
#' FT metrics discussed in Choo et al, 2024 and other works.
#' Can serve as the `prep` input to `runPPC()`.
#'
#' @inheritParams timeRescaling
#' @export

FTMetrics_prep <- function(S, C, duration, pars,
                      hazard, capthist){


  # Get Expected counts across entire survey
  expected <- Integrals(S, C, duration, pars, hazard, nrow(S), capthist)

  return(expected)
}

#' FT metrics
#'
#' FT metrics discussed in Choo et al, 2024 and other works.
#' Can serve as the `metric` input to `runPPC()`.
#'
#' `duration` is not required for this metric
#' @inheritParams chainMetrics
#' @export

FTMetrics <- function(data, capthist, N, K, duration = NULL){

  # Make table
  # calculate expectations

  output <- numeric(3)

  ids <- capthist[,2]
  cameras <- capthist[,3]

  I <- length(ids)

  # Observed counts
  tab <- matrix(0, nrow = N, ncol = K)
  for(ani in 1:N){
    for(k in 1:K){
      tab[ani, k] <- length(which(ids == ani & cameras == k))
    }
  }

  output[1] <- FT_Stat(tab, data)
  output[2] <- FT_Stat(rowSums(tab), rowSums(data))
  output[3] <- FT_Stat(colSums(tab), colSums(data))

  return(output)
}

#' FT Transition Metric Preparation
#'
#' For standard SCR only. Transition matrix metric discussed in Appendix A
#' Can serve as the `prep` input to `runPPC()`.
#'
#' @inheritParams timeRescaling
#' @export

transMetric_prep_SCR <- function(S, C, duration, pars,
                        hazard, capthist){

  # Constants
  N <- nrow(S)
  K <- nrow(C)

  # Get expected counts by ind and trap
  expected <- Integrals(S, C, duration, pars, hazard, N, capthist)

  # Create matrix of expected transition counts
  trans_exp <- matrix(0, nrow = K, ncol = K)
  for(ani in 1:N){
    counts = expected[ani,]
    sum_c = sum(counts)
    et = (sum_c - 1) + exp(-sum_c)
    for(cam in 1:K){
      for(cam2 in 1:K){
        trans_exp[cam,cam2] = trans_exp[cam,cam2] + (et * expected[ani, cam] * expected[ani, cam2]) / sum_c^2
      }
    }
  }

  return(trans_exp)

}

#' FT Transition Metric Preparation
#'
#' For SESCR only. Transition matrix metric discussed in Appendix A
#' Can serve as the `prep` input to `runPPC()`.
#'
#' @inheritParams timeRescaling
#' @export

transMetric_prep_SESCR <- function(S, C, duration = NULL, pars,
                                 hazard, capthist){

  # Constants
  N <- nrow(S)
  K <- nrow(C)

  # Unpack pars
  lambda0 <- pars[1]
  sigma <- pars[2]
  beta <- pars[3]
  r <- pars[4]

  # Get expected counts by ind and trap
  expected <- Integrals(S, C, duration, pars, hazard, N, capthist)

  trans_ind <- list()
  for(ani in 1:N){
    trans_ind[[ani]] = trans_matrix(matrix(S[ani,], nrow=1), C, lambda0, beta, r, sigma)
  }
  # Transition matrix between cameras
  trans_exp = matrix(0, nrow = K, ncol = K)
  for(ani in 1:N){
    counts = expected[ani,]
    sum_c = sum(counts)
    et = (sum_c - 1) + exp(-sum_c)
    for(cam in 1:K){
      for(cam2 in 1:K){
        trans_exp[cam,cam2] = trans_exp[cam,cam2] + (counts[cam] * et * trans_ind[[ani]][cam,cam2] / sum_c)
      }
    }
  }

  return(trans_exp)

}

#' FT Transition Metric Preparation
#'
#' As discussed in Appendix A. Can serve as the `metric` input to `runPPC()`.
#'
#' `N` and duration` are not required for this metric
#' @inheritParams chainMetrics
#' @export

transMetric <- function(data, capthist, N = NULL, K, duration = NULL){

  cameras <- capthist[,3]
  ids <- capthist[,2]
  I <- length(cameras)
  m <- max(ids)

  counts <- matrix(0, nrow = K, ncol = K)

  for(ani in 1:m){
    cameras_ani = cameras[ids == ani]
    Ij = length(cameras_ani)
    if (Ij > 1){
      for(i in 2:Ij){
        counts[cameras_ani[i],cameras_ani[(i-1)]] = counts[cameras_ani[i],cameras_ani[(i-1)]] + 1
      }
    }
  }

  output <- FT_Stat(counts, data)
  return(output)

}
