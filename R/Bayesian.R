#' Run Bayesian PPCs for SCR Models
#'
#' Runs a Bayesian PPC after MCMC has been run and a `metric` indicated
#'
#' Regardless of `mode`, the metric is computed on both simulated and the actual `capthist` \cr
#' When `mode=0`, the metric depends only on a `capthist` object. \cr
#' When `mode=1` or `mode=2`, the metric depends on `capthist` and the output of `prep` \cr
#' When `mode=2`, `prep` requires a further simulated `capthist`
#'
#' @inheritParams timeRescaling
#' @param S matrix of activity center estimates. X
#' @param samples matrix of parameter estimates
#' @param Z matrix of data augmentation binary variables
#' @param simulate function that simulates model and returns a `capthist` object
#' @param metric function that returns one or more scalar test statistics
#' @param prep auxillary function needed for modes 1 and 2
#' @param mode Type of PPC metric, see details
#' @param n_metrics number of scalar metrics to compute
#' @export

runPPC <- function(S, C, duration, samples, Z, hazard, simulate, metric, prep,
                   mode = 0, n_metrics = 1, capthist){

  times <- capthist[,1]
  ids <- capthist[,2]
  cameras <- capthist[,3]

  # Constants
  I <- nrow(samples)
  m <- max(ids)
  K <- nrow(C)


  # pre-concatenate Z if missing the observed individuals
  M = ncol(S)/2
  if (M %% 2 != 0){
    stop("S must have 2 columns for every member of the superpopulation")
  }

  if (ncol(Z) == (M - m)){
    Z = cbind(matrix(1,nrow=I,ncol=m), Z)
  } else if (ncol(Z) != I){
    stop("Z must either have 2M or 2(M-m) columns")
  }


  # Setup storage for metrics
  metrics_sim <- metrics_data <- matrix(0, nrow = I, ncol = n_metrics)

  info_sum <- data_sum <- numeric(I)

  if (mode == 0){
    metrics_data[1:I, ] <- metric(data = NULL, capthist)
  }

  for(i in 1:I){
    # Get pars
    S_i <- matrix(S[i, ], nrow = M, ncol = 2)
    S_i <- S_i[which(Z[i,] == 1),]
    N_i <- sum(Z[i,])
    pars_i <- samples[i,]

    capthist_sim <- simulate(S_i, C, duration, pars_i)

    if (mode == 0){

      metrics_sim[i,] <- metric(data = NULL, capthist_sim, N = NULL, K = NULL,
                                 duration = NULL)
    } else if (mode == 1){
      info <- prep(S_i, C, duration, pars_i,
                  hazard, capthist = NULL)
      info_sum[i] <- sum(info)
      data_sum[i] <- nrow(capthist_sim)
      metrics_sim[i, ] <- metric(info, capthist_sim, N_i, K, duration)
      metrics_data[i, ] <- metric(info, capthist, N_i, K, duration)
    } else if (mode == 2){
      capthist_ref <- simulate(S_i, C, duration, pars_i)
      info <- prep(S_i, C, duration, pars_i,
                   hazard, capthist = capthist_ref)
      info_sum[i] <- sum(info)
      data_sum[i] <- nrow(capthist_ref)
      metrics_sim[i, ] <- metric(info, capthist_sim, N_i, K, duration)
      metrics_data[i, ] <- metric(info, capthist, N_i, K, duration)

    }

    if (i %% 100 == 0){
      print(i)
    }
  }

  BPV <- numeric(n_metrics)
  for(j in 1:n_metrics){
    BPV[j] <- mean(metrics_sim[,j] > metrics_data[,j])
    # Ties count as half above, half below
    BPV[j] <- BPV[j] + 0.5*mean(metrics_sim[,j] == metrics_data[,j])
  }

  print(mean(info_sum)); print(mean(data_sum))
  print(summary(info_sum)); print(summary(data_sum))

  return(list(P_values = BPV, metrics_sim = metrics_sim, metrics_data = metrics_data))
}
