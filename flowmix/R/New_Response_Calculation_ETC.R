# New Response Calculation ETC
# Author: Riddhiman Bhattacharya
# Date: 2025-10-14

# --- Preamble ---
library(mvtnorm)
library(Rcpp)
library(assertthat)

# --- New Response Generation ---

#' Computing the new responses conditioned on the observed censoring.
#' @param y response list containing the censored responses. List of length TT; each element is an nt x dim matrix.
#' @param X the matrix of covariates with dimensions TT x d.
#' @param censor_indicator_left list (length TT) of nt x dim matrices indicating left censoring (1) or NA.
#' @param censor_indicator_right list (length TT) of nt x dim matrices indicating right censoring (1) or NA.
#' @param cens_lim_l_vec Lower censoring limits (length = dim).
#' @param cens_lim_u_vec Upper censoring limits (length = dim).
#' @param numclust number of clusters.
#' @param mu_list list of means by cluster; each is a TT x d matrix.
#' @param sigma_list list of covariance matrices by cluster; each is d x d.
#' @return a nested list: new_responses, means, second_moments.
Estep_y <- function(y, X, censor_indicator_left, 
                    censor_indicator_right, cens_lim_l_vec, 
                    cens_lim_u_vec, numclust, mu_list, sigma_list) {
  
  # Basic Checks
  stopifnot(is.list(y))
  stopifnot(is.list(censor_indicator_left))
  stopifnot(is.list(censor_indicator_right))
  stopifnot(is.list(mu_list))
  stopifnot(is.list(sigma_list))
  stopifnot(length(cens_lim_l_vec) == length(cens_lim_u_vec))
  
  # Main Body
  TT     <- length(y)
  dimdat <- ncol(y[[1]])
  ntlist <- sapply(y, nrow)
  
  out_list <- list(new_responses = list(), means = list(), second_moments = list())
  
  res_list <- lapply(seq_len(TT), function(tt) {
    lapply(seq_len(ntlist[tt]), function(ii) {
      res <- cens_cond_normal(ii, tt, y, X,
                              censor_indicator_left, censor_indicator_right,
                              cens_lim_l_vec, cens_lim_u_vec,
                              numclust, mu_list, sigma_list)
      list(new_responses = res$new_response_list,
           mean          = res$all_conditional_means,
           second_moment = res$all_conditional_second_moment_list)
    })
  })
  
  censored_y_out            <- lapply(res_list, function(tt) lapply(tt, `[[`, "new_responses"))
  censored_means_out        <- lapply(res_list, function(tt) lapply(tt, `[[`, "mean"))
  censored_second_moment_out<- lapply(res_list, function(tt) lapply(tt, `[[`, "second_moment"))
  
  out_list$new_responses <- lapply(seq_len(numclust), function(iclust) {
    lapply(seq_len(TT), function(tt) {
      pieces <- lapply(censored_y_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      sapply(pieces, cbind)
    })
  })
  
  out_list$means <- lapply(seq_len(numclust), function(iclust) {
    lapply(seq_len(TT), function(tt) {
      pieces <- lapply(censored_means_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      sapply(pieces, cbind)
    })
  })
  
  out_list$second_moments <- lapply(seq_len(numclust), function(iclust) {
    lapply(seq_len(TT), function(tt) {
      pieces <- lapply(censored_second_moment_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      simplify2array(pieces)
    })
  })
  
  new_response_list <- lapply(out_list$new_responses, function(iclust_block) {
    lapply(iclust_block, function(tt_block) {
      if (is.list(tt_block)) {
        t(do.call(rbind, tt_block))
      } else if (is.matrix(tt_block)) {
        t(tt_block)
      } else {
        matrix(tt_block, nrow = 1)
      }
    })
  })
  
  list(new_responses = new_response_list,
       means         = out_list$means,
       second_moments= out_list$second_moments)
}

# --- Mixture responsibilities (Estep) utilities & call scaffold ---

# The following relies on an external C++ file providing `dmvnorm_arma_fast`.
# Adjust the path to your local file if needed:
# sourceCpp('C://Users//16128//Documents//Sangwon_Project//Codes//dmvnorm.cpp')

#' Calculates responsibilities (posterior membership probabilities).
#' @param mn Array TT x dimdat x numclust of means.
#' @param sigma Array numclust x dimdat x dimdat of covariances.
#' @param prob Matrix TT x numclust of component weights.
#' @param ylist List length TT; each element nt x dimdat data matrix.
#' @param censor_indicator_left List length TT; each element nt x dimdat.
#' @param censor_indicator_right List length TT; each element nt x dimdat.
#' @param cens_lim_l_vec lower limits (length dimdat).
#' @param cens_lim_u_vec upper limits (length dimdat).
#' @param numclust Number of clusters.
#' @param denslist_by_clust Optional precomputed densities.
#' @param first_iter Logical, if TRUE compute densities fresh.
#' @param eps Small constant to stabilize normalization.
#' @param countslist Optional list of per-point weights per time.
#' @return list of TT responsibility matrices (nt x numclust).
Estep <- function(mn, sigma, prob, ylist = NULL,
                  censor_indicator_left,
                  censor_indicator_right,
                  cens_lim_l_vec,
                  cens_lim_u_vec,
                  numclust,
                  denslist_by_clust = NULL,
                  first_iter = FALSE,
                  eps = 1E-20,
                  countslist = NULL) {
  
  TT     <- length(ylist)
  ntlist <- sapply(ylist, nrow)
  dimdat <- dim(mn)[2]
  
  assertthat::assert_that(dim(mn)[1] == length(ylist))
  assertthat::assert_that(length(censor_indicator_left)  == length(ylist))
  assertthat::assert_that(length(censor_indicator_right) == length(ylist))
  assertthat::assert_that(length(cens_lim_l_vec) == dimdat)
  assertthat::assert_that(length(cens_lim_u_vec) == dimdat)
  
  calculate_dens <- function(iclust, tt, y,
                             mn, sigma,
                             censor_indicator_left_t,
                             censor_indicator_right_t,
                             cens_lim_l_vec,
                             cens_lim_u_vec,
                             denslist_by_clust,
                             first_iter) {
    
    mu   <- mn[tt, , iclust]
    nt   <- nrow(y)
    dens <- numeric(nt)
    Sigma<- sigma[iclust, , ]
    
    if (first_iter) {
      for (ii in seq_len(nt)) {
        if (dimdat == 1) {
          dens[ii] <- stats::dnorm(y[ii,], mu, sd = sqrt(sigma[iclust,,])) *
            ((censor_indicator_left_t[ii,]  != 1 | is.na(censor_indicator_left_t[ii,])) &
               (censor_indicator_right_t[ii,] != 1 | is.na(censor_indicator_right_t[ii,]))) +
            pnorm(cens_lim_l_vec, mu, sd = sqrt(sigma[iclust,,])) * (censor_indicator_left_t[ii,] == 1) +
            (1 - pnorm(cens_lim_u_vec, mu, sd = sqrt(sigma[iclust,,]))) * (censor_indicator_right_t[ii,] == 1)
          # NOTE: uses sd (sqrt of variance)
        } else {
          left_cens_index  <- which((censor_indicator_left_t[ii,]  == 1) & (1 - is.na(censor_indicator_left_t[ii,])))
          right_cens_index <- which((censor_indicator_right_t[ii,] == 1) & (1 - is.na(censor_indicator_right_t[ii,])))
          uncensored_index <- which((censor_indicator_left_t[ii,]  != 1 | is.na(censor_indicator_left_t[ii,])) &
                                      (censor_indicator_right_t[ii,] != 1 | is.na(censor_indicator_right_t[ii,])))
          lower_limits <- cens_lim_l_vec[left_cens_index]
          upper_limits <- cens_lim_u_vec[right_cens_index]
          
          res_cond          <- cond_mean_var_func(y[ii,], mu, Sigma, uncensored_index)
          mu_conditional    <- res_cond$mu_conditional
          mu_observed       <- mu[uncensored_index]
          Sigma_conditional <- res_cond$Sigma_conditional
          Sigma_observed    <- Sigma[uncensored_index, uncensored_index, drop = FALSE]
          y_observed        <- y[ii, uncensored_index, drop = TRUE]
          
          if ((length(upper_limits) > 0) || (length(lower_limits) > 0)) {
            p_lower_limit <- numeric(dimdat)
            p_upper_limit <- numeric(dimdat)
            p_lower_limit[left_cens_index]  <- -Inf
            p_lower_limit[uncensored_index] <- -Inf
            p_lower_limit[right_cens_index] <- upper_limits
            p_upper_limit[left_cens_index]  <- lower_limits
            p_upper_limit[uncensored_index] <- Inf
            p_upper_limit[right_cens_index] <- Inf
            
            dens[ii] <-
              dmvnorm_arma_fast(matrix(y_observed, 1, length(y_observed)),
                                mu_observed, as.matrix(Sigma_observed), FALSE) *
              my_pmvnorm(p_lower_limit[sort(c(left_cens_index, right_cens_index))],
                         p_upper_limit[sort(c(left_cens_index, right_cens_index))],
                         mean = as.vector(mu_conditional),
                         sigma = Sigma_conditional)[1]
          } else {
            dens[ii] <- dmvnorm_arma_fast(y[ii, , drop = FALSE],
                                          mu, as.matrix(sigma[iclust, , ]), FALSE)
          }
        }
      }
    } else {
      dens <- unlist(denslist_by_clust[[iclust]][[tt]])
    }
    dens
  }
  
  ncol.prob <- ncol(prob)
  resp <- lapply(seq_len(TT), function(tt) {
    ylist_tt               <- ylist[[tt]]
    censor_indicator_left_t  <- censor_indicator_left[[tt]]
    censor_indicator_right_t <- censor_indicator_right[[tt]]
    ntt <- ntlist[tt]
    
    if (nrow(ylist_tt) == 0) return(ylist_tt)
    
    densmat <- sapply(seq_len(numclust),
                      calculate_dens,
                      tt, ylist_tt, mn, sigma,
                      censor_indicator_left_t,
                      censor_indicator_right_t,
                      cens_lim_l_vec,
                      cens_lim_u_vec,
                      denslist_by_clust, first_iter)
    
    wt.densmat <- matrix(prob[tt,], nrow = ntlist[tt], ncol = ncol.prob, byrow = TRUE) * densmat
    wt.densmat <- wt.densmat + eps
    wt.densmat <- wt.densmat / rowSums(wt.densmat)
    
    if (!is.null(countslist)) {
      wt.densmat <- wt.densmat * countslist[[tt]]
    }
    wt.densmat
  })
  
  resp
}

# --- Example scaffolding (adjust or comment out as needed) ---

# Example responsibility/prior matrices shape setup
# numclust <- 2
# TT <- 30
# pmat <- matrix(NA, TT, numclust)
# for (iclust in 1:numclust) {
#   pmat[1:(TT/2 + iclust), iclust] <- 1/2
#   pmat[(TT/2 + iclust + 1):TT, iclust] <- 1/3
# }
# mn_array    <- simplify2array(mn_list)                 # requires mn_list
# sigma_array <- aperm(simplify2array(Sig_list), c(3,1,2))  # requires Sig_list

# Example: compile external C++ density (edit path)
# sourceCpp('C://Users//16128//Documents//Sangwon_Project//Codes//dmvnorm.cpp')

# Example call (ensure objects exist: mn_array, sigma_array, pmat, ylist, cens_ind_left/right, cens_vec_1/2)
# Estep(mn_array, sigma_array, pmat, ylist,
#       censor_indicator_left = cens_ind_left,
#       censor_indicator_right = cens_ind_right,
#       cens_lim_l_vec = cens_vec_1,
#       cens_lim_u_vec = cens_vec_2,
#       numclust = 2,
#       denslist_by_clust = NULL,
#       first_iter = TRUE,
#       eps = 1e-20,
#       countslist = NULL)
