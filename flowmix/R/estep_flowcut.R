##' posterieri membership probabilities).
##'
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
#'
#' @export
Estep_flowcut <- function(mn, sigma, prob, ylist = NULL,
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
    #browser()
    mu   <- mn[tt, , iclust]
    nt   <- nrow(y)
    dens <- numeric(nt)
    Sigma<- sigma[iclust, , ]
    
    if (first_iter) {
      dens <- sapply(seq_len(nt), function(ii) {
        if (dimdat == 1) {
          stats::dnorm(y[ii,], mu, sd = sqrt(sigma[iclust,,])) *
            ((censor_indicator_left_t[ii,]  != 1 | is.na(censor_indicator_left_t[ii,])) &
               (censor_indicator_right_t[ii,] != 1 | is.na(censor_indicator_right_t[ii,]))) +
            pnorm(cens_lim_l_vec, mu, sd = sqrt(sigma[iclust,,])) * (censor_indicator_left_t[ii,] == 1) +
            (1 - pnorm(cens_lim_u_vec, mu, sd = sqrt(sigma[iclust,,]))) * (censor_indicator_right_t[ii,] == 1)
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
            
            dmvnorm_arma_fast(matrix(y_observed, 1, length(y_observed)),
                              mu_observed, as.matrix(Sigma_observed), FALSE) *
              my_pmvnorm(p_lower_limit[sort(c(left_cens_index, right_cens_index))],
                         p_upper_limit[sort(c(left_cens_index, right_cens_index))],
                         mean = as.vector(mu_conditional),
                         sigma = Sigma_conditional)[1]
          } else {
            dmvnorm_arma_fast(y[ii, , drop = FALSE],
                              mu, as.matrix(sigma[iclust, , ]), FALSE)
          }
        }
      })
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


# --- New Response Generation ---

#' Computing the new responses conditioned on the observed censoring.
#'@param y response list containing the censored responses. List containing TT elements with each element as a matrix of size ntxdim.
#'@param X the matrix of covariates with dimensions tt x d.
#'@param censor_indicator_left A list checking whether an observation is left censored. The list if of size TT with each entry being a ntxdim  matrix.
#'@param censor_indicator_right A list checking whether an observation is right censored. The list if of size TT with each entry being a ntxdim  matrix.
#'@param cens_lim_l_vec Lower censoring limits.
#'@param cens_lim_u_vec Upper censoring limits.
#'@param numclust number of clusters
#'@param mu_list list of means by cluster. Each cluster has a matrix of dim TT x d.
#'@param sigma_list list of covariance matrices by cluster. Each cluster has a matrix of dim d x d.
#'
#'@return a nested list containing the responses per time point,observation,cluster in that order of nesting.
#'
#' @export

Estep_y_flowcut <- function(y, X, censor_indicator_left, 
                    censor_indicator_right, cens_lim_l_vec, 
                    cens_lim_u_vec, numclust, mu_list,sigma_list){
  #browser()
  ##Basic Checks
  stopifnot(is.list(y))
  stopifnot(is.list(censor_indicator_left))
  stopifnot(is.list(censor_indicator_right))
  stopifnot(is.list(mu_list))
  stopifnot(is.list(sigma_list))
  stopifnot(length(cens_lim_l_vec)==length(cens_lim_u_vec))
  
  ##Main Body
  TT = length(y)
  dimdat = ncol(y[[1]])
  ntlist = sapply(y,nrow)
  out_list=list()
  out_list$new_responses = list()
  out_list$means = list()
  out_list$second_moments = list()
  
  
  #debug(cens_cond_normal)  
  
  res_list <- lapply(seq_len(TT), function(tt) {
    n_tt <- ntlist[tt]
    lapply(seq_len(n_tt), function(ii) {
      #if((debug_global==9)&&(tt==3)&&(ii==71)) debug(cens_cond_normal) 
      res <- cens_cond_normal(
        ii, tt, y, X,
        censor_indicator_left, censor_indicator_right,
        cens_lim_l_vec, cens_lim_u_vec,
        numclust, mu_list, sigma_list
      )
      # directly return res; don't rebuild mini-lists
      list(
        new_responses  = res$new_response_list,
        mean           = res$all_conditional_means,
        second_moment  = res$all_conditional_second_moment_list
      )
    })
  })
  
  censored_y_out <- lapply(res_list, function(tt) lapply(tt, `[[`, "new_responses"))
  censored_means_out <- lapply(res_list, function(tt) lapply(tt, `[[`, "mean"))
  censored_second_moment_out <- lapply(res_list, function(tt) lapply(tt, `[[`, "second_moment"))
  
  out_list$new_responses <- lapply(1:numclust, function(iclust) {
    lapply(1:TT, function(tt) {
      # pick iclust-th component from each ii in this time tt
      pieces <- lapply(censored_y_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      sapply(pieces, cbind)
    })
  })
  
  out_list$means <- lapply(1:numclust, function(iclust) {
    lapply(1:TT, function(tt) {
      # pick iclust-th component from each ii in this time tt
      pieces <- lapply(censored_means_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      sapply(pieces, cbind)
    })
  })
  
  out_list$second_moments <- lapply(1:numclust, function(iclust) {
    lapply(1:TT, function(tt) {
      pieces <- lapply(censored_second_moment_out[[tt]], function(ii_elem) ii_elem[[iclust]])
      simplify2array(pieces)
    })
  })
  
  new_response_list = lapply(out_list$new_responses, function(iclust){
    lapply(iclust,function(tt){
      if (is.list(tt)) {
        t(do.call(rbind, tt))
      } else {
        # If it's already a matrix, return as is or handle appropriately
        if (is.matrix(tt)) {
          t(tt)
        } else {
          # If it's a vector, make it a 1-row matrix
          matrix(tt, nrow = 1)
        }
      }
    })
  })
  
  return(list(new_responses = new_response_list,
              means = out_list$means,
              second_moments = out_list$second_moments))
}

