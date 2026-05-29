#' Documentation needed.
#'
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
    ## Setup
    mu   <- mn[tt, , iclust]
    nt   <- nrow(y)
    dens <- numeric(nt)
    Sigma<- sigma[iclust, , ]


    ## Main loop
    if (first_iter) {
      dens <- sapply(seq_len(nt), function(ii) {
        if (dimdat == 1) {
          stats::dnorm(y[ii,], mu, sd = sqrt(sigma[iclust,,])) *
            ((censor_indicator_left_t[ii,]  != 1 | is.na(censor_indicator_left_t[ii,])) &
               (censor_indicator_right_t[ii,] != 1 | is.na(censor_indicator_right_t[ii,]))) +
            pnorm(cens_lim_l_vec, mu, sd = sqrt(sigma[iclust,,])) * ((censor_indicator_left_t[ii,] == 1)&(1 - is.na(censor_indicator_left_t[ii,]))) +
            (1 - pnorm(cens_lim_u_vec, mu, sd = sqrt(sigma[iclust,,]))) * ((censor_indicator_right_t[ii,] == 1)&(1 - is.na(censor_indicator_right_t[ii,])))
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
#'@param mu_list list of means by cluster. Each cluster has a matrix of dim TT x d. (TODO: update to mn)
#'@param sigma_list list of covariance matrices by cluster. Each cluster has a matrix of dim d x d. (TODO: update to mn)
#'
#'@return a nested list containing the responses per time point,observation,cluster in that order of nesting.
#'
#' @export

Estep_y_flowcut <- function(ylist_cens,
                            ylist_uncens,
                            X,
                            ## censor_indicator_left,
                            ## censor_indicator_right,
                            left_cens_list,
                            right_cens_list,
                            cens_lim_l_vec,
                            cens_lim_u_vec,
                            numclust,
                            mn,
                            sigma){

  ## Basic Checks
  stopifnot(length(cens_lim_l_vec)==length(cens_lim_u_vec))
  

  ## make censored conditional normal (1) mean and 2nd moment, and (2)
  ## imputations of the response.
  res = cens_cond_normal_new(ylist_cens, X, left_cens_list,
                             right_cens_list, cens_lim_l_vec, cens_lim_u_vec,
                             numclust, mn, sigma)

  ## Return the result
  stopifnot(all(dim(res$new_responses[[1]]) == dim(ylist_cens[[1]])))
  return(res)
}
