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
#' @param eps Small constant to stabilize normalization.
#' @param countslist Optional list of per-point weights per time.
#' @return list of TT responsibility matrices (nt x numclust).
#'
#' @export
Estep_flowcut <- function(mn, sigma, prob,
                          ylist,
                          ylist_cens,
                          ylist_uncens,
                          left_cens_list, 
                          right_cens_list,
                          idx_cens_list,
                          cens_lim_l_vec,
                          cens_lim_u_vec,
                          numclust,
                          eps = 1E-20,
                          countslist = NULL) {
  
  TT     <- length(ylist)
  ntlist <- sapply(ylist, nrow)
  dimdat <- dim(mn)[2]
  
  assertthat::assert_that(dim(mn)[1] == length(ylist))
  assertthat::assert_that(length(left_cens_list)  == length(ylist_cens))
  assertthat::assert_that(length(right_cens_list) == length(ylist_cens))
  assertthat::assert_that(length(cens_lim_l_vec) == dimdat)
  assertthat::assert_that(length(cens_lim_u_vec) == dimdat)

  calculate_dens <- function(iclust, tt, y,
                             y_uncens,
                             y_cens,
                             mn, sigma,
                             left_cens_list_t, 
                             right_cens_list_t,
                             integ_lim_l_t,
                             integ_lim_u_t,
                             cens_mat,
                             idx_cens_list_t,
                             cens_lim_l_vec,
                             cens_lim_u_vec) {
    ## Setup
    ##browser()
    mu   <- mn[tt, , iclust]
    nt   <- nrow(y)
    dens <- numeric(nt)
    Sigma<- sigma[iclust, , ]
    
    ##Main logic
    ##for non censored
    if(nrow(y_uncens)==0){
      densvec_uncensored = NULL
    }else{
      densvec_uncensored <- dmvnorm_arma_fast(y_uncens,mu,Sigma)
    }
    
    nt = nrow(y_cens)
    if(nt == 0){
      densvec_censored =NULL
    }else{
      densvec_censored <- sapply(1:nt, function(ii){
        #browser()
        #print(ii)
        uncensored_index <- which(cens_mat[ii,]==0)
        res_cond          <- cond_mean_var_func(y_cens[ii,], mu, Sigma, uncensored_index)
        mu_conditional    <- res_cond$mu_conditional
        mu_observed       <- mu[uncensored_index]
        Sigma_conditional <- res_cond$Sigma_conditional
        Sigma_observed    <- Sigma[uncensored_index,uncensored_index]
        y_observed        <- y_cens[ii, uncensored_index, drop = TRUE]
        
        return(dmvnorm_arma_fast(matrix(y_observed, 1, length(y_observed)),
                                 mu_observed, as.matrix(Sigma_observed), FALSE) *
                 my_pmvnorm(integ_lim_l_t[ii,-uncensored_index],
                            integ_lim_u_t[ii,-uncensored_index],
                            mean = as.vector(mu_conditional),
                            sigma = Sigma_conditional)[1])
      })
      
    }

    dens <- numeric(nrow(y))
    
    if(length(idx_cens_list_t)>0){
      dens[idx_cens_list_t] = densvec_censored
      dens[-idx_cens_list_t] = densvec_uncensored
    }else{
      dens = densvec_uncensored
    }
    
    return(dens)
  }
  
  
  temp_fn_right <- function(M) {
    if (nrow(M) == 0){
      prod_u = rep(-Inf, dimdat)
    }else{
      
      u <- matrix(cens_lim_u_vec, nrow = nrow(M),
                  ncol = length(cens_lim_u_vec), byrow = TRUE)
      
      prod_u <- M * u
      prod_u[is.nan(prod_u)] <- -Inf
      
    } 
    
    prod_u
  }
  
  temp_fn_left <- function(M) {
    if (nrow(M) == 0){
      prod_l = rep(Inf,dimdat)
    }else{
      l <- matrix(cens_lim_l_vec, nrow = nrow(M),
                  ncol = length(cens_lim_l_vec), byrow = TRUE)
      
      prod_l <- M * l
      prod_l[is.nan(prod_l)] <- Inf
    }
   
    prod_l
  }
  
  integration_upper_limits <- lapply(1:length(left_cens_list), 
                                     function(tt){temp_fn_left(left_cens_list[[tt]])})
  
  integration_lower_limits <- lapply(1:length(right_cens_list), 
                                     function(tt){temp_fn_right(right_cens_list[[tt]])})
  
  cens_list <- lapply(1:length(right_cens_list), function(tt){
    cens_tt = left_cens_list[[tt]] + right_cens_list[[tt]] - left_cens_list[[tt]] * right_cens_list[[tt]]
  })
  
  ncol.prob <- ncol(prob)
  resp <- lapply(seq_len(TT), function(tt) {
    #if(tt == 13) browser()
    #print(tt)
    ylist_tt<- ylist[[tt]]
    y_cens_t <- ylist_cens[[tt]]
    y_uncens_t <- ylist_uncens[[tt]]
    left_cens_list_t  <- left_cens_list[[tt]]
    right_cens_list_t <- right_cens_list[[tt]]
    cens_list_t <- cens_list[[tt]]
    integ_lim_l_t <- integration_lower_limits[[tt]]
    integ_lim_u_t <- integration_upper_limits[[tt]]
    idx_cens_list_t = idx_cens_list[[tt]]
    
    #ntt <- ntlist[tt]
    
    if (nrow(ylist_tt) == 0) return(ylist_tt)
    
    #debug(calculate_dens)
    densmat <- sapply(seq_len(numclust),
                      calculate_dens,
                      tt, ylist_tt,
                      y_uncens_t,
                      y_cens_t,
                      mn, sigma,
                      left_cens_list_t, 
                      right_cens_list_t,
                      integ_lim_l_t,
                      integ_lim_u_t,
                      cens_list_t,
                      idx_cens_list_t,
                      cens_lim_l_vec,
                      cens_lim_u_vec)
  
    
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
