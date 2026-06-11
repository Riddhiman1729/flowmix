#' Impute censored responses using a fitted flowcut object.
#'
#' @param flowcutobj A flowcut class object obtained from running the flowcut function.
#' @param datobj A data object containing a list of responses, left and right censor indicators,
#' and, if required, a list of covariates.
#'
#' @return A list of imputed responses and indexes which are censored.
#' @export
imputation_prediction <- function(flowcutobj, datobj) {
  
  stopifnot("ylist" %in% names(datobj))
  stopifnot("censor_indicator_left" %in% names(datobj))
  stopifnot("censor_indicator_right" %in% names(datobj))
  
  ylist <- datobj$ylist
  censor_indicator_left <- datobj$censor_indicator_left
  censor_indicator_right <- datobj$censor_indicator_right
  cens_lim_l_vec <- datobj$cens_lim_l_vec
  cens_lim_u_vec <- datobj$cens_lim_u_vec
  numclust <- datobj$numclust
  
  out_obj_dat <- flowmix:::my_convert(
    ylist,
    censor_indicator_left,
    censor_indicator_right
  )
  
  ylist_cens <- out_obj_dat$ylist_cens
  ylist_uncens <- out_obj_dat$ylist_uncens
  left_cens_list <- out_obj_dat$left_cens_list
  right_cens_list <- out_obj_dat$right_cens_list
  idx_cens_list <- out_obj_dat$idx_cens_list
  idx_uncens_list <- out_obj_dat$idx_uncens_list
  
  TT <- length(ylist)
  ntlist <- sapply(ylist, nrow)
  dimdat <- ncol(ylist[[1]])
  
  sigmaarray <- array(NA_real_, dim = c(numclust, dimdat, dimdat))
  
  for (iclust in seq_len(numclust)) {
    sigmaarray[iclust, , ] <- flowcutobj$sigma[iclust, , ]
  }
  
  pimat <- flowcutobj$prob
  
  mnarray <- array(NA_real_, dim = c(TT, dimdat, numclust))
  
  if ("X" %in% names(datobj)) {
    
    X <- datobj$X
    Xa <- cbind(1, X)
    
    for (iclust in seq_len(numclust)) {
      mnarray[, , iclust] <- Xa %*% flowcutobj$beta[[iclust]]
    }
    
    estep_yy_out <- flowmix::Estep_y_flowcut(
      ylist_cens,
      ylist_uncens,
      X,
      left_cens_list,
      right_cens_list,
      cens_lim_l_vec,
      cens_lim_u_vec,
      numclust,
      mnarray,
      sigmaarray
    )
    
  } else {
    
    for (iclust in seq_len(numclust)) {
      mnarray[, , iclust] <- flowcutobj$mn[, , iclust]
    }
    
    estep_yy_out <- flowtrend::Estep_y_cut(
      ylist_cens,
      ylist_uncens,
      left_cens_list,
      right_cens_list,
      cens_lim_l_vec,
      cens_lim_u_vec,
      numclust,
      mnarray,
      sigmaarray
    )
  }
  
  new_responses <- estep_yy_out$new_responses %>% purrr::list_transpose()
  
  
  ylist_imputed_temp <- vector("list", numclust)
  
  for (iclust in seq_len(numclust)) {
    ylist_imputed_temp[[iclust]] <- flowmix:::my_reassemble(
      ylist_uncens,
      new_responses[[iclust]], 
      idx_cens_list,
      idx_uncens_list
    )
  }
  
  weighted_ylist_imputed_temp <- lapply(seq_len(numclust), function(iclust) {
    lapply(seq_len(TT), function(tt) {
      pimat[tt, iclust] * ylist_imputed_temp[[iclust]][[tt]]
    })
  })
  
  ylist_imputed <- lapply(seq_len(TT), function(tt) {
    Reduce("+", lapply(seq_len(numclust), function(iclust) {
      weighted_ylist_imputed_temp[[iclust]][[tt]]
    }))
  })
  
  return(list(
    ylist_imputed = ylist_imputed,
    imputed_responses = new_responses,
    idx_cens_list = idx_cens_list,
    idx_uncens_list = idx_uncens_list
  ))
}

