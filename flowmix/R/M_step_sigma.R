# M_step_Sigma
# Author: Riddhiman Bhattacharya
# Date: 2025-10-22

# --- Libraries ---
library(abind)

# --- Sigma_hat_calculation_fn -------------------------------------------------
#' Function to assist in calculating Sigma_hat for the censored terms
#' @param iclust cluster index
#' @param resp responsibilities list (length TT) of nt x K matrices
#' @param ylist list of length K; each [[iclust]] is a list length TT of nt x d matrices (censored y)
#' @param second_moment_list list (length K) -> [[iclust]][[tt]][[ii]] = d x d 2nd moment (cond.)
#' @param mn fitted means array (TT x d x K)
#' @param cens_mat matrix with columns: (tt, ii, d censor flags) where 1 indicates censored
#' @return d x d matrix for censored residual contribution divided by sum of resp for cluster
Sigma_hat_calculation_fn <- function(iclust, resp, ylist, second_moment_list, mn, cens_mat){
  ylist <- ylist[[iclust]]
  TT <- length(ylist)
  ntlist <- sapply(ylist, nrow)
  dimdat <- ncol(ylist[[1]])
  numclust <- ncol(resp[[1]])
  
  mat_calc_func <- function(iclust, tt, ii, resp, ylist, second_moment_list, cens_mat){
    mat_out <- matrix(0, dimdat, dimdat)
    condition <- ((cens_mat[,1] == tt) & (cens_mat[,2] == ii))
    cens_vec <- cens_mat[condition, , drop = FALSE]
    if (nrow(cens_vec) == 0) return(mat_out)
    if (sum(cens_vec[1, 3:(dimdat+2)]) == 0) {
      return(mat_out)
    } else {
      indices_temp <- which(cens_vec[1, 3:(dimdat+2)] == 1)
      mean_temp <- ylist[[tt]][ii, indices_temp]
      mat_out[indices_temp, indices_temp] <-
        resp[[tt]][ii, iclust] * (second_moment_list[[iclust]][[tt]][[ii]] - mean_temp %*% t(mean_temp))
      return(mat_out)
    }
  }
  
  out <- lapply(seq_len(TT), function(tt) {
    lapply(seq_len(ntlist[tt]), function(ii) {
      mat_calc_func(iclust, tt, ii, resp, ylist, second_moment_list, cens_mat)
    })
  })
  
  resp.long <- do.call(rbind, resp)
  resp.sum  <- apply(resp.long, 2, sum)[iclust]
  
  flatten_out <- unlist(out, recursive = FALSE)
  out_mat <- Reduce(`+`, flatten_out)
  out_mat / resp.sum
}

# --- Mstep_sigma --------------------------------------------------------------
#' M-step for covariance matrix Sigma_k for clusters k = 1..K
#'
#' @param resp responsibilities (list length TT of nt x K)
#' @param ylist list length K; [[iclust]] is list length TT of nt x d (y_new)
#' @param cens_mat censoring matrix with columns (tt, ii, flags...)
#' @param first_moment_list list of first conditional moments (by cluster/time/obs)
#' @param second_moment_list list of second conditional moments (by cluster/time/obs)
#' @param mn fitted means array (TT x d x K)
#' @return (K x d x d) array of covariance estimates
Mstep_sigma <- function(resp, ylist, cens_mat, first_moment_list, second_moment_list, mn){
  TT      <- length(ylist[[1]])
  ntlist  <- sapply(ylist[[1]], nrow)
  dimdat  <- ncol(ylist[[1]][[1]])
  numclust<- ncol(resp[[1]])
  
  # Map row indices of mn to long form
  irows <- rep(seq_len(nrow(mn)), times = ntlist)
  
  vars <- vector("list", numclust)
  
  for (iclust in 1:numclust) {
    ylong <- do.call(rbind, ylist[[iclust]])
    
    resp.thisclust <- lapply(resp, function(myresp) myresp[, iclust, drop = TRUE])
    resp.long <- do.call(c, resp.thisclust)
    
    mnlong <- mn[irows, , iclust, drop = FALSE]
    if (is.vector(mnlong)) mnlong <- cbind(mnlong)
    
    # estepC is assumed to be provided by your C++ code; it returns d x d matrix:
    vars[[iclust]] <-
      estepC(ylong, mnlong, sqrt(resp.long), sum(resp.long)) +
      Sigma_hat_calculation_fn(iclust, resp, ylist, second_moment_list, mn, cens_mat)
  }
  
  sigma_array <- array(NA, dim = c(numclust, dimdat, dimdat))
  for (iclust in 1:numclust) sigma_array[iclust, , ] <- vars[[iclust]]
  
  stopifnot(all(dim(sigma_array) == c(numclust, dimdat, dimdat)))
  sigma_array
}

# --- Matrix inverse square-root helper ---------------------------------------
#' Given a PSD matrix a, compute a^{-1/2}
mtsqrt_inv <- function(a){
  a.eig <- eigen(a)
  vec <- 1 / sqrt(a.eig$values)
  mat <- if (length(vec) == 1) vec else diag(vec)
  a.eig$vectors %*% mat %*% t(a.eig$vectors)
}

# --- Build censoring matrix from left/right indicators -----------------------
# (Assumes objects: cens_ind_left, cens_ind_right exist in the environment)
# Produces: cens_mat and cens_ind
build_cens_mat <- function(cens_ind_left, cens_ind_right){
  TT <- length(cens_ind_left)
  left_cens_index_temp  <- vector("list", TT)
  right_cens_index_temp <- vector("list", TT)
  
  for (tt in 1:TT) {
    left_cens_index_temp[[tt]]  <- cbind(rep(tt, nrow(cens_ind_left[[tt]])),
                                         seq_len(nrow(cens_ind_left[[tt]])),
                                         cens_ind_left[[tt]])
    right_cens_index_temp[[tt]] <- cbind(rep(tt, nrow(cens_ind_right[[tt]])),
                                         seq_len(nrow(cens_ind_right[[tt]])),
                                         cens_ind_right[[tt]])
  }
  
  temp_mat_left  <- do.call(rbind, left_cens_index_temp)
  temp_mat_right <- do.call(rbind, right_cens_index_temp)
  
  temp_mat_left[is.na(temp_mat_left)]   <- 0
  temp_mat_right[is.na(temp_mat_right)] <- 0
  
  cens_mat <- pmax(temp_mat_left, temp_mat_right)
  cens_ind <- 1 * (apply(cens_mat[, -c(1,2), drop = FALSE], 1, sum) != 0)
  list(cens_mat = cens_mat, cens_ind = cens_ind)
}

# --- Example wiring (comment/uncomment as needed) ----------------------------
# The following expects that you have already created:
#   - estepy_saved (from Estep_y)
#   - y_new (from Estep_y$new_responses)
#   - resp  (from Estep)
#   - mn_array, sigma_array, mn_list / Sig_list
#   - cens_ind_left, cens_ind_right
#
# Example:
# first_moment_list  <- estepy_saved$means
# second_moment_list <- estepy_saved$second_moments
#
# tmp <- build_cens_mat(cens_ind_left, cens_ind_right)
# cens_mat <- tmp$cens_mat
#
# # Compute Sigma (requires estepC from your C++ code)
# Mstep_sigma(resp, y_new, cens_mat, first_moment_list, second_moment_list, mn_array)
