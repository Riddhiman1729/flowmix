# Normal Probability_Calculations_ver_3
# Author: Riddhiman Bhattacharya
# Date: 2025-09-29

# ---- Preamble ----
library(mvtnorm)

# ---- Generate data ----
## TODO: Replace the entire censoring structure with a censoring matrix

## Helper for generating mixture of means
## @param tt time point 
## @param nt vector of number of particles at each time  point.
## @param mnlist list of means
## @param pilist list of cluster probabilities by time.
get_mixture_at_timepoint <- function(tt, nt, mnlist, pilist, sigma = NULL, sigmalist = NULL){
  ## Check dimensions
  stopifnot(length(mnlist) == length(pilist))
  if (!is.null(sigma)) stopifnot(nrow(sigma) == ncol(mnlist[[1]]))
  dimdat <- ncol(mnlist[[1]])
  
  ## Draw data from randomly chosen mean.
  prob <- sapply(pilist, function(prob) prob[[tt]])
  prob <- prob / sum(prob)
  stopifnot(sum(prob) == 1)
  mns <- lapply(seq_along(mnlist), function(ii) { mnlist[[ii]][tt,] })
  numclust <- length(prob)
  
  ## Samples nt memberships out of 1:numclust according to the probs in prob.
  draws <- sample(seq_len(numclust), size = nt, replace = TRUE, prob = prob)
  
  ## Randomly chosen means according to pi
  if (!is.null(sigmalist)) {
    dat <- vector("list", numclust)
    for (iclust in seq_len(numclust)) {
      ntk <- sum(draws == iclust)
      if (ntk == 0) next
      membership <- rep(iclust, ntk)
      dat[[iclust]] <- cbind(
        MASS::mvrnorm(n = ntk, mu = mns[[iclust]], Sigma = sigmalist[[iclust]]),
        membership
      )
    }
    datapoints <- do.call(rbind, dat)
  } else {
    ## Add noise to the means.
    means <- mns[draws]
    means <- do.call(rbind, means)
    datapoints <- means + MASS::mvrnorm(n = nt, mu = rep(0, dimdat), Sigma = sigma)
  }
  return(datapoints)
}

## Data generator
## @return list containing ylist, X, and censoring-related outputs
generate_data_generic <- function(p = 5, TT = 50, fac = 1, nt = 1000, dimdat = 2,
                                  cens_lim_l_vec, cens_lim_u_vec,
                                  seed = NULL, prob1 = 3/4){
  if (!is.null(seed)) set.seed(seed)
  stopifnot(dimdat %in% c(2,3))
  if (length(cens_lim_l_vec) != length(cens_lim_u_vec)) stop("Error - censoring lengths are not equal.")
  if (length(cens_lim_l_vec) != dimdat) stop("Error - censoring does not have the correct dimension.")
  
  ## Generate covariates.
  X <- matrix(stats::rnorm(TT * p), ncol = p, nrow = TT)
  X[,1] <- sin((1:TT) / TT * 4 * pi)
  X[,2] <- c(rep(0, TT/2), rep(1, TT/2))
  Xa <- cbind(rep(1, TT), X)
  
  ## Coefficients (dimdat x p) via stacked columns in beta*
  beta11 <- cbind(c(1,0), c(1,0))
  beta00 <- cbind(c(0,0), c(0,0))
  beta10 <- cbind(c(1,0), c(0,0))
  intercept1 <- c(0,0)
  intercept2 <- c(1,1)
  intercept3 <- c(0,5)
  intercept4 <- c(5,5)
  beta1 <- rbind(intercept1, beta11/3)
  beta2 <- rbind(intercept2, beta10/3)
  beta3 <- rbind(intercept3, beta00)
  beta4 <- rbind(intercept4, beta00)
  beta1 <- rbind(beta1, matrix(0, nrow = p-2, ncol = 2))
  beta2 <- rbind(beta2, matrix(0, nrow = p-2, ncol = 2))
  beta3 <- rbind(beta3, matrix(0, nrow = p-2, ncol = 2))
  beta4 <- rbind(beta4, matrix(0, nrow = p-2, ncol = 2))
  
  if (dimdat == 3) {
    beta1 <- cbind(beta1, c(0, rep(0, p)))
    beta2 <- cbind(beta2, c(1, rep(0, p)))
    beta3 <- cbind(beta3, c(0, rep(0, p)))
    beta4 <- cbind(beta4, c(5, rep(0, p)))
  }
  betalist <- list(beta1, beta2, beta3, beta4)
  
  ## Means per component
  mn1 <- Xa %*% beta1
  mn2 <- Xa %*% beta2
  mn3 <- Xa %*% beta3
  mn4 <- Xa %*% beta4
  mnlist <- list(mn1, mn2, mn3, mn4)
  
  ## Mixture weights
  pi1 <- rep(prob1, TT)
  pi2 <- rep((1-prob1)/2, TT)
  pi3 <- rep((1-prob1)/2, TT)
  pi4 <- c(rep(0, TT/2), rep(1/8, TT/2))
  pilist <- list(pi1, pi2, pi3, pi4)
  
  pimat <- do.call(cbind, pilist)
  pimat <- pimat / rowSums(pimat)
  
  ## Particle counts
  ntlist <- c(rep(nt, TT/2), rep(nt*9/8, TT/2))
  ntlist <- apply(ntlist * cbind(pi1, pi2, pi3, pi4), 1, sum)
  
  ## Covariances
  sigma1 <- diag(c(1,1,1))[1:dimdat, 1:dimdat]
  sigma2 <- diag(c(10,1,1))[1:dimdat, 1:dimdat]
  sigma3 <- matrix(c(3,1.5,0, 1.5,3,0, 0,0,1), ncol = 3)[1:dimdat, 1:dimdat]
  sigma4 <- diag(c(1,1,1))[1:dimdat, 1:dimdat]
  sigmalist <- list(sigma1, sigma2, sigma3, sigma4)
  sigmalist <- lapply(sigmalist, function(a) a/3 * fac)
  stopifnot(all(unlist(lapply(sigmalist, dim)) == dimdat))
  
  ## Generate data per time
  datapoints <- lapply(1:TT, function(tt) {
    get_mixture_at_timepoint(tt, ntlist[[tt]], mnlist, pilist, sigmalist = sigmalist)
  })
  
  ## Reformat
  ylist <- lapply(datapoints, cbind)
  classlist <- lapply(ylist, function(a) a[, "membership"])
  `%ni%` <- function(x, y) !(x %in% y)
  ylist <- lapply(ylist, function(a) a[, which(colnames(a) %ni% "membership")])
  
  ylist_cens <- list()
  cens_ind_left <- list()
  cens_ind_right <- list()
  for (tt in 1:TT) {
    cens_ind_left[[tt]]  <- matrix(NA, nrow(ylist[[tt]]), ncol(ylist[[tt]]))
    cens_ind_right[[tt]] <- matrix(NA, nrow(ylist[[tt]]), ncol(ylist[[tt]]))
    ylist_cens[[tt]]     <- matrix(NA, nrow(ylist[[tt]]), ncol(ylist[[tt]]))
    for (coord in 1:dimdat) {
      cens_ind_left[[tt]][, coord]  <- ifelse(ylist[[tt]][, coord] <= cens_lim_l_vec[coord], 1, NA)
      cens_ind_right[[tt]][, coord] <- ifelse(ylist[[tt]][, coord] >= cens_lim_u_vec[coord], 1, NA)
      
      ylist_cens[[tt]][, coord] <-
        cens_lim_l_vec[coord] * ifelse(is.na(cens_ind_left[[tt]][, coord]), 0, 1) +
        cens_lim_u_vec[coord] * ifelse(is.na(cens_ind_right[[tt]][, coord]), 0, 1) +
        ylist[[tt]][, coord] * (
          1 - ifelse(is.na(cens_ind_left[[tt]][, coord]), 0, 1) -
            ifelse(is.na(cens_ind_right[[tt]][, coord]), 0, 1) +
            ifelse(is.na(cens_ind_left[[tt]][, coord]), 0, 1) *
            ifelse(is.na(cens_ind_right[[tt]][, coord]), 0, 1)
        )
    }
  }
  
  ## Return
  list(
    ylist = ylist_cens,
    classlist = classlist,
    censor_indicator_left  = cens_ind_left,
    censor_indicator_right = cens_ind_right,
    X = X,
    Xa = Xa,
    sigmalist = sigmalist,
    betalist = betalist,
    ntlist = ntlist
  )
}

## Example data
set.seed(0)
dimdatt <- 3
cens_vec_1 <- rep(-.1, dimdatt)
cens_vec_2 <- rep(.7,  dimdatt)
datobj <- generate_data_generic(
  p = dimdatt, TT = 30, fac = .2, nt = 200, dimdat = dimdatt,
  cens_lim_l_vec = cens_vec_1, cens_lim_u_vec = cens_vec_2
)
ylist <- datobj$ylist
X <- datobj$X
cens_ind_left  <- datobj$censor_indicator_left
cens_ind_right <- datobj$censor_indicator_right

# ---- Probability helper ----
# Computes the normal probabilities giving a small lower bound to very low probability events
my_pmvnorm <- function(lower, upper, mean, sigma){
  prob <- tryCatch({
    mvtnorm::pmvnorm(lower, upper, mean = mean, sigma = sigma)
  }, error = function(e) {
    NaN
  })
  if (is.nan(prob) || prob < 1e-12) 1e-12 else mvtnorm::pmvnorm(lower, upper, mean = mean, sigma = sigma)
}

# ---- Conditional mean/variance of MVN given observed dims ----
cond_mean_var_func <- function(y, mu, Sigma, observed_dims){
  dimdat <- length(y)
  stopifnot(length(observed_dims) == 0 || all(observed_dims %in% seq_len(dimdat)))
  stopifnot(length(mu) == dimdat)
  stopifnot(all(dim(Sigma) == dimdat))
  
  y_observed_dims <- y[observed_dims]
  muu <- mu[-observed_dims]
  muo <- mu[observed_dims]
  Sigmaoo <- Sigma[observed_dims, observed_dims]
  Sigmauu <- Sigma[-observed_dims, -observed_dims]
  Sigmauo <- Sigma[-observed_dims, observed_dims]
  Sigmaou <- Sigma[observed_dims, -observed_dims]
  
  if (any(y_observed_dims == -Inf) || any(y_observed_dims == Inf)) {
    mu_out <- muu
    Sigma_out <- Sigmauu
  } else {
    if ((length(observed_dims) == length(mu)) || (length(observed_dims) == 0)) {
      mu_out <- mu
      Sigma_out <- Sigma
    }
    if (length(observed_dims) == 1) {
      Sigma_out <- Sigmauu - (1 / Sigmaoo) * Sigmauo %*% t(Sigmaou)
      mu_out <- muu + Sigmauo * (y_observed_dims - muo) / Sigmaoo
    }
    if ((length(observed_dims) > 1) & (length(observed_dims) < length(mu))) {
      Sigmaoo_svd <- svd(Sigmaoo)
      Uoo <- Sigmaoo_svd$u
      Doo <- Sigmaoo_svd$d
      Doo_inv <- diag(1 / (Sigmaoo_svd$d))
      Sigmaooinv <- Uoo %*% Doo_inv %*% t(Uoo)
      mu_out <- muu + Sigmauo %*% Sigmaooinv %*% (y_observed_dims - muo)
      Sigma_out <- Sigmauu - Sigmauo %*% Sigmaooinv %*% Sigmaou
    }
  }
  
  list(mu_conditional = mu_out, Sigma_conditional = Sigma_out)
}

# ---- Moments of centered truncated Gaussian ----
moment_cal_func_centered <- function(a_vec, b_vec, Sigma){
  dimdat <- length(a_vec)
  stopifnot(all(dim(Sigma) == dimdat))
  stopifnot(length(a_vec) == length(b_vec))
  
  prob_andd_density_function_F <- function(x, a_vec, b_vec, pos, Sigma){ # F_k(x)
    if (length(a_vec) >= 2) {
      a_vec_min_pos <- a_vec[-pos]
      b_vec_min_pos <- b_vec[-pos]
      x_vec <- a_vec; x_vec[pos] <- x
      sigma <- Sigma[pos, pos]
      Sigma_cond_pos <- cond_mean_var_func(x_vec, rep(0, length(x_vec)), Sigma, c(pos))$Sigma_conditional
      mu_cond_pos <- cond_mean_var_func(x_vec, rep(0, length(x_vec)), Sigma, c(pos))$mu_conditional
      out <- dnorm(x, 0, sigma) * my_pmvnorm(a_vec_min_pos, b_vec_min_pos, mean = as.vector(mu_cond_pos), sigma = Sigma_cond_pos)
    } else {
      sigma <- as.numeric(Sigma)
      out <- dnorm(x, 0, sigma)
    }
    out
  }
  
  prob_andd_density_function_FF <- function(x, y, a_vec, b_vec, pos1, pos2, Sigma){ # F_{k,q}(x,y)
    if (length(a_vec) >= 3) {
      if (pos1 == pos2) stop("positions cannot be same")
      pos <- c(pos1, pos2)
      a_vec_min_pos <- a_vec[-pos]
      b_vec_min_pos <- b_vec[-pos]
      x_vec <- a_vec; x_vec[pos] <- c(x, y)
      sigma <- Sigma[pos, pos]
      Sigma_cond_pos <- cond_mean_var_func(x_vec, rep(0, length(x_vec)), Sigma, c(pos))$Sigma_conditional
      mu_cond_pos <- cond_mean_var_func(x_vec, rep(0, length(x_vec)), Sigma, c(pos))$mu_conditional
      
      if (length(a_vec_min_pos) == 1) {
        if (any(c(x,y) == -Inf) || any(c(x,y) == Inf)) {
          out <- 0
        } else {
          out <- mvtnorm::dmvnorm(c(x,y), rep(0,2), sigma) *
            (pnorm(b_vec_min_pos, mu_cond_pos, Sigma_cond_pos) - pnorm(a_vec_min_pos, mu_cond_pos, Sigma_cond_pos))
        }
      } else {
        if (any(c(x,y) == -Inf) || any(c(x,y) == Inf) ||
            any(a_vec_min_pos == -Inf) || any(a_vec_min_pos ==  Inf) ||
            any(b_vec_min_pos == -Inf) || any(b_vec_min_pos ==  Inf)) {
          out <- 0
        } else {
          out <- mvtnorm::dmvnorm(c(x,y), rep(0,2), sigma) *
            my_pmvnorm(a_vec_min_pos, b_vec_min_pos, as.vector(mu_cond_pos), Sigma_cond_pos)
        }
      }
    } else {
      if (any(c(x,y) == -Inf) || any(c(x,y) == Inf)) {
        out <- 0
      } else {
        sigma <- Sigma
        out <- mvtnorm::dmvnorm(c(x,y), rep(0,2), sigma)
      }
    }
    out
  }
  
  alpha <- my_pmvnorm(a_vec, b_vec, mean = rep(0, length(a_vec)), sigma = Sigma)
  mean_vec <- numeric(dimdat)
  Sigma_raw_mat <- matrix(0, dimdat, dimdat)
  
  temp_vec_1 <- numeric(dimdat)
  temp_vec_2 <- numeric(dimdat)
  temp_vec_3 <- numeric(dimdat)
  temp_vec_4 <- numeric(dimdat)
  
  for (i in 1:dimdat) {
    temp_vec_1[i] <- prob_andd_density_function_F(a_vec[i], a_vec, b_vec, i, Sigma)
    temp_vec_2[i] <- prob_andd_density_function_F(b_vec[i], a_vec, b_vec, i, Sigma)
    
    if (a_vec[i] == -Inf || a_vec[i] == Inf) {
      temp_vec_3[i] <- 0
    } else {
      temp_vec_3[i] <- a_vec[i] * prob_andd_density_function_F(a_vec[i], a_vec, b_vec, i, Sigma)
    }
    
    if (b_vec[i] == -Inf || b_vec[i] == Inf) {
      temp_vec_4[i] <- 0
    } else {
      temp_vec_4[i] <- b_vec[i] * prob_andd_density_function_F(a_vec[i], a_vec, b_vec, i, Sigma)
    }
  }
  
  for (i in 1:dimdat) {
    mean_vec[i] <- Sigma[i,] %*% (temp_vec_1 - temp_vec_2) / alpha
  }
  
  temp_mat_1 <- matrix(0, dimdat, dimdat)
  temp_mat_2 <- matrix(0, dimdat, dimdat)
  for (k in 1:dimdat) {
    for (q in 1:dimdat) {
      if (k == q) {
        temp_mat_1[k,q] <- 0
        temp_mat_2[k,q] <- 0
      } else {
        temp_mat_1[k,q] <- prob_andd_density_function_FF(a_vec[k], a_vec[q], a_vec, b_vec, k, q, Sigma) +
          prob_andd_density_function_FF(b_vec[k], b_vec[q], a_vec, b_vec, k, q, Sigma)
        temp_mat_2[k,q] <- prob_andd_density_function_FF(a_vec[k], b_vec[q], a_vec, b_vec, k, q, Sigma) +
          prob_andd_density_function_FF(b_vec[k], a_vec[q], a_vec, b_vec, k, q, Sigma)
      }
    }
  }
  
  temp_fn <- function(i, j, dimdat, temp_mat_1, temp_mat_2, Sigma){
    out_k <- 0
    for (k in 1:dimdat) {
      out_q <- 0
      for (q in 1:dimdat) {
        if (q != k) {
          out_q <- out_q + (Sigma[j,q] - Sigma[k,q] * Sigma[j,k] / Sigma[k,k]) *
            (temp_mat_1[k,q] - temp_mat_2[k,q])
        }
      }
      out_k <- out_k + Sigma[i,k] * out_q
    }
    out_k
  }
  
  Sigma_diag_inv <- matrix(0, dimdat, dimdat)
  diag(Sigma_diag_inv) <- (temp_vec_3 - temp_vec_4) / diag(Sigma)
  
  for (i in 1:dimdat) {
    for (j in 1:dimdat) {
      Sigma_raw_mat[i,j] <- Sigma[i,j] +
        Sigma[i,] %*% Sigma_diag_inv %*% Sigma[,j] / alpha +
        temp_fn(i, j, dimdat, temp_mat_1, temp_mat_2, Sigma) / alpha
    }
  }
  
  list(mean_vec = mean_vec, Sigma_raw_mat = Sigma_raw_mat)
}

# ---- Censoring dictionary (note: as written in Rmd this has issues; keep as-is) ----
## Maps a vector of size 3^{dim} -> censoring code per coordinate
cens_dict <- function(cens_vec){
  stopifnot(is.numeric(cens_vec), length(y) > 0)   # NOTE: 'y' not defined here in original; preserved
  dimdat <- log(length(cens_vec), base = 3)
  idx0 <- which(cens_vec != 0) - 1
  x <- numeric(d)                                  # NOTE: 'd' not defined in original; preserved
  for (index in d:1) {
    x[i] <- idx0 %% 3                              # NOTE: 'i' not defined in original; preserved
    idx0 <- idx0 %/% 3
  }
  x
}

# ---- Censoring conversion ----
## @return list(censor_indicator_left, censor_indicator_right, cens_lim_l_vec, cens_lim_u_vec)
cens_conversion <- function(cens_mat, y){
  cens_mat_new <- cbind(cens_mat[, c(1,2)], apply(cens_mat[, -c(1,2)], cens_dict, 1))
  cens_ind_left <- list()
  cens_ind_right <- list()
  y_mat <- do.call(rbind, y)
  conversion_fn <- function(num, index){
    stopifnot(num %in% c(0,1,2))
    if (index == 1) 1 * (num == 1) else 1 * (num == 2)
  }
  for (tt in 1:TT) {
    cens_ind_left[[tt]]  <- conversion_fn(cens_mat_new[(cens_mat_new[,1] == tt), -c(1,2)], 1)
    cens_ind_right[[tt]] <- conversion_fn(cens_mat_new[(cens_mat_new[,1] == tt), -c(1,2)], 2)
  }
  y_lim_upper <- numeric(dimdat)
  y_lim_lower <- numeric(dimdat)
  for (dd in 1:dimdat) {
    y_lim_lower[dd] <- mean(y_mat[(cens_mat_new[, dd+2] == 1), dd])
    y_lim_upper[dd] <- mean(y_mat[(cens_mat_new[, dd+2] == 2), dd])
  }
  cens_list <- list(
    censor_indicator_left  = cens_ind_left,
    censor_indicator_right = cens_ind_right,
    cens_lim_l_vec = y_lim_lower,
    cens_lim_u_vec = y_lim_upper
  )
  return(cens_list)
}

# ---- Censored conditional normal per (ii, tt) ----
cens_cond_normal <- function(ii, tt, y, X, censor_indicator_left, 
                             censor_indicator_right, cens_lim_l_vec,
                             cens_lim_u_vec, numclust, mu_list, sigma_list){
  dimdat <- ncol(y[[1]])
  ntlist <- sapply(y, nrow)
  
  ## Basic checks
  stopifnot(is.list(y))
  stopifnot(is.list(censor_indicator_left))
  stopifnot(is.list(censor_indicator_right))
  stopifnot(is.list(mu_list))
  stopifnot(is.list(sigma_list))
  stopifnot(length(cens_lim_l_vec) == length(cens_lim_u_vec))
  
  mu_conditional <- vector("list", numclust)
  Sigma_conditional <- vector("list", numclust)
  prob_cens_conditional <- vector("list", numclust)
  conditional_response <- vector("list", numclust)
  all_cond_mean_list <- vector("list", numclust)
  all_cond_second_moment_list <- vector("list", numclust)
  new_response_list <- vector("list", numclust)
  y_vec <- y[[tt]][ii,]
  
  if ((sum(censor_indicator_left[[tt]][ii,], na.rm = TRUE) != dimdat) |
      (sum(censor_indicator_right[[tt]][ii,], na.rm = TRUE) != dimdat)) {
    
    left_cens_index  <- which(censor_indicator_left[[tt]][ii,]  == 1)
    right_cens_index <- which(censor_indicator_right[[tt]][ii,] == 1)
    uncensored_index <- which((censor_indicator_left[[tt]][ii,] != 1 | is.na(censor_indicator_left[[tt]][ii,])) &
                                (censor_indicator_right[[tt]][ii,] != 1 | is.na(censor_indicator_right[[tt]][ii,])))
    
    lower_limits <- cens_lim_l_vec[left_cens_index]
    upper_limits <- cens_lim_u_vec[right_cens_index]
    
    for (iclust in 1:numclust) {
      mu <- mu_list[[iclust]][tt,]
      Sigma <- sigma_list[[iclust]]
      
      mu_conditional[[iclust]] <- cond_mean_var_func(y_vec, mu, Sigma, uncensored_index)$mu_conditional
      Sigma_conditional[[iclust]] <- cond_mean_var_func(y_vec, mu, Sigma, uncensored_index)$Sigma_conditional
      
      if ((length(upper_limits) > 0) & (length(lower_limits) > 0)) {
        p_lower_limit <- numeric(dimdat)
        p_upper_limit <- numeric(dimdat)
        p_lower_limit[left_cens_index]  <- -Inf
        p_lower_limit[uncensored_index] <- -Inf
        p_lower_limit[right_cens_index] <- upper_limits
        p_upper_limit[left_cens_index]  <- lower_limits
        p_upper_limit[uncensored_index] <- Inf
        p_upper_limit[right_cens_index] <- Inf
        
        if (length(uncensored_index) == 0) {
          prob_cens_conditional[[iclust]] <- my_pmvnorm(
            p_lower_limit[sort(c(left_cens_index, right_cens_index))],
            p_upper_limit[sort(c(left_cens_index, right_cens_index))],
            mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]]
          )[1]
          
          astar_vec <- as.vector(p_lower_limit[sort(c(left_cens_index, right_cens_index))] - mu_conditional[[iclust]])
          bstar_vec <- as.vector(p_upper_limit[sort(c(left_cens_index, right_cens_index))] - mu_conditional[[iclust]])
          
          all_cond_mean_list[[iclust]] <- mu_conditional[[iclust]] +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec
          
          all_cond_second_moment_list[[iclust]] <-
            mu_conditional[[iclust]] %*% t(mu_conditional[[iclust]]) +
            mu_conditional[[iclust]] %*% t(moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec) +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec %*% t(mu_conditional[[iclust]]) +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$Sigma_raw_mat
        } else {
          prob_cens_conditional[[iclust]] <- my_pmvnorm(
            p_lower_limit[-uncensored_index],
            p_upper_limit[-uncensored_index],
            mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]]
          )[1]
          
          astar_vec <- as.vector(p_lower_limit[-uncensored_index] - mu_conditional[[iclust]])
          bstar_vec <- as.vector(p_upper_limit[-uncensored_index] - mu_conditional[[iclust]])
          
          all_cond_mean_list[[iclust]] <- mu_conditional[[iclust]] +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec
          
          all_cond_second_moment_list[[iclust]] <-
            mu_conditional[[iclust]] %*% t(mu_conditional[[iclust]]) +
            mu_conditional[[iclust]] %*% t(moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec) +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec %*% t(mu_conditional[[iclust]]) +
            moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$Sigma_raw_mat
        }
      }
      
      if ((length(upper_limits) == 0) & (length(lower_limits) > 0)) {
        prob_cens_conditional[[iclust]] <- my_pmvnorm(rep(-Inf, length(lower_limits)), lower_limits,
                                                      mean = as.vector(mu_conditional[[iclust]]),
                                                      sigma = Sigma_conditional[[iclust]])[1]
        astar_vec <- rep(-Inf, length(lower_limits))
        bstar_vec <- as.vector(lower_limits - mu_conditional[[iclust]])
        
        all_cond_mean_list[[iclust]] <- mu_conditional[[iclust]] +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec
        
        all_cond_second_moment_list[[iclust]] <-
          mu_conditional[[iclust]] %*% t(mu_conditional[[iclust]]) +
          mu_conditional[[iclust]] %*% t(moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec) +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec %*% t(mu_conditional[[iclust]]) +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$Sigma_raw_mat
      }
      
      if ((length(upper_limits) > 0) & (length(lower_limits) == 0)) {
        prob_cens_conditional[[iclust]] <- my_pmvnorm(upper_limits, rep(Inf, length(upper_limits)),
                                                      mean = as.vector(mu_conditional[[iclust]]),
                                                      sigma = Sigma_conditional[[iclust]])
        
        astar_vec <- as.vector(upper_limits - mu_conditional[[iclust]])
        bstar_vec <- rep(Inf, length(upper_limits))
        
        all_cond_mean_list[[iclust]] <- mu_conditional[[iclust]] +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec
        
        all_cond_second_moment_list[[iclust]] <-
          mu_conditional[[iclust]] %*% t(mu_conditional[[iclust]]) +
          mu_conditional[[iclust]] %*% t(moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec) +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$mean_vec %*% t(mu_conditional[[iclust]]) +
          moment_cal_func_centered(astar_vec, bstar_vec, Sigma_conditional[[iclust]])$Sigma_raw_mat
      }
      
      if ((length(upper_limits) == 0) & (length(lower_limits) == 0)) {
        prob_cens_conditional[[iclust]] <- NA
        all_cond_second_moment_list[[iclust]] <- NA
        all_cond_mean_list[[iclust]] <- mu_conditional[[iclust]]
      }
      
      new_response_list[[iclust]] <- numeric(dimdat)
      new_response_list[[iclust]][uncensored_index] <- y_vec[uncensored_index]
      new_response_list[[iclust]][-uncensored_index] <- all_cond_mean_list[[iclust]]
    }
  }
  
  list(
    new_response_list = new_response_list,
    all_conditional_means = all_cond_mean_list,
    all_conditional_second_moment_list = all_cond_second_moment_list
  )
}

# ---- Initializations for checks ----
mn_list <- list()
Sig_list <- list()
numclust <- 2
for (iclust in 1:numclust) {
  mn_list[[iclust]] <- X %*% t(matrix(rnorm(dimdatt * ncol(X), iclust/2, 1), dimdatt, ncol(X)))
  Sig_list[[iclust]] <- (iclust/2) * diag(dimdatt)
}

cens_vec_1 <- rep(-2, dimdatt)
cens_vec_2 <- rep(6.8, dimdatt)

## Example debug hooks (commented)
## debug(moment_cal_func_centered)
## debug(cens_cond_normal)
## cens_cond_normal(ii = 73, tt = 1, y = ylist, X,
##                  censor_indicator_left = cens_ind_left,
##                  censor_indicator_right = cens_ind_right,
##                  cens_lim_l_vec = cens_vec_1,
##                  cens_lim_u_vec = cens_vec_2,
##                  numclust = 2, mu_list = mn_list, sigma_list = Sig_list)
