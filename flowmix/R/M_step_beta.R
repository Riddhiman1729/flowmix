# M_step_beta
# Author: Riddhiman Bhattacharya
# Date: 2025-10-21

# --- Libraries ---
library(dplyr)
library(tibble)
library(mvtnorm)
library(Rcpp)
library(assertthat)
library(Matrix)

# --- Mstep_alpha --------------------------------------------------------------

#' Solves the M step for estimating the coefficients for the cluster
#' probabilities (multinomial logit with lasso via solve_multinom()).
#'
#' Note, the estimated alphas are identifiable only up to a constant shift.
#' @param resp Responsibilities; list length T of nt x K matrices OR a T x K matrix of sums.
#' @param X Covariate matrix (T x p).
#' @param lambda Regularization parameter.
#' @param zerothresh Values below \code{zerothresh} are set to zero.
#' @param numclust Number of clusters.
#' @return list(prob = TT x K, alpha = K x (p+1))
Mstep_alpha <- function(resp, X, numclust, lambda, zerothresh = 1e-8) {
  TT <- nrow(X)
  p  <- ncol(X)
  
  ## If resp provided as list-of-matrices, collapse to T x K
  if (is.list(resp)) {
    resp.sum <- t(sapply(resp, colSums))    # (T x K)
  } else {
    resp.sum <- as.matrix(resp)             # already (T x K)
  }
  resp.sum <- as.matrix(resp.sum)
  stopifnot(all(dim(resp.sum) == c(TT, numclust)))
  
  ## Fit (user-provided) multinomial with l1 penalty
  alpha <- solve_multinom(resp.sum, X, lambda)
  
  ## Threshold tiny values
  alpha[which(abs(alpha) < zerothresh, arr.ind = TRUE)] <- 0
  alpha <- t(as.matrix(alpha))
  stopifnot(all(dim(alpha) == c(numclust, p + 1)))
  
  ## Fitted probabilities
  Xa <- cbind(1, X)
  probhatmat <- as.matrix(exp(Xa %*% t(alpha)))
  probhat    <- probhatmat / rowSums(probhatmat)
  
  stopifnot(all(dim(probhat) == c(TT, numclust)))
  if (any(is.na(probhat)) || !all(probhat >= 0)) stop("probhat was erroneous")
  
  list(prob = probhat, alpha = alpha)
}

# --- Mstep_beta (CVXR/backup) -------------------------------------------------

#' Backup M step for beta via lasso (slower than ADMM version).
Mstep_beta <- function(resp, ylist, X,
                       mean_lambda = 0,
                       sigma,
                       maxdev = NULL,
                       sigma_eig_by_clust = NULL,
                       first_iter = FALSE,
                       cvxr_ecos_thresh = 1e-8,
                       cvxr_scs_eps = 1e-5,
                       zerothresh = 1e-8) {
  
  TT       <- length(ylist)
  numclust <- ncol(resp[[1]])
  dimdat   <- ncol(ylist[[1]])
  ntlist   <- sapply(ylist, nrow)
  resp.sum <- t(sapply(resp, colSums))   # (T x K)
  resp.sum <- as.matrix(resp.sum)
  N        <- sum(resp.sum)
  p        <- ncol(X)
  Xa       <- cbind(1, X)
  
  if (!is.null(sigma)) {
    assertthat::assert_that(all.equal(dim(sigma), c(numclust, dimdat, dimdat)) == TRUE)
  }
  
  ## Setup
  manip_obj <- manip(ylist, Xa, resp, sigma, numclust,
                     sigma_eig_by_clust = sigma_eig_by_clust,
                     first_iter = first_iter)
  Xtildes <- manip_obj$Xtildes
  yvecs   <- manip_obj$yvecs
  
  ## Intercepts excluded from penalty
  exclude.from.penalty <- (0:(dimdat - 1)) * (ncol(Xa)) + 1
  
  results <- lapply(seq_len(numclust), function(iclust) {
    betahat <- cvxr_lasso(
      X  = Xtildes[[iclust]],
      Xorig = X,
      y  = yvecs[[iclust]],
      lambda = mean_lambda,
      exclude.from.penalty = exclude.from.penalty,
      maxdev = maxdev,
      dimdat = dimdat,
      N  = N,
      ecos_thresh = cvxr_ecos_thresh,
      scs_eps = cvxr_scs_eps
    )
    betahat[which(abs(betahat) < zerothresh, arr.ind = TRUE)] <- 0
    yhat <- Xa %*% betahat
    assertthat::assert_that(all(dim(betahat) == c(p + 1, dimdat)))
    list(betahat = betahat, yhat = as.matrix(yhat))
  })
  
  betahats     <- lapply(results, \(a) a$betahat)
  yhats        <- lapply(results, \(a) a$yhat)
  yhats_array  <- array(NA, dim = c(TT, dimdat, numclust))
  for (k in seq_len(numclust)) yhats_array[, , k] <- yhats[[k]]
  
  list(beta = betahats, mns = yhats_array)
}

# --- manip helper -------------------------------------------------------------

#' Prepares Xtildes / Ytildes / yvecs per cluster.
manip <- function(ylist, X, resp, sigma, numclust,
                  sigma_eig_by_clust = NULL,
                  first_iter = FALSE) {
  
  ntlist  <- sapply(ylist, nrow)
  dimdat  <- ncol(ylist[[1]])
  TT      <- nrow(X)
  resp.sum <- t(sapply(resp, colSums))   # (T x K)
  sigma.inv.halves <- array(NA, dim = dim(sigma))
  
  if (first_iter) {
    for (k in seq_len(numclust)) sigma.inv.halves[k, , ] <- mtsqrt_inv(sigma[k, , ])
  } else {
    for (k in seq_len(numclust)) sigma.inv.halves[k, , ] <- sigma_eig_by_clust[[k]]$inverse_sigma_half
  }
  
  emptymat <- matrix(0, nrow = dimdat, ncol = TT)
  Ytildes  <- vector("list", numclust)
  for (k in seq_len(numclust)) {
    for (tt in 1:TT) emptymat[, tt] <- colSums(resp[[tt]][, k] * ylist[[tt]])
    Ytildes[[k]] <- emptymat
  }
  Xtildes <- lapply(seq_len(numclust), function(k) {
    sigma.inv.halves[k, , ] %x% (sqrt(resp.sum[, k]) * X)
  })
  
  yvecs <- lapply(seq_len(numclust), function(k) {
    yvec <- (1 / sqrt(resp.sum[, k]) * t(Ytildes[[k]])) %*% sigma.inv.halves[k, , ]
    as.vector(yvec)
  })
  
  list(Xtildes = Xtildes, Ytildes = Ytildes, yvecs = yvecs)
}

# --- multinom objective (used by solve_multinom) ------------------------------

multinom_objective <- function(alpha, x, y, lambda, N, exclude.from.penalty = NULL) {
  n <- nrow(x); p <- ncol(x); L <- ncol(y)
  eta <- x %*% alpha
  v <- 1:p
  if (!is.null(exclude.from.penalty)) {
    stopifnot(all(exclude.from.penalty %in% (1:p)))
    v <- (1:p)[-exclude.from.penalty]
  }
  ys <- rowSums(y)
  (1 / N) * (sum(eta * y) - sum(ys * log(rowSums(exp(eta))))) -
    lambda * sum(abs(alpha[v, ]))
}

# --- LA-ADMM wrapper ----------------------------------------------------------

la_admm_oneclust <- function(K, ...) {
  args <- list(...)
  p <- args$p; TT <- args$TT; dimdat <- args$dimdat
  
  if (isTRUE(args$first_iter)) {
    beta <- matrix(0, nrow = p + 1, ncol = dimdat)
    Z    <- matrix(0, nrow = TT, ncol = dimdat)
    W    <- matrix(0, nrow = p, ncol = dimdat)
    U    <- matrix(0, nrow = TT + p, ncol = dimdat)
    args$beta <- beta; args$Z <- Z; args$W <- W; args$U <- U
  }
  
  objectives <- c()
  for (kk in 1:K) {
    if (kk > 1) {
      args$beta <- beta; args$Z <- Z; args$W <- W; args$U <- U
      args$rho  <- rho
    }
    args$outer_iter <- kk
    
    ## original code used eval(as.call(...)); keep that behavior:
    argn <- lapply(names(args), as.name); names(argn) <- names(args)
    call <- as.call(c(list(as.name("admm_oneclust")), argn))
    res  <- eval(call, args)
    
    objectives <- c(objectives, res$fit)
    padding <- 1e-12
    objectives <- objectives + padding
    
    if (outer_converge(stats::na.omit(objectives)) || res$converge) break
    
    rho  <- args$rho * 2
    beta <- res$beta; Z <- res$Z; W <- res$W; U <- res$U
  }
  
  res$kk <- kk
  res
}

outer_converge <- function(objectives) {
  consec <- 4
  if (length(objectives) < consec) return(FALSE)
  mytail <- utils::tail(objectives, consec)
  rel_diffs <- mytail[1:(consec - 1)] / mytail[2:consec]
  all(abs(rel_diffs - 1) < 1e-5)
}

# --- Core ADMM one cluster ----------------------------------------------------

admm_oneclust <- function(iclust, niter, Xtilde, yvec, p,
                          TT, N, dimdat, maxdev,
                          Xa,
                          rho,
                          rhoinit,
                          Xinv,
                          schurA,
                          schurB,
                          term3,
                          sigmainv,
                          Xaug,
                          ybar,
                          Q,
                          lambda,
                          resp,
                          resp.sum,
                          ylist, X, tX, err_rel, err_abs,
                          zerothresh,
                          beta, Z, W, U,
                          first_iter, outer_iter,
                          local_adapt,
                          sigma,
                          sigma_eig_by_clust,
                          space = 20) {
  
  resid_mat <- matrix(NA, nrow = ceiling(niter/5), ncol = 4)
  colnames(resid_mat) <- c("primresid", "primerr", "dualresid", "dualerr")
  zrows <- 1:TT
  wrows <- TT + (1:p)
  rhofac <- rho / rhoinit
  
  fits <- rep(NA, ceiling(niter / space))
  converge <- FALSE
  
  schurA <- myschur(schurA$orig * rhofac)
  TA <- schurA$T; TB <- schurB$T
  UA <- schurA$Q; UB <- schurB$Q
  tUA <- schurA$tQ; tUB <- schurB$tQ
  
  for (iter in 1:niter) {
    syl_C <- prepare_sylC_const3(U, Xaug, rho, Z, X, W, term3, as.matrix(sigma[iclust,,]), Xinv)
    F     <- (-1) * tUA %*% syl_C %*% UB
    beta  <- UA %*% matrix_function_solve_triangular_sylvester_barebones(TA, TB, F) %*% tUB
    beta  <- t(beta)
    if (any(is.nan(beta))) browser()
    
    Xbeta <- X %*% beta
    Z     <- Z_updateC(Xbeta, U[zrows,, drop=FALSE], maxdev, rho, dimdat, TT)
    W     <- W_update(beta, U[wrows,, drop=FALSE], lambda, rho)
    if (is.vector(W)) W <- cbind(W)
    U     <- U_update(U, rho, Xaug, beta, Z, W)
    
    if (iter > 1  && iter %% 5 == 0) {
      obj <- converge(beta, rho, W, Z, W_prev, Z_prev,
                      Uw = U[wrows,], Uz = U[zrows,], tX = tX,
                      Xbeta1 = Xbeta, err_rel = err_rel, err_abs = err_abs)
      
      jj <- (iter / 5)
      resid_mat[jj, ] <- c(norm(obj$primal_resid, "F"), obj$primal_err,
                           norm(obj$dual_resid, "F"), obj$dual_err)
      
      if (obj$converge) { converge <- TRUE; break }
    }
    
    W_prev <- W; Z_prev <- Z
  }
  
  beta <- W
  beta[which(abs(beta) < zerothresh, arr.ind = TRUE)] <- 0
  beta0 <- intercept(resp, resp.sum, ylist, beta, X, N, iclust, ybar)
  betafull <- rbind(beta0, beta)
  yhat <- Xa %*% betafull
  
  Q2  <- (1 / 2) * Q %*% (beta %*% sigmainv %*% t(beta))
  M   <- t(term3) %*% t(beta)
  fit <- sum(diag(Q2)) - sum(diag(M)) + lambda * sum(abs(beta) > zerothresh)
  
  list(beta = betafull, yhat = yhat, resid_mat = resid_mat, fits = fits,
       converge = converge, fit = fit, Z = Z, W = W, U = U)
}

U_update <- function(U, rho, Xaug, beta, Z, W) {
  U + rho * (Xaug %*% beta - rbind(Z, W))
}

W_update <- function(beta, Uw, lambda, rho) {
  soft_thresh(beta + Uw / rho, lambda / rho)
}

soft_thresh <- function(a, b) {
  sign(a) * pmax(0, abs(a) - b)
}

intercept <- function(resp, resp.sum, ylist, beta, X, N, iclust, ybar) {
  dimdat <- ncol(ylist[[1]])
  TT     <- length(ylist)
  resp.sum.thisclust <- sum(resp.sum[, iclust])
  mn <- X %*% beta
  yhat <- (resp.sum[, iclust] / resp.sum.thisclust) * mn
  ybar - colSums(yhat)
}

Z_update <- function(Xbeta, Uz, C, rho) {
  mat <- Xbeta + Uz / rho
  projCmat(mat, C)
}

projCmat <- function(mat, C) {
  if (!is.null(C)) {
    vlens <- sqrt(rowSums(mat * mat))
    inds  <- which(vlens > C)
    if (length(inds) > 0) mat[inds, ] <- mat[inds, ] * C / vlens[inds]
  }
  mat
}

center_X <- function(iclust, resp.sum, X) {
  resp.sum.thisclust <- sum(resp.sum[, iclust])
  Xtilde  <- colSums(resp.sum[, iclust] * X) / resp.sum.thisclust
  sweep(X, 2, Xtilde, check.margin = FALSE)
}

# --- Schur decomposition helper ----------------------------------------------

myschur <- function(mat) {
  stopifnot(nrow(mat) == ncol(mat))
  if (is.numeric(mat) && length(mat) == 1) mat <- as.matrix(mat)
  obj <- Matrix::Schur(mat)
  obj$tQ <- t(obj$Q)
  obj$orig <- mat
  obj
}

# --- ADMM convergence metrics -------------------------------------------------

converge <- function(beta1, rho, w, Z, w_prev, Z_prev, Uw, Uz, tX, Xbeta1,
                     err_rel = 1e-4, err_abs = 0) {
  
  prim1 <- rbind(beta1, Xbeta1)
  prim2 <- -rbind(w, Z)
  primal_resid <- prim1 + prim2
  dual_resid   <- -rho * ((w - w_prev) + (tX %*% (Z - Z_prev)))
  tAU <- Uw + tX %*% Uz
  
  primal_err <- sqrt(length(primal_resid)) * err_abs +
    err_rel * max(norm(prim1, "F"), norm(prim2, "F"))
  dual_err <- sqrt(length(dual_resid)) * err_abs +
    err_rel * norm(tAU, "F")
  
  primal_resid_size <- norm(primal_resid, "F")
  dual_resid_size   <- norm(dual_resid, "F")
  primal_converge <- (primal_resid_size <= primal_err)
  dual_converge   <- (dual_resid_size   <= dual_err)
  
  assertthat::assert_that(is.numeric(primal_resid_size))
  assertthat::assert_that(is.numeric(primal_err))
  assertthat::assert_that(is.numeric(dual_resid_size))
  assertthat::assert_that(is.numeric(dual_err))
  
  list(primal_resid = primal_resid,
       primal_err   = primal_err,
       dual_resid   = dual_resid,
       dual_err     = dual_err,
       converge     = (primal_converge & dual_converge))
}

# --- Per-cluster objective (optional) -----------------------------------------

objective_per_cluster <- function(beta, ylist, Xa, resp, lambda, N, dimdat,
                                  iclust, sigma, iter, zerothresh, first_iter,
                                  sigma_eig_by_clust = NULL, rcpp = FALSE) {
  TT <- length(ylist)
  p  <- ncol(Xa) - 1
  mn <- Xa %*% beta
  
  if (first_iter) {
    sigma_half <- sigma_half_from_eig(eigendecomp_sigma_barebones(sigma[iclust,,]))
  } else {
    sigma_eig  <- sigma_eig_by_clust[[iclust]]
    sigma_half <- sigma_eig$inverse_sigma_half
  }
  
  if (!rcpp) {
    wt_resids_list <- lapply(1:TT, function(tt) {
      y <- ylist[[tt]]
      mumat <- matrix(mn[tt, ], ncol = ncol(y), nrow = nrow(y), byrow = TRUE)
      sqrt(resp[[tt]][, iclust]) * (y - mumat)
    })
    transformed_resids <- do.call(rbind, wt_resids_list) %*% sigma_half
    resid.prods <- rowSums(transformed_resids * transformed_resids)
  } else {
    wt_resids_list <- lapply(1:TT, function(tt) {
      subtractC3(sqrt(resp[[tt]][, iclust]), ylist[[tt]], mn[tt, ])
    })
    transformed_resids <- do.call(rbind, wt_resids_list) %*% sigma_half
    resid.prods <- rowSumsC2_arma(transformed_resids)
  }
  
  grand.sum <- sum(resid.prods)
  obj <- (1/(2 * N)) * grand.sum + lambda * sum(abs(beta[-1, ]) > zerothresh)
  obj
}

# --- RUN SECTION (from the Rmd "Run above this") ------------------------------

## Estep_y() and cens_cond_normal(), my_pmvnorm(), cond_mean_var_func() should
## be defined/loaded beforehand (from your other script).

# Example:
# estepy_saved <- Estep_y(
#   y = ylist, X,
#   censor_indicator_left  = cens_ind_left,
#   censor_indicator_right = cens_ind_right,
#   cens_lim_l_vec = cens_vec_1,
#   cens_lim_u_vec = cens_vec_2,
#   numclust = 2,
#   mu_list = mn_list,
#   sigma_list = Sig_list
# )
# y_new <- estepy_saved$new_responses
#
# resp <- Estep(
#   mn_array, sigma_array, pmat, ylist,
#   censor_indicator_left  = cens_ind_left,
#   censor_indicator_right = cens_ind_right,
#   cens_lim_l_vec = cens_vec_1,
#   cens_lim_u_vec = cens_vec_2,
#   numclust = 2,
#   denslist_by_clust = NULL,
#   first_iter = TRUE,
#   eps = 1e-20,
#   countslist = NULL
# )

## Compile C++ helpers (adjust paths if needed)
# Rcpp::sourceCpp("admm.cpp", verbose = TRUE, showOutput = TRUE, rebuild = TRUE)
# Rcpp::sourceCpp("syl.cpp",  verbose = TRUE, showOutput = TRUE, rebuild = TRUE)
# Rcpp::sourceCpp("estep.cpp", verbose = TRUE, showOutput = TRUE, rebuild = TRUE)

## Run Mstep_beta_admm (requires objects created above)
# Mstep_beta_admm(
#   resp,
#   ylist, ynew = y_new,
#   X,
#   mean_lambda = 0,
#   sigma = aperm(simplify2array(Sig_list), perm = c(3,1,2)),
#   sigma_eig_by_clust = NULL,
#   first_iter = TRUE,
#   censor_indicator_left  = cens_ind_left,
#   censor_indicator_right = cens_ind_right,
#   cens_lim_l_vec = cens_vec_1,
#   cens_lim_u_vec = cens_vec_2,
#   betas = NULL, Zs = NULL, Ws = NULL, Us = NULL,
#   maxdev = dimdatt,
#   niter = 1e4, rho = 100,
#   err_rel = 1e-3, err_abs = 0,
#   zerothresh = 1e-6,
#   local_adapt = FALSE,
#   local_adapt_niter = 5,
#   space = 50
# )
