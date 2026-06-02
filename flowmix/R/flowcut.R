##' Main function for our method. Repeats the EM algorithm with |nrep| restarts
##' (5 by default).
##'
##' @param ... Arguments for \code{flowcut_once()}.
##' @param nrep Number of restarts.
##'
##' @return The |flowcut| class object that had the best likelihood.
##' @export
flowcut <- function(..., nrep = 5, mc.cores = 1){
  
  dots <- list(...)
  if("verbose" %in% names(dots)){
    if(dots$verbose){
      cat("EM will restart", nrep, "times", fill=TRUE)
    }
  }
  
  ## Don't do many restarts if warmstart-able mean is provided
  if(!is.null(dots$mn)) nrep = 1
  if(!is.null(dots$seed)) stop("Can't provide seed for flowcut()! Only for flowcut_once().")

  ## Main loop
  reslist <- parallel::mclapply(1:nrep, function(irep){
    if("verbose" %in% names(dots)){
      if(dots$verbose){ cat("EM restart:", irep, fill=TRUE) }
    }
    if("verbose" %in% names(dots)){
      if(dots$verbose){
        cat(fill=TRUE)
      }
    }
    flowcut_once(...)
  }, mc.cores = mc.cores)

  
  ## Pick the best one and return
  objlist = lapply(reslist, function(res){ res$obj[-1]})
  ii = which.min(sapply(objlist, min))
  final_model = reslist[[ii]]
  
  ## Also save /all/ the objectives
  final_model$all_objectives =
    lapply(1:nrep, function(irep){
      one_model = reslist[[irep]]
      data.frame(objective = one_model$objectives) %>%
        mutate(iter=row_number(), irep = irep) %>%
        select(irep, iter, objective)
    }) %>% bind_rows()
  
  return(final_model)
}

##' Main function for running the EM algorithm once.
##' @noRd
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }


###CHECK IF THE OBJECTIVE IS CORRECT. IF NOT CORRECT IT.
##' Main function for running the EM algorithm once.
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   d), which contains coordinates of the d-variate particles, organized over
##'   time (T) and with (nt) particles at every time point.
##' @param countslist Multiplicity for particles in \code{ylist}.
##' @param numclust Number of clusters
##' @param X Matrix of size (T x p+1)
##' @param cens_left_list A T-lengthed list of d-column, nt-row boolean
##'   matrices; the (i,j)'th entry of each matrix marks whether the $i$-th row
##'   of the ylist[[tt]] is left-censored or not.
##' @param cens_right_list A T-lengthed list of d-column, nt-row boolean
##'   matrices; the (i,j)'th entry of each matrix marks whether the $i$-th row
##'   of the ylist[[tt]] is right-censored or not.
##' @param mean_lambda lambda for lasso for the mean.
##' @param prob_lambda lambda for lasso for probabilities.
##' @param tol_em Relative tolerance for EM convergence. Defaults to 1E-4.
##' @param zero_stabilize Defaults to FALSE. If TRUE, the EM is only run until
##'   the pattern of zeros in the coefficients stabilizes over EM iterations.
##' @param seed Seven integers that is called directly before \code{init_mn()}
##'   so it can be assigned to \code{.Random.seed} for setting the random state
##'   for the initial mean generation.
##' @param niter Number of EM iterations.
##' @param mn Initial means to use. Defaults to NULL.
##' @param verbose TRUE for loudness (e.g. printing EM iterations).
##' @param maxdev Radius for maximum deviation of cluster means over time.
##' @param admm_rho Step size for ADMM
##' @param admm_err_rel Relative error threshold for stopping ADMM.
##' @param admm_err_abs Absolute error threshold for stopping ADMM.
##' @param admm_local_adapt_niter Absolute error threshold for stopping ADMM.
##' @param admm_niter Number of ADMM iterations.
##' @param admm_local_adapt TRUE if locally adaptive ADMM (LA-ADMM) is to be
##'   used. If so, \code{admm_niter} becomes the inner number of iterations, and
##'   \code{admm_local_adapt_niter} becomes the number of outer iterations.
##' @param admm_local_adapt_niter Number of inner iterations in LA ADMM.
##' @param CVXR If TRUE, use CVXR instead of ADMM. Slow, and meant to be used
##'   only for sanity checking during code development.
##' @param flatX_thresh Threshold for detecting if any covariates are flat (low
##'   variance). These flat coefficients will have be set to zero and excluded
##'   from estimation altogether.
##' @param sigma_fac Defaults to 1, and governs how big the initial covariance
##'   matrices should be in size.
##' @param countslist_overwrite Mainly used by \code{cv.flowmix()}; not to be
##'   modified by user.
##' @param zerothresh Alpha coefficient values below \code{zerothresh} are set
##'   to zero.
##'
##' @return List containing fitted parameters and means and mixture weights,
##'   across algorithm iterations. \code{beta} is a list of (p+1 x dimdat)
##'   arrays. \code{alpha} is a (numclust x (p+1)) array.
##'
##' @export
flowcut_once <- function(ylist, X,
                         ## cens_left_list,
                         ## cens_right_list,
                         cens_indicator_list_left,
                         cens_indicator_list_right,
                         cens_lim_vec_lower,
                         cens_lim_vec_upper,
                         countslist = NULL,
                         numclust, niter = 1000,
                         mn = NULL, prob_lambda,
                         mean_lambda, verbose = FALSE,
                         sigma_fac = 1, tol_em = 1E-4,
                         maxdev = NULL,
                         countslist_overwrite = NULL,
                         zero_stabilize  = FALSE,
                         zerothresh = 1E-6,
                         ## beta Mstep (ADMM) settings
                         admm_rho = 0.01,
                         admm_err_rel = 1E-3,
                         admm_err_abs = 1E-4,
                         ## beta M step (Locally Adaptive ADMM) settings
                         admm_local_adapt = TRUE,
                         admm_local_adapt_niter = 10,
                         admm_niter = (if(admm_local_adapt)1E3 else 1E4),
                         CVXR =FALSE, ## temporary
                         seed = NULL,
                         flatX_thresh = 1e-5
){

  
  . = NULL ## Fixing check()
  
  ## Capture all arguments once
  call <- sys.call();
  call[[1]] <- as.name('list');
  args <- eval.parent(call)
  
  ## Basic checks
  if(!is.null(maxdev)){
    assertthat::assert_that(maxdev!=0)
  } else {
    maxdev = 1E10 ## Some large number
  }
  assertthat::assert_that(!(is.data.frame(X)))
  assertthat::assert_that(sum(is.na(X)) == 0)
  assertthat::assert_that(length(ylist) == nrow(X))
  assertthat::assert_that(numclust > 1)
  assertthat::assert_that(niter > 1)
  
  ## assert_that(!(is.data.frame(ylist[[1]])))
  ## assertthat::assert_that(prob_lambda > 0)
  ## assertthat::assert_that(all(sapply(ylist, nrow) == sapply(countslist, length)))
  
  
  ## Detect if any covariates are flat (low variance)
  ## If so, remove that, run flowmix_once(), and add back zeros in the coefficients.
  flatX = which(apply(X, 2, stats::sd) < flatX_thresh)
  if(length(flatX) == 0) rm(args)
  if(length(flatX) > 0){
    
    warning("Some covariates are flat. The coefficients for these variables will be set to zero.")
    orig_names = colnames(X)
    
    ## Recursive call of flowmix_once
    args$X = X[, -flatX]
    argn <- lapply(names(args), as.name)
    names(argn) <- names(args)
    call <- as.call(c(list(as.name("flowmix_once")), argn))
    res = eval(call, args)
    
    ## Alter the call so that those coefficients are all zero.
    res$beta <- alter_beta(res$beta, flatX, orig_names)
    res$alpha <- alter_alpha(res$alpha, flatX, orig_names)
    
    ## Alter some other things
    res$p = res$p + 1
    res$X = X
    
    ## Return the result
    return(res)
  }
  
  ## Setup
  TT = length(ylist)
  dimdat = ncol(ylist[[1]])
  p = ncol(X)
  if(!is.null(seed)){
    assertthat::assert_that(all((seed %>% sapply(., class)) == "integer"))
    assertthat::assert_that(length(seed) == 7)
  }
  if(is.null(mn)) mn = init_mn(ylist, numclust, TT, dimdat, countslist, seed)
  ntlist = sapply(ylist, nrow)
  N = sum(ntlist)


  ## For censoring
  any_censoring = TRUE ## TODO incorporate as a flowcut() function input
  if(any_censoring){

    ## ## Basic checks
    ## assertthat::assert_that(is.null(cens_lim_vec_lower)==F)
    ## assertthat::assert_that(is.null(cens_lim_vec_upper)==F)
    ## assertthat::assert_that(length(cens_indicator_list_left)==nrow(X))
    ## assertthat::assert_that(length(cens_indicator_list_right)==nrow(X))
    ## assertthat::assert_that(length(cens_lim_vec_lower)==ncol(ylist[[1]]))
    ## assertthat::assert_that(length(cens_lim_vec_upper)==ncol(ylist[[1]]))

    ## convert to the main new objects we'lll use (separated!)
    convert_obj = my_convert(ylist, cens_indicator_list_left, cens_indicator_list_right)
    ylist_cens = convert_obj$ylist_cens
    ylist_uncens = convert_obj$ylist_uncens

    left_cens_list = convert_obj$left_cens_list
    right_cens_list = convert_obj$right_cens_list

    idx_cens_list = convert_obj$idx_cens_list
    idx_uncens_list = convert_obj$idx_uncens_list

  } else {

    ylist_cens = NULL
    ylist_uncens = ylist

    left_cens_list = NULL
    right_cens_list = NULL

    idx_cens_list = NULL
    idx_uncens_list = NULL
  }

  ## Initialize some objects
  prob = matrix(1/numclust, nrow = TT, ncol = numclust) ## Initialize to all 1/K.
  denslist_by_clust <- NULL
  objectives = c(+1E20, rep(NA, niter-1))
  sigma = init_sigma(ylist, numclust, sigma_fac) ## (T x numclust x dimdat x dimdat)
  sigma_eig_by_clust = NULL
  zero.betas = zero.alphas = list()
  admm_niters = list()
  
  ## Warm startable variables
  betas = NULL
  Zs = NULL
  wvecs = NULL
  uws = NULL
  Uzs = NULL
  
  ## New ADMM parameters.
  Zs = NULL
  Ws = NULL
  Us = NULL
  
  ## The least elegant solution I can think of.. used only for blocked cv
  if(!is.null(countslist_overwrite)) countslist = countslist_overwrite
  if(!is.null(countslist)) check_trim(ylist, countslist)


  ## Main loop
  start.time = Sys.time()
  for(iter in 2:niter){
    if(verbose){
      print(cat("Iteration", iter - 1, "of", niter - 1, "EM iterations.\n"))
    }
    
    ## Estep for flowcut
    resp <- Estep_flowcut(mn, sigma, prob,
                          ylist,
                          ylist_cens,
                          ylist_uncens,
                          left_cens_list, 
                          right_cens_list,
                          idx_cens_list,
                          cens_lim_vec_lower,
                          cens_lim_vec_upper,
                          numclust,
                          countslist = countslist)

    ## Conditional means of y and yy^T (given censored particles)
    estepy_saved <- Estep_y_flowcut(ylist_cens, ylist_uncens, X,
                                    left_cens_list, right_cens_list,
                                    cens_lim_vec_lower,
                                    cens_lim_vec_upper,
                                    numclust, mn,
                                    sigma)

    ## browser()
    new_responses = estepy_saved$new_responses %>% purrr::list_transpose()
    second_moments = estepy_saved$second_moments %>% purrr::list_transpose()


    ## iclust=1
    ## tt=1
    ## ylist_uncens[[1]] %>% plot(ylim = c(-1,2))
    ## ylist_cens[[1]] %>% points(col='green')
    ## iclust=1
    ## new_responses[[iclust]][[tt]] %>% points(col='blue')
    ## iclust=2
    ## new_responses[[iclust]][[tt]] %>% points(col='blue')
    ## cens_lim_vec_upper
    ## cens_lim_vec_lower

    ## mn %>% .[1,,]
    ## ylist_cens[[1]]
    ## new_responses[[iclust]][[tt]]
    ## [[iclust]][[2]] %>% plot()

    ## M step (three parts)
    ## 1. Alpha
    res.alpha = Mstep_alpha(resp, X,
                            numclust, lambda = prob_lambda,
                            zerothresh = zerothresh)
    
    prob = res.alpha$prob
    alpha = res.alpha$alpha
    rm(res.alpha)
    
    ## 2. Beta
    res.beta = Mstep_beta_admm(resp,
                               ylist_uncens,
                               ## New for flowcut
                               new_responses,
                               idx_cens_list,
                               idx_uncens_list,
                               ## End of new
                               X,
                               mean_lambda = mean_lambda,
                               first_iter = (iter == 2),
                               sigma_eig_by_clust = sigma_eig_by_clust,
                               sigma = sigma, maxdev = maxdev, rho = admm_rho,
                               betas = betas,
                               Zs = Zs,
                               Ws = Ws,
                               Us = Us,
                               err_rel = admm_err_rel,
                               err_abs = admm_err_abs,
                               niter = admm_niter,
                               local_adapt = admm_local_adapt,
                               local_adapt_niter = admm_local_adapt_niter)

    admm_niters[[iter]] = unlist(res.beta$admm_niters)
    
    ## Harvest means
    mn = res.beta$mns
    betas = beta = res.beta$beta
    
    ## Harvest other things for next iteration's ADMM.
    Zs = res.beta$Zs
    Ws = res.beta$Ws
    Us = res.beta$Us
    ## rm(res.beta)
    
    ## Check if the number of zeros in the alphas and betas have stabilized.
    zero.betas[[iter]] = lapply(beta, function(mybeta) which(mybeta==0))
    zero.alphas[[iter]] = which(alpha == 0)
    if(zero_stabilize & iter >= 30){ ## If 5 is to low, try 10 instead of 5.
      if(check_zero_stabilize(zero.betas, zero.alphas, iter)) break
    }
    
    ## 3. Sigma (TODO: revamp this.)
    ntlist = sapply(ylist, nrow)
    sigma = Mstep_sigma_flowcut(resp = resp,
                                ntlist = ntlist,
                                ylist = ylist_uncens,
                                new_responses = new_responses,
                                idx_cens_list = idx_cens_list,
                                idx_uncens_list = idx_uncens_list,
                                left_cens_list = left_cens_list,
                                right_cens_list = right_cens_list,
                                ntlist_orig = ntlist,
                                second_moment_list = second_moments, ## TODO; change argument name
                                mn = mn)

      ## 3. (Continue) Decompose the sigmas.
      sigma_eig_by_clust <- eigendecomp_sigma_array(sigma)
      denslist_by_clust <- make_denslist_eigen(ylist, mn, TT, dimdat, numclust,
                                               sigma_eig_by_clust)
      
      ## Calculate the objectives
      objectives[iter] = objective(mn, prob, sigma, ylist,
                                   prob_lambda = prob_lambda,
                                   mean_lambda = mean_lambda,
                                   alpha = alpha, beta = beta,
                                   denslist_by_clust = denslist_by_clust,
                                   countslist = countslist)
      

        
    ## Check convergence
    converged = check_converge_rel(objectives[iter-1],
                                   objectives[iter],
                                   tol = tol_em)
    if((iter > 10) & converged){
      break()
    }
    ## if(objectives[iter] > objectives[iter-1] * 1.01 ) break # Additional stopping
    ## of the likelihood
    ## increasing more
    ## than 1%.
  }
  
  ## Measure time
  lapsetime = difftime(Sys.time(), start.time, units = "secs")
  time_per_iter = lapsetime / (iter-1)
  
  
  ## Also calculate per-cytogram likelihoods (NOT divided by nt)
  loglikelihoods = objective(mn, prob, sigma, ylist,
                             prob_lambda = prob_lambda,
                             mean_lambda = mean_lambda,
                             alpha = alpha, beta = beta,
                             denslist_by_clust = denslist_by_clust,
                             countslist = countslist,
                             each = TRUE)
  ## loglikelihoods_particle = objective(mn, prob, sigma, ylist,
  ##                            prob_lambda = prob_lambda,
  ##                            mean_lambda = mean_lambda,
  ##                            alpha = alpha, beta = beta,
  ##                            denslist_by_clust = denslist_by_clust,
  ##                            countslist = countslist,
  ##                            each = FALSE,
  ##                            sep=TRUE)
  
  ## Also reformat the coefficients
  obj <- reformat_coef(alpha, beta, numclust, dimdat, X)
  alpha = obj$alpha
  beta = obj$beta
  
  return(structure(list(alpha = alpha,
                        beta = beta,
                        mn = mn,
                        prob = prob,
                        sigma = sigma,
                        ## denslist_by_clust = denslist_by_clust,
                        objectives = objectives[2:iter],
                        final.iter = iter,
                        time_per_iter = time_per_iter,
                        total_time = lapsetime,
                        loglikelihoods = loglikelihoods,
                        ## loglikelihoods_particle = loglikelihoods_particle,
                        ## Above is output, below are data/algorithm settings.
                        dimdat = dimdat,
                        TT = TT,
                        N = N,
                        p = p,
                        numclust = numclust,
                        X = X,
                        prob_lambda = prob_lambda,
                        mean_lambda = mean_lambda,
                        maxdev=maxdev,
                        niter = niter,
                        admm_niters = admm_niters
  ), class = "flowmix"))
  
  
  
}

## ## Some tests to add
## ## object is the result of having run flowmix() or flowmix_once().
## check_size <- function(obj){
##   assert_that(check_beta_size(res$beta, p, dimdat, numclust))
##   assert_that(check_alpha_size(res$alpha, p, dimdat))
## }

## check_beta_size <- function(beta, p, dimdat, numclust){
##   all.equal(dim(beta), c(p+1, dimdat, numclust))
## }
## check_alpha_size <- function(alpha, p, dimdat){
##   all.equal(dim(alpha), c(dimdat, p+1))
## }


#### NEED TO CHANGE THIS: DOING THIS NEXT...
##' Prediction: Given new covariates X's, generate a set of predicted cluster
##' means and probs (and return the same Sigma).
##'
##' @param object Object returned from \code{flowmix()}.
##' @param logits Logical: should the function return the logits of the cluster probs (TRUE)  
##'   or the cluster probs (FALSE)?
##' @param ... Make sure to provide new covariate values, a row vector in the
##'   same format and column names as the original X used in \code{object}, in
##'   \code{newx}.
##'
##' @return List containing mean, prob, and sigma. If 
##' user specifies logits = TRUE, prob will still be 
##' called prob, but it will actually contain the linear 
##' predictor of prob (the multinomial analogue of the 
##' logit function). 
##'
##' @export
##'
predict.flowmix <- function(object, logits = FALSE, ...){
  
  ## Basic checks
  ## stopifnot(ncol(new.x) == ncol(object$X))
  ## newx = X[1,,drop=FALSE]
  . = NULL ## Fixing check()
  rest = list(...)
  newx = rest$newx
  if(is.null(newx)){
    newx = object$X
  }
  
  ## Check if the variable names are the same.
  cnames = object$X %>% colnames()
  cnames_new = newx %>% colnames()
  stopifnot(all(cnames == cnames_new))
  
  ## Augment it with a dummy variable 1
  if(nrow(newx)>1){
    newx.a = cbind(rep(1, nrow(newx)), newx)
  } else {
    newx.a = c(1, newx)
  }
  
  TT = nrow(newx) ## This used to be nrow(X)..
  numclust = object$numclust
  dimdat = object$dimdat
  if(is.null(dimdat)) dimdat = object %>%.$mn %>% dim() %>% .[2] ## for back=compatibility
  
  ## Predict the means (manually).
  newmn = lapply(1:numclust, function(iclust){
    newx.a %*% object$beta[[iclust]]
  })
  newmn_array = array(NA, dim=c(TT, dimdat, numclust))
  for(iclust in 1:numclust){ newmn_array[,,iclust] = newmn[[iclust]] }
  
  ## Predict the probs.
  ## newprob = predict(object$alpha.fit, newx=newx, type='response')[,,1]
  probhatmat = as.matrix(tcrossprod(cbind(1, newx), object$alpha))
  if(logits) {
    newprob = probhatmat
  } else {
    probhatmat = exp(probhatmat)
    newprob = probhatmat / rowSums(probhatmat)
  }
  # probhatmat = as.matrix(exp(cbind(1,newx) %*% t(object$alpha)))
  # newprob = probhatmat / rowSums(probhatmat)
  
  ## predict(fit, newx=X, type="response")[,,1]
  stopifnot(all(dim(newprob) == c(TT,numclust)))
  stopifnot(logits | all(newprob >= 0)) # stop if not logits and any probs < 0 
  
  ## Return all three things
  return(list(mn = newmn_array,
              prob = newprob,
              pie = newprob, ## Just a copy of prob, for back-compatibility
              alpha = object$alpha,
              beta = object$beta,
              sigma = object$sigma,
              TT = object$TT,
              N = object$N,
              numclust = object$numclust,
              X = newx))
}



##' Helper for making list of densities. Returns list by cluster then time
##' e.g. access by \code{denslist_by_clust[[iclust]][[tt]]}
##'
##' @param ylist T-length list each containing response matrices of size (nt x
##'   3), which contains coordinates of the 3-variate particles, organized over
##'   time (T) and with (nt) particles at every time.
##' @param mu (T x dimdat x numclust) array.
##' @param dimdat dimension of data.
##' @param numclust number of clusters.
##' @param TT number of time points
##' @param sigma_eig_by_clust Result of running
##'   \code{eigendecomp_sigma_array(sigma.list[[iter]])}.
##'
##' @return numclust-lengthed list of TT-lengthed.
##'
##' @noRd
make_denslist_eigen <- function(ylist, mu,
                                TT, dimdat, numclust,
                                sigma_eig_by_clust){
  
  ## Basic checks
  assertthat::assert_that(!is.null(sigma_eig_by_clust))
  
  ## Calculate densities (note to self: nested for loop poses no problems)
  lapply(1:numclust, function(iclust){
    mysigma_eig <- sigma_eig_by_clust[[iclust]]
    lapply(1:TT, function(tt){
      mn = mu[tt,,iclust]
      sgm = mysigma_eig$sigma
      if(dimdat>1){
        ## return(dmvnorm_fast(ylist[[tt]],
        ##                     mu[tt,,iclust],
        ##                     sigma_eig=mysigma_eig))

        #if(dimdat == 1){
        #  mn = as.matrix(mn)
        #  sgm = sgm %>% as.matrix()
        #}
        return(dmvnorm_arma_fast(ylist[[tt]],
                                 mn,
                                 sgm))
      }else{
        return(dnorm(ylist[[tt]],
                                 mn,
                                 sgm))
      }

    })
  })
}



##' Functions to check convergence.
##' @noRd
check_converge_rel <- function(old, new, tol=1E-6){ return(abs((old-new)/old) < tol )  }




#' Convert the old format to the new.
#'
#' @param ylist original ylist
#' @param cens_indicator_list_left original left censor indicator
#' @param cens_indicator_list_right original right censor indicator.
#'
#' @return Separate list of responses, uncensored and censored. Also, two list
#'   objects that mark, in which (of the dimdat) dimension the censoring has
#'   happened, left or right.
my_convert <- function(ylist,
                    cens_indicator_list_left,
                    cens_indicator_list_right){

  ## Basic setup
  TT <- length(ylist)

  ## Make a cens_left_list and cens_right_list the same size as ylist, but with
  ## 0/1s.
  cens_left_list =
    lapply(cens_indicator_list_left, function(one_mat){
      one_mat[is.na(one_mat)] <- 0
      return(one_mat)
    })

  cens_right_list =
    lapply(cens_indicator_list_right, function(one_mat){
      one_mat[is.na(one_mat)] <- 0
      return(one_mat)
    })

  ## Do the conversion (from ylist to ylist_uncen and ylist_cen)
  ylist_cens = ylist_uncens = list()
  left_cens_list = right_cens_list = list()
  idx_cens_list = idx_uncens_list = list() ## track rows

  for(tt in 1:TT){
    one_y_full = ylist[[tt]]
    nt = nrow(one_y_full)
    one_left_bool = cens_left_list[[tt]]
    one_right_bool = cens_right_list[[tt]]

    left_cens_rows = which(rowSums(one_left_bool) > 0)
    right_cens_rows = which(rowSums(one_right_bool) > 0)

    # Unique censored row numbers (to avoid duplicates if a row is left AND right censored)
    all_censored_rownums = unique(c(left_cens_rows, right_cens_rows))
    all_uncensored_rownums = setdiff(1:nt, all_censored_rownums)

    one_y_cens = one_y_full[all_censored_rownums, , drop=FALSE]
    one_y_uncens = one_y_full[all_uncensored_rownums, , drop=FALSE]

    ylist_cens[[tt]] = one_y_cens
    ylist_uncens[[tt]] = one_y_uncens

    # Store the index tracking
    idx_cens_list[[tt]] = all_censored_rownums
    idx_uncens_list[[tt]] = all_uncensored_rownums

    left_cens_list[[tt]] = one_left_bool[all_censored_rownums, , drop=FALSE]
    right_cens_list[[tt]] = one_right_bool[all_censored_rownums, , drop=FALSE]
  }

  return(list(ylist_cens = ylist_cens,
              ylist_uncens = ylist_uncens,
              left_cens_list = left_cens_list,
              right_cens_list = right_cens_list,
              idx_cens_list = idx_cens_list,
              idx_uncens_list = idx_uncens_list))
}


#' An opposite operation from the separation in myconvert().
#'
my_reassemble <- function(## converted_res
                          ylist_uncens,
                          ylist_cens,
                          idx_cens_list,
                          idx_uncens_list
                          ) {
  ## ylist_cens = converted_res$ylist_cens
  ## ylist_uncens = converted_res$ylist_uncens
  ## idx_cens_list = converted_res$idx_cens_list
  ## idx_uncens_list = converted_res$idx_uncens_list

  TT <- length(ylist_cens)
  ylist_recovered <- list()

  for(tt in 1:TT) {

    idx_cens <- idx_cens_list[[tt]]
    idx_uncens <- idx_uncens_list[[tt]]

    # Calculate total rows and columns for this time step
    total_rows <- length(idx_cens) + length(idx_uncens)
    total_cols <- ncol(ylist_cens[[tt]])

    # Create an empty matrix of the original size
    recovered_mat <- matrix(NA, nrow = total_rows, ncol = total_cols)

    # Map the rows back using their saved index positions
    if(length(idx_cens) > 0) {
      recovered_mat[idx_cens, ] <- ylist_cens[[tt]]
    }
    if(length(idx_uncens) > 0) {
      recovered_mat[idx_uncens, ] <- ylist_uncens[[tt]]
    }

    # Preserve original column names if they existed
    colnames(recovered_mat) <- colnames(ylist_cens[[tt]])

    ylist_recovered[[tt]] <- recovered_mat
  }

  return(ylist_recovered)
}
