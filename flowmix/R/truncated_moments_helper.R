###################

#' Computes the normal probabilities giving a small lower bound to very low probability events
#' 
#' @param lower d-dimensional vector of lower entries. Can be thought of as the lowest edge in the cube considered.
#' @param upper d-dimensional vector of upper entries. Can be thought of as the highest edge in the cube considered.
#' @param mean the d-dimensional vector of means for the multivariate normal.
#' @param sigma the dxd covariance matrix for the multivariate normal.
#' @return a probability value. This value is lower bounded by 1e-12.

my_pmvnorm <- function(lower, upper, mean, sigma){
  prob <- tryCatch({
    mvtnorm::pmvnorm(lower, upper, mean = mean, sigma = sigma)
  }, error = function(e) {
    NaN
  })
  if (prob < 1e-12 || is.nan(prob)) {
    return(1e-12)
  }else{
    return(mvtnorm::pmvnorm(lower,upper, mean = mean, sigma = sigma))
  }
  
}



##################


#' Compute conditional (lower-dimensional) multivariate normal (TODO: if needed, y can turned into a matrix).
#'
#' @param y d-dimensional observation that may be censored in some dimensions.
#' @param mu d-dimensionl Normal mean.
#' @param Sigma Normal dxd covariance matrix.
#' @param observed_dims A subset of 1,2,..,d on which we condition
#'
#' @return list containing mean and variance of the conditional multivariate
#'   normal.


cond_mean_var_func <- function(y,mu,Sigma,observed_dims){
  
  ## Basic check
  dimdat = length(y)
  stopifnot(length(observed_dims) == 0 || all(observed_dims %in% seq_len(dimdat)))
  stopifnot(length(mu)==dimdat)
  stopifnot(all(dim(Sigma)==dimdat))
  
  ##Main Body
  y_observed_dims = y[observed_dims]
  
  muu = mu[-observed_dims]
  muo = mu[observed_dims]
  
  Sigmaoo = Sigma[observed_dims,observed_dims]
  Sigmauu = Sigma[-observed_dims,-observed_dims]
  Sigmauo = Sigma[-observed_dims,observed_dims]
  Sigmaou = Sigma[observed_dims,-observed_dims]
  
  
  
  if(any(y_observed_dims==-Inf)|any(y_observed_dims==Inf)){
    mu_out = muu
    Sigma_out = Sigmauu
  }else{
    
    if((length(observed_dims)==length(mu))||(length(observed_dims)== 0)){
      mu_out=mu
      Sigma_out=Sigma
    }
    if(length(observed_dims)==1){
      Sigma_out = Sigmauu - (1/Sigmaoo) * Sigmauo %*% t(Sigmaou)
      mu_out = muu+Sigmauo * (y_observed_dims-muo)/Sigmaoo
    }
    
    
    if((length(observed_dims)>1)&(length(observed_dims)<length(mu))){
      
      Sigmaoo_svd = svd(Sigmaoo)
      
      Uoo = Sigmaoo_svd$u
      Doo = Sigmaoo_svd$d
      Doo_inv = diag(1/(Sigmaoo_svd$d))
      
      Sigmaooinv = Uoo %*% Doo_inv %*% t(Uoo)
      
      
      mu_out = muu+Sigmauo %*% Sigmaooinv %*% (y_observed_dims-muo)
      Sigma_out = Sigmauu - Sigmauo %*% Sigmaooinv %*% Sigmaou
    }
    
  }     
  
  
  return(list(mu_conditional=mu_out,
              Sigma_conditional=Sigma_out))
  
}




###################

#'Compute the moments of a d-dimensional truncated Gaussian with zero mean with given truncation limits.
#'@param a_vec the lower truncation limits; is a vector of length d.
#'@param b_vec the upper truncation limits; is a vector of length d.
#'@param Sigma the covariance matrix for the normal of size dxd.
#'
#'@return a list containing the first and second moments.



moment_cal_func_centered <- function(a_vec,b_vec,Sigma){
  
  #browser()
  ## Basic checks
  dimdat=length(a_vec) 
  stopifnot(all(dim(Sigma)==dimdat))
  stopifnot(length(a_vec)==length(b_vec))
  
  
  prob_andd_density_function_F <- function(x,a_vec,b_vec,pos,Sigma){ #F_k(x) in Lee 2012
    if(length(a_vec)>=2){
      a_vec_min_pos = a_vec[-pos]
      b_vec_min_pos = b_vec[-pos]
      
      x_vec = a_vec
      x_vec[pos] = x
      
      sigma = Sigma[pos,pos]
      Sigma_cond_pos = cond_mean_var_func(x_vec, rep(0,length(x_vec)), Sigma, c(pos))$Sigma_conditional
      mu_cond_pos = cond_mean_var_func(x_vec, rep(0,length(x_vec)), Sigma, c(pos))$mu_conditional
      
      out = dnorm(x,0,sqrt(sigma)) * my_pmvnorm(a_vec_min_pos,b_vec_min_pos, mean = as.vector(mu_cond_pos), sigma = Sigma_cond_pos)
    }else{
      sigma = as.numeric(Sigma)
      out = dnorm(x,0,sqrt(sigma))
    }
    
    
    return(out) 
  }
  
  
  
  #this function calculates the density/probabilities required to calculate the truncated second moment 
  prob_andd_density_function_FF <- function(x, y, a_vec,b_vec, pos1, pos2, Sigma){ # F_{k,q}(x,y) as in Lee 2012
    if(length(a_vec)>=3){
      if(pos1==pos2){
        print("positions cannot be same")
      }else{
        pos=c(pos1,pos2)
        a_vec_min_pos = a_vec[-pos]
        b_vec_min_pos = b_vec[-pos]
        
        x_vec = a_vec
        
        x_vec[pos] = c(x,y)
        
        sigma = Sigma[pos,pos]
        
        Sigma_cond_pos = cond_mean_var_func(x_vec, rep(0,length(x_vec)), Sigma, c(pos))$Sigma_conditional
        mu_cond_pos = cond_mean_var_func(x_vec, rep(0,length(x_vec)), Sigma, c(pos))$mu_conditional
        
        if(length(a_vec_min_pos)==1){
          if(any(c(x,y)==-Inf) || any(c(x,y)==Inf)){
            out=0
          }else{
            out = mvtnorm::dmvnorm(c(x,y),rep(0,2),sigma) * (pnorm(b_vec_min_pos, mu_cond_pos, sqrt(Sigma_cond_pos))-pnorm(a_vec_min_pos, mu_cond_pos, sqrt(Sigma_cond_pos)))
          }
          
        }else{
          if(any(c(x,y)==-Inf) || any(c(x,y)==Inf) || any(a_vec_min_pos ==-Inf) || any(a_vec_min_pos ==Inf) || any(b_vec_min_pos ==-Inf) || any(b_vec_min_pos ==Inf)){
            out=0
          }else{
            out = mvtnorm::dmvnorm(c(x,y),rep(0,2),sigma) * my_pmvnorm(a_vec_min_pos,b_vec_min_pos, as.vector(mu_cond_pos), Sigma_cond_pos)
          }
          
        }
      }
      
    }else{
      out=0
    }  
    
    
    return(out) 
  }
  
  
  #debug(prob_andd_density_function_F )
  #debug(prob_andd_density_function_FF)
  
  alpha = my_pmvnorm(a_vec, b_vec, mean = rep(0,length(a_vec)), sigma = Sigma)
  mean_vec = numeric(dimdat)
  Sigma_raw_mat= matrix(0, dimdat,dimdat)
  
  temp_vec_1 = numeric(dimdat)
  temp_vec_2 = numeric(dimdat)
  temp_vec_3 = numeric(dimdat)
  temp_vec_4 = numeric(dimdat)
  
  for(i in 1:dimdat){
    #        debug(prob_andd_density_function_F)
    temp_vec_1[i] = prob_andd_density_function_F(a_vec[i],a_vec,b_vec,i,Sigma)
    temp_vec_2[i] = prob_andd_density_function_F(b_vec[i],a_vec,b_vec,i,Sigma)
    
    if(a_vec[i] == -Inf || a_vec[i] == Inf) {
      temp_vec_3[i] = 0
    } else {
      temp_vec_3[i] = a_vec[i] * prob_andd_density_function_F(a_vec[i], a_vec, b_vec, i, Sigma)
    }
    
    if(b_vec[i] == -Inf || b_vec[i] == Inf) {
      temp_vec_4[i] = 0
    } else {
      temp_vec_4[i] = b_vec[i] * prob_andd_density_function_F(b_vec[i], a_vec, b_vec, i, Sigma)
    }
    
  }
  
  for(i in 1:dimdat){
    mean_vec[i] = Sigma[i,] %*% (temp_vec_1-temp_vec_2)/alpha
  }
  
  #      debug(prob_andd_density_function_FF)
  temp_mat_1 = matrix(0,dimdat,dimdat)
  temp_mat_2 = matrix(0,dimdat,dimdat)
  for(k in 1:dimdat){
    for(q in 1:dimdat){
      if(k==q){
        temp_mat_1[k,q]=0
        temp_mat_2[k,q]=0
      }else{
        temp_mat_1[k,q] = prob_andd_density_function_FF(a_vec[k], a_vec[q], a_vec, b_vec,k, q, Sigma) + 
          prob_andd_density_function_FF(b_vec[k], b_vec[q], a_vec, b_vec, k, q, Sigma)
        
        temp_mat_2[k,q] = prob_andd_density_function_FF(a_vec[k], b_vec[q] , a_vec,b_vec, k, q, Sigma) + 
          prob_andd_density_function_FF(b_vec[k], a_vec[q] , a_vec,b_vec, k, q, Sigma)
      }
      
      
    }
    
  }
  
  
  temp_fn <- function(i,j,dimdat,temp_mat_1,temp_mat_2,Sigma){
    out_k=0
    for(k in 1:dimdat){
      out_q=0
      for(q in 1:dimdat){
        if(q==k){
          out_q=out_q
        }else{
          out_q=out_q+(Sigma[j,q]-(Sigma[k,q] * Sigma[j,k])/Sigma[k,k]) * (temp_mat_1[k,q]- temp_mat_2[k,q]) 
        }
      }
      out_k =out_k +Sigma[i,k] *out_q
    }
    return(out_k)
  }
  
  Sigma_diag_inv = matrix(0,dimdat,dimdat)
  
  diag(Sigma_diag_inv) = (temp_vec_3-temp_vec_4)/diag(Sigma)
  
  for(i in 1:dimdat){
    for(j in 1:dimdat){
      Sigma_raw_mat[i,j] = Sigma[i,j] + Sigma[i,]%*% Sigma_diag_inv %*% Sigma[,j]/alpha +
        temp_fn(i,j,dimdat,temp_mat_1,temp_mat_2,Sigma)/alpha
    }
    
  }
  
  for(i in 1:dimdat){
    for(j in 1:dimdat){
      if ((!is.na(Sigma_raw_mat[i,j])) && 
          (!is.na(Sigma_raw_mat[j,i])) && 
          ((abs(Sigma_raw_mat[i,j] - Sigma_raw_mat[j,i]) / 
            ifelse(abs(Sigma_raw_mat[i,j] + Sigma_raw_mat[j,i]) == 0, 1, abs(Sigma_raw_mat[i,j] + Sigma_raw_mat[j,i]))) < 1e-6)){
        Sigma_raw_mat[i,j] = Sigma_raw_mat[j,i]
      }else{
        print(Sigma_raw_mat)
        #browser()
        stop("Censored second moment matrix is not symmetric")
      }
    }
  } 
  
  for(i in 1:dimdat){
    if(Sigma_raw_mat[i,i]<=0){
      print(a_vec)
      print(b_vec)
      print(Sigma)
      print(Sigma_raw_mat)
      #browser()
      stop("Variances are negative")
    }
  }
  
  #stopifnot(all(Sigma_raw_mat == t(Sigma_raw_mat)))
  
  return(list(mean_vec=mean_vec,
              Sigma_raw_mat=Sigma_raw_mat))  
  
}


#'Compute the moments of a d-dimensional truncated Gaussian with zero mean with given truncation limits.
#'@param ii the particle index for which the censored 
#'@param tt the time point in question
#'@param y the response. Is a list with a matrix at each time point of dimension nt x dimension.
#'@param X the matrix of covariates with dimensions tt x d.
#'@param censor_indicator_left A list checking whether an observation is left censored. The list if of size TT with each entry being a ntxdim  matrix.
#'@param censor_indicator_right A list checking whether an observation is right censored. The list if of size TT with each entry being a ntxdim  matrix.
#'@param cens_lim_l_vec Lower censoring limits.
#'@param cens_lim_u_vec Upper censoring limits.
#'@param numclust number of clusters
#'@param mu_list list of means by cluster. Each cluster has a matrix of dim TT x d.
#'@param sigma_list list of covariance matrices by cluster. Each cluster has a matrix of dim d x d.
#'
#'@return a list for the means and variances of the censored/uncensored observation ii at time tt per each cluster.
#'
#'

cens_cond_normal <- function(ii,tt,y,X,censor_indicator_left, 
                             censor_indicator_right,
                             cens_lim_l_vec,
                             cens_lim_u_vec,
                             numclust, mu_list,sigma_list){
  
  #print(c(tt,ii))
  #if((tt==3)&(ii==71)) browser()
  
  dimdat = ncol(y[[1]])
  ntlist = sapply(y, nrow)
  
  ##Basic Checks
  stopifnot(is.list(y))
  stopifnot(is.list(censor_indicator_left))
  stopifnot(is.list(censor_indicator_right))
  stopifnot(is.list(mu_list))
  stopifnot(is.list(sigma_list))
  stopifnot(length(cens_lim_l_vec)==length(cens_lim_u_vec))
  
  ##Main body
  
  ##if((debug_global==9)&&(tt==3)&&(ii==71)) debug(moment_cal_func_centered)
  
  mu_conditional = list()
  Sigma_conditional = list()
  prob_cens_conditional = list()
  conditional_response = list()
  all_cond_mean_list = list()
  all_cond_second_moment_list = list()
  new_response_list = list()
  y_vec = y[[tt]][ii,]
  
  left_cens_index = which(censor_indicator_left[[tt]][ii,]==1)
  right_cens_index = which(censor_indicator_right[[tt]][ii,]==1)
  uncensored_index=which((censor_indicator_left[[tt]][ii,]!=1 |is.na(censor_indicator_left[[tt]][ii,]))&(censor_indicator_right[[tt]][ii,]!=1|is.na(censor_indicator_right[[tt]][ii,])))
    lower_limits = cens_lim_l_vec[left_cens_index]
    upper_limits = cens_lim_u_vec[right_cens_index]
    for(iclust in 1:numclust){
      #print(iclust)
      #if(iclust==2) debug(moment_cal_func_centered)
      mu = mu_list[[iclust]][tt,]
      Sigma = sigma_list[[iclust]]
      
      mu_conditional[[iclust]] = cond_mean_var_func(y_vec, mu,Sigma,uncensored_index)$mu_conditional
      Sigma_conditional[[iclust]] = cond_mean_var_func(y_vec, mu,Sigma,uncensored_index)$Sigma_conditional
      
      p_lower_limit = numeric(dimdat)  
      p_upper_limit = numeric(dimdat)
      p_lower_limit[left_cens_index] = -Inf
      p_lower_limit[uncensored_index] = -Inf
      p_lower_limit[right_cens_index] = upper_limits
      p_upper_limit[left_cens_index] = lower_limits
      p_upper_limit[uncensored_index] = Inf
      p_upper_limit[right_cens_index] = Inf
      
      
      
      if((length(upper_limits)==0)&(length(lower_limits)==0)){
        prob_cens_conditional[[iclust]] = NA
        all_cond_second_moment_list[[iclust]] = NA
        all_cond_mean_list[[iclust]] = NA

      }  
    
      
      if((length(upper_limits)>0)&(length(lower_limits)>0)){
        if(length(uncensored_index)==0){
          prob_cens_conditional[[iclust]] = my_pmvnorm(p_lower_limit[sort(c(left_cens_index,right_cens_index))], p_upper_limit[sort(c(left_cens_index,right_cens_index))], mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]])[1]
          
          astar_vec = as.vector(p_lower_limit[sort(c(left_cens_index,right_cens_index))] -mu_conditional[[iclust]])
          bstar_vec = as.vector(p_upper_limit[sort(c(left_cens_index,right_cens_index))] -mu_conditional[[iclust]])
          
          saved_moments<-moment_cal_func_centered(astar_vec,bstar_vec,Sigma_conditional[[iclust]])
          #print(saved_moments)
          saved_mean <- saved_moments$mean_vec
          #print(saved_mean)
          saved_second_moment_raw <- saved_moments$Sigma_raw_mat
          #print(saved_second_moment_raw)
          
          
          all_cond_mean_list[[iclust]] = mu_conditional[[iclust]] + saved_mean
          
          
          all_cond_second_moment_list[[iclust]] = mu_conditional[[iclust]]%*% t(mu_conditional[[iclust]])+ mu_conditional[[iclust]]%*% t(saved_mean)+ saved_mean %*% t(mu_conditional[[iclust]]) +saved_second_moment_raw
        }else{
          prob_cens_conditional[[iclust]] = my_pmvnorm(p_lower_limit[-uncensored_index], p_upper_limit[-uncensored_index], mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]])[1]
          
          astar_vec = as.vector(p_lower_limit[-uncensored_index] -mu_conditional[[iclust]])
          bstar_vec = as.vector(p_upper_limit[-uncensored_index] -mu_conditional[[iclust]])
          
          saved_moments<-moment_cal_func_centered(astar_vec,bstar_vec,Sigma_conditional[[iclust]])
          
          saved_mean <- saved_moments$mean_vec
          saved_second_moment_raw <- saved_moments$Sigma_raw_mat
          
          all_cond_mean_list[[iclust]] = mu_conditional[[iclust]] +saved_mean
          
          
          all_cond_second_moment_list[[iclust]] = mu_conditional[[iclust]]%*% t(mu_conditional[[iclust]])+ mu_conditional[[iclust]]%*% t(saved_mean)+ saved_mean %*% t(mu_conditional[[iclust]]) +saved_second_moment_raw
        }
      }            
      
      
      if((length(upper_limits)==0)&(length(lower_limits)>0)){
        prob_cens_conditional[[iclust]] = my_pmvnorm(rep(-Inf,length(lower_limits)), lower_limits, mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]])[1]
        astar_vec = rep(-Inf,length(lower_limits))
        bstar_vec = as.vector(lower_limits -mu_conditional[[iclust]])
        
        saved_moments<-moment_cal_func_centered(astar_vec,bstar_vec,Sigma_conditional[[iclust]])
        
        saved_mean <- saved_moments$mean_vec
        saved_second_moment_raw <- saved_moments$Sigma_raw_mat
        
        all_cond_mean_list[[iclust]] = mu_conditional[[iclust]]+saved_mean
        
        all_cond_second_moment_list[[iclust]] = mu_conditional[[iclust]]%*% t(mu_conditional[[iclust]])+ mu_conditional[[iclust]]%*% t(saved_mean)+
          saved_mean %*% t(mu_conditional[[iclust]]) +
          saved_second_moment_raw
        
      }
      
      if((length(upper_limits)>0)&(length(lower_limits)==0)){
        prob_cens_conditional[[iclust]] = my_pmvnorm(upper_limits,rep(Inf,length(upper_limits)), mean = as.vector(mu_conditional[[iclust]]), sigma = Sigma_conditional[[iclust]])
        
        astar_vec = as.vector(upper_limits -mu_conditional[[iclust]])
        bstar_vec = rep(Inf,length(upper_limits))
        
        saved_moments<-moment_cal_func_centered(astar_vec,bstar_vec,Sigma_conditional[[iclust]])
        
        saved_mean <- saved_moments$mean_vec
        saved_second_moment_raw <- saved_moments$Sigma_raw_mat
        
        all_cond_mean_list[[iclust]] = mu_conditional[[iclust]]+saved_mean
        
        all_cond_second_moment_list[[iclust]] = mu_conditional[[iclust]]%*% t(mu_conditional[[iclust]])+ mu_conditional[[iclust]]%*% t(saved_mean)+ saved_mean %*% t(mu_conditional[[iclust]]) + saved_second_moment_raw
        
      }
      new_response_list[[iclust]] = numeric(dimdat)
      new_response_list[[iclust]][uncensored_index] = y_vec[uncensored_index]
      new_response_list[[iclust]][-uncensored_index] = all_cond_mean_list[[iclust]]
    }
    


  #  lapply(1:numclust, function(iclust) {
  #    mat <- all_cond_second_moment_list[[iclust]]
  
  # Skip IF any NA
  #    if (!any(is.na(mat))) {
  #      if (!all(mat == t(mat))) {
  #        print(paste("Cluster number, time, particle number in order is :", 
  #                    paste(c(iclust, tt, ii), collapse = " ")))
  #        stop("Non-symmetric second moment found. Stopping execution.")
  #      }
  #    }
  #  })
  
  return(list(new_response_list = new_response_list,
              all_conditional_means = all_cond_mean_list,
              all_conditional_second_moment_list = all_cond_second_moment_list))
}







  