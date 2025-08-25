library(coda)
library(dplyr)

summarize <- function(results, n_iter, burn_in, rep, ind = FALSE){
  
  if (ind == FALSE) {
    ind <- 1:length(results)
  }
  
  N_mean <- c() 
  N_mode <- c()
  N_sd <- c()
  N_median <- c()
  N_CI_95_lower <- c()
  N_CI_95_upper <- c()
  N_ESS <- c()
  N_ESS_per_s <- c()
  
  sigma_mean <- c()
  sigma_mode <- c()
  sigma_sd <- c()
  sigma_median <- c()
  sigma_CI_95_lower <- c()
  sigma_CI_95_upper <- c()
  sigma_ESS <- c()
  sigma_ESS_per_s <- c()
  
  lam0_mean <- c()
  lam0_mode <- c()
  lam0_sd <- c()
  lam0_median <- c()
  lam0_CI_95_lower <- c()
  lam0_CI_95_upper <- c()
  lam0_ESS <- c()
  lam0_ESS_per_s <- c()
  
  
  for (i in ind){
    
    mc_samples <- results[[i]]$out
    sigma_chain <- mc_samples[burn_in:n_iter,1] #sigma 
    lam0_chain <- mc_samples[burn_in:n_iter,2] #lambda
    N_chain <- mc_samples[burn_in:n_iter,3] #N
    
    #append list
    N_mean <- append(N_mean, mean(N_chain))
    N_mode <- append(N_mode, as.numeric(names(table(N_chain))[which.max(table(N_chain))]))
    N_sd <- append(N_sd, sd(N_chain))
    N_median <- append(N_median, median(N_chain))
    N_CI_95_lower <- append(N_CI_95_lower, quantile(N_chain, 0.025))
    N_CI_95_upper <- append(N_CI_95_upper, quantile(N_chain, 0.975))
    N_ESS <- append(N_ESS, effectiveSize(N_chain))
    N_ESS_per_s <- append(N_ESS_per_s, effectiveSize(N_chain)/(n_iter-burn_in))
    
    sigma_mean <- append(sigma_mean, mean(sigma_chain))
    sigma_mode <- append(sigma_mode, density(sigma_chain)$x[which.max(density(sigma_chain)$y)])
    sigma_sd <- append(sigma_sd, sd(sigma_chain))
    sigma_median <- append(sigma_median, median(sigma_chain))
    sigma_CI_95_lower <- append(sigma_CI_95_lower, quantile(sigma_chain, 0.025))
    sigma_CI_95_upper <- append(sigma_CI_95_upper, quantile(sigma_chain, 0.975))
    sigma_ESS <- append(sigma_ESS, effectiveSize(sigma_chain))
    sigma_ESS_per_s <- append(sigma_ESS_per_s, effectiveSize(sigma_chain)/(n_iter-burn_in))
    
    lam0_mean <- append(lam0_mean, mean(lam0_chain))
    lam0_mode <- append(lam0_mode, density(lam0_chain)$x[which.max(density(lam0_chain)$y)])
    lam0_sd <- append(lam0_sd, sd(lam0_chain))
    lam0_median <- append(lam0_median, median(lam0_chain))
    lam0_CI_95_lower <- append(lam0_CI_95_lower, quantile(lam0_chain, 0.025))
    lam0_CI_95_upper <- append(lam0_CI_95_upper, quantile(lam0_chain, 0.975))
    lam0_ESS <- append(lam0_ESS, effectiveSize(lam0_chain))
    lam0_ESS_per_s <- append(lam0_ESS_per_s, effectiveSize(lam0_chain)/(n_iter-burn_in))
    
  }
  
  
  result_summary <- data.frame(
  N_mean = N_mean,
  N_mode = N_mode,
  N_sd = N_sd,
  N_median = N_median,
  N_CI_95_lower = N_CI_95_lower,
  N_CI_95_upper = N_CI_95_upper,
  N_ESS = N_ESS,
  N_ESS_per_s = N_ESS_per_s,
  sigma_mean = sigma_mean,
  sigma_mode = sigma_mode,
  sigma_sd = sigma_sd,
  sigma_median = sigma_median,
  sigma_CI_95_lower = sigma_CI_95_lower,
  sigma_CI_95_upper =  sigma_CI_95_upper,
  sigma_ESS = sigma_ESS,
  sigma_ESS_per_s = sigma_ESS_per_s,
  lam0_mean = lam0_mean,
  lam0_mode = lam0_mode,
  lam0_sd= lam0_sd,
  lam0_median = lam0_median,
  lam0_CI_95_lower = lam0_CI_95_lower,
  lam0_CI_95_upper = lam0_CI_95_upper,
  lam0_ESS= lam0_ESS,
  lam0_ESS_per_s = lam0_ESS_per_s
  )
  
  summary <- result_summary %>% mutate(group = (1:dim(result_summary)[1]-1) %/% rep) %>% group_by(group) %>% summarise(across(everything(), mean)) 
  summary <- summary %>% select(N_mean ,N_mode,  N_sd, N_median, N_CI_95_lower, N_CI_95_upper, sigma_mean, sigma_mode,  sigma_sd, sigma_median, sigma_CI_95_lower, sigma_CI_95_upper,lam0_mean, lam0_mode, lam0_sd, lam0_median, lam0_CI_95_lower, lam0_CI_95_upper)
  return(summary)
}
  
