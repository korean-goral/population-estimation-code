#===========================================================================================================
# implement RJMCMC <- fixed proposal distribution method

rjmcmc <- function(n, X, N_prior, sigma_prior, niters, xlims, ylims, thres = 1, burn_in = 10000, init_delta = 10, init_tune=c(0.1, 0.1, 2), N_updates = 15, log, stochastic = FALSE, lam0_prior = FALSE, check_point = FALSE, monitorS=FALSE)
{
  
  start <- proc.time()
  
  # initial setting
  K <- ncol(n) # number of periods
  sigma_ups <- 0
  lam0_ups <- 0
  N_ups <- 0
  
  
  
  sigma_tune <- init_tune[1]
  lam0_tune <- init_tune[2]
  S_tune <- init_tune[3]
  delta <- init_delta
  
  
  
  sigma_acc_count <- 0
  lam0_acc_count <- 0
  S_acc_count <- 0
  S_try_count <- 0
  N_acc_count <- 0
  N_try_count <- 0
  
  # set initial value : start with out checkpoints
  if (isFALSE(check_point)){ 
    N <- sample(1:500 , 1)
    sigma <-runif(1, 0, 5)
    lam0 <- runif(1, .1, 1)
    ## generate distance matrix
    
    S <- cbind(runif(N, xlims[1], xlims[2]),
               runif(N, ylims[1], ylims[2]))
    init_ind <- 1
  }
  
  else {
    N <- check_point$last$N
    sigma <- check_point$last$sigma
    lam0 <- check_point$last$lam0
    S <- check_point$last$S
    init_ind <- dim(check_point$out)[1]
  }
  
  
  D <- e2dist1(S, X)
  lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
  lam <- (lam * (D < thres))
  
  ## Big lambda (parameter for each camera)
  lamv.curr <- colSums(lam)
  
  # sampler
  ## just in case the first sigma is rejected
  ll <- sum(dpois(n, lamv.curr, log=TRUE)) # ll : log-likelihood
  
  ## matrix to hold samples
  out <- matrix(NA, nrow=niters, ncol=3)
  colnames(out) <- c("sigma", "lam0", "N")
  
  Sout <- NULL
  if(monitorS)
    Sout <- array(NA, c(N, 2, niters))
  
  log_file <- paste0("results/logs/", log) #주소 본인 코드에 맞게 바꾸기.
  sink(log_file, append = TRUE)
  
  cat("spNrj : sigma ~ unif, N ~ Negbin(),  lambda ~ unif, s ~ unif", "\n") # model description.
  cat("\ninitial values =", c(sigma, lam0, N), "\n\n")  
  
  
  #start iteration
  for(iter in 1:niters) {
    ## adaptive tuning and progress notice
    
    if(iter %% 100 == 0){
      sigma_acc_rate <- sigma_acc_count / 100
      lam0_acc_rate <- lam0_acc_count / 100
      N_acc_rate <- N_acc_count / N_try_count
      S_acc_rate <- S_acc_count / S_try_count
      
      
      ## print acceptance rate
      if(iter %% 1000 == 0) {
        cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
        cat("current =", out[iter-1,], "\n")
        cat("  Acceptance rates\n")
        cat("    S =", S_acc_rate, "\n") # ?
        cat("    N =", N_acc_rate, "\n") # ?
        cat("    lam0 =", lam0_acc_rate, "\n") # ?
        cat("    sigma =", sigma_acc_rate, "\n") # ?
      }
      
      # tune parameter
      if(iter <= burn_in){
        if(sigma_acc_rate > 0.4) sigma_tune <- sigma_tune * 1.1
        if(sigma_acc_rate < 0.2) sigma_tune <- sigma_tune * 0.9
        if(lam0_acc_rate > 0.4) lam0_tune <- lam0_tune * 1.1
        if(lam0_acc_rate < 0.2) lam0_tune <- lam0_tune * 0.9
        if(N_acc_rate > 0.4) delta <- delta + 1
        if(N_acc_rate < 0.2) delta <- max(1, delta - 1)
        if(S_acc_rate > 0.4) S_tune <- S_tune * 1.1
        if(S_acc_rate < 0.2) S_tune <- S_tune * 0.9
      }
      
      sigma_acc_count <- 0
      lam0_acc_count <- 0
      S_acc_count <- 0
      S_try_count <- 0
      N_acc_count <- 0
      N_try_count <- 0
    }
    ## update sigma
    sigma.cand <- rnorm(1, sigma, sigma_tune)
    if(sigma.cand > 1e-4) {
      lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand)) #if prior for sigma is uniform
      lam.cand <- (lam.cand * (D < thres))
      lamv.cand <- colSums(lam.cand)
      lamv.curr <- colSums(lam)
      ll<- sum(dpois(n, lamv.curr, log=TRUE))
      llcand<- sum(dpois(n, lamv.cand, log=TRUE))
      
      # prior for sigma : uniform
      if (sigma_prior[1] == 'uniform'){
        if(runif(1) < exp(llcand  - ll)){
          ll <- llcand
          lamv.curr <- lamv.cand
          lam <- lam.cand
          sigma <- sigma.cand
          sigma_acc_count <- sigma_acc_count + 1
        }
      }
      # prior for sigma : gamma
      if (sigma_prior[1] == 'gamma'){
        if(runif(1) < exp(llcand  - ll + dgamma(sigma.cand, as.double(sigma_prior[2]), 1/as.double(sigma_prior[3]), log = TRUE) - dgamma(sigma, as.double(sigma_prior[2]), 1/as.double(sigma_prior[3]), log = TRUE))){
          ll <- llcand
          lamv.curr <- lamv.cand
          lam <- lam.cand
          sigma <- sigma.cand
          sigma_acc_count <- sigma_acc_count + 1
        }
      }
    }
    
    
    # update lam0
    lam0.cand <- rnorm(1, lam0, lam0_tune)
    if((lam0.cand>0) & (lam0.cand < 1)) {
      lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
      lam.cand <- (lam.cand * (D < thres))
      lamv.cand <- colSums(lam.cand)
      llcand <- sum(dpois(n, lamv.cand, log=TRUE))
      
      prior <- dbeta(lam0, lam0_prior[1], lam0_prior[2], log=TRUE)
      prior.cand <- dbeta(lam0.cand, lam0_prior[1], lam0_prior[2], log=TRUE)
      
      if(runif(1) < exp((llcand + prior.cand) - (ll + prior))) { 
        ll <- llcand
        lamv.curr<-lamv.cand
        lam0<-lam0.cand
        lam<-lam.cand
        lam0_acc_count <- lam0_acc_count + 1
      }
    }
    
    # update S
    Sups <- 0
    for(i in 1:N) {
      S_try_count <- S_try_count + 1
      Scand <-c(rnorm(1, S[i,1], S_tune), rnorm(1, S[i,2], S_tune)) # S_i candidate
      inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] & Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
      if(!inbox)
        next
      
      dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2) # S_i ~ X_j 거리 : J차원 벡터
      lam.cand <- lam # IxJ matrix
      lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma)) # lambda_cand i th row update
      lam.cand[i,] <- lam.cand[i,]*(dtmp < thres)
      lamv.cand <- colSums(lam.cand) # sum each individual per camera
      llcand <- sum(dpois(n, lamv.cand, log=TRUE))
      
      if(runif(1)< exp(llcand - ll)) {
        ll <- llcand
        lamv.curr <- lamv.cand
        S[i,] <- Scand
        lam <- lam.cand
        D[i,] <- dtmp
        S_acc_count <- S_acc_count + 1
      }
    }
    
    # update N -> fixed or stochastic version
    eps <- c(c(-delta:-1), c(1:delta)) # for proposal of N
    
    if (N_prior[1] == 'negbin'){
      for (i in seq(1:N_updates)){
      N.cand <- N + sample(eps, size = 1)
      N_try_count <- N_try_count + 1
      if((0 < N.cand) & (N.cand <= 1600)) { # truncated negbin prior.
        
        # modify S
        if (N.cand >= N ){ # append new S
          if (stochastic == FALSE){
          S.cand <- rbind(S, cbind(runif(N.cand - N, xlims[1], xlims[2]), runif(N.cand - N, ylims[1], ylims[2])))
          D.cand <- e2dist1(S.cand, X)
          } else {
            sample_idx <- sample(1:nrow(S))
            S.cand <- S[sample_idx, , drop = FALSE]
            S.cand <- rbind(S.cand, cbind(runif(N.cand - N, xlims[1], xlims[2]), runif(N.cand - N, ylims[1], ylims[2])))
            D.cand <- e2dist1(S.cand, X)
          }
        } 
        else { # N.cand < N
          if (stochastic == FALSE){
          # delete existing S elements
          S.cand <- S[1:N.cand, , drop = FALSE]
          D.cand <- e2dist1(S.cand, X)
          }
          else {
            sample_idx <- sample(1:nrow(S))
            S.cand <- S[sample_idx,]
            S.cand <- S.cand[1:N.cand, ,drop = FALSE]
            D.cand <- e2dist1(S.cand, X)
          }
        }
        
        lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma)) 
        lam.cand <- (lam.cand * (D.cand < thres))
        lamv.cand <- colSums(lam.cand)
        llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
        
        prior <- dnbinom(N, as.numeric(N_prior[2]), as.numeric(N_prior[3]), log=TRUE)
        prior.cand <- dnbinom(N.cand, as.numeric(N_prior[2]), as.numeric(N_prior[3]), log=TRUE)
        
        #if accepted
        if(runif(1) < (exp((llcand + prior.cand) - (ll + prior)))){
          ll <- llcand
          lamv.curr<-lamv.cand
          lam<-lam.cand
          S <- S.cand
          D <- D.cand
          N <- N.cand
          N_acc_count <- N_acc_count + 1
        
        }
      }
    }
    } 
    
    if(N_prior[1] == 'uniform'){
      for (i in seq(1:N_updates)){
        N.cand <- N + sample(eps, size = 1)
        N_try_count <- N_try_count + 1
        if((as.numeric(N_prior[2]) < N.cand) & (N.cand <= as.numeric(N_prior[3]))) {
          
          # modify S
          if (N.cand >= N ){ # append new S
            if (stochastic == FALSE){
              S.cand <- rbind(S, cbind(runif(N.cand - N, xlims[1], xlims[2]), runif(N.cand - N, ylims[1], ylims[2])))
              D.cand <- e2dist1(S.cand, X)
            } else {
              sample_idx <- sample(1:nrow(S))
              S.cand <- S[sample_idx, , drop = FALSE]
              S.cand <- rbind(S.cand, cbind(runif(N.cand - N, xlims[1], xlims[2]), runif(N.cand - N, ylims[1], ylims[2])))
              D.cand <- e2dist1(S.cand, X)
            }
          } 
          else { # N.cand < N
            if (stochastic == FALSE){
              # delete existing S elements
              S.cand <- S[1:N.cand, , drop = FALSE]
              D.cand <- e2dist1(S.cand, X)
            }
            else {
              sample_idx <- sample(1:nrow(S))
              S.cand <- S[sample_idx,]
              S.cand <- S.cand[1:N.cand, ,drop = FALSE]
              D.cand <- e2dist1(S.cand, X)
            }
          }
          
          lam.cand <- lam0*exp(-(D.cand*D.cand)/(2*sigma*sigma)) 
          lam.cand <- (lam.cand * (D.cand < thres))
          lamv.cand <- colSums(lam.cand)
          llcand <- sum(dpois(n, lamv.cand, log=TRUE))
          
          #if accepted
          if(runif(1) < (exp(llcand  - ll))){
            ll <- llcand
            lamv.curr<-lamv.cand
            lam<-lam.cand
            S <- S.cand
            D <- D.cand
            N <- N.cand
            N_acc_count <- N_acc_count + 1
            
          }
        }
      }
      
    }
    
    
    
    out[iter,] <- c(sigma,lam0,N)
    #if(monitorS)
    #  Sout[1:sum(w),,iter] <- S[w==1,]
  }
  last <- list(N=N,lam0=lam0, sigma=sigma, S = S)
  cat("end of the sampling,", "\n")
  end <- proc.time()
  
  sink()
  list(out=out, last=last, Sout=Sout, time = end-start)
}
