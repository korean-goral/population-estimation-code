##############################################################################
# File defining the basic required functions
##############################################################################

##############################################################################
# update of m (number of living wildlife in the grid)
##############################################################################

update_s <- function(X, Y, s_old, lambda = 0, max_iter = 30) {
  if (any(!is.finite(X))) stop("X contains non-finite values!")
  p <- ncol(X)
  idx_y_zero <- which(Y == 0)
  
  mu_old <- as.numeric(X %*% s_old)
  mu_old <- pmax(mu_old, 1e-6)
  ll_old <- sum(Y * log(mu_old) - mu_old)
  
  for (iter in 1:max_iter) {
    s_body <- s_old
    X_body <- X
    
    mu <- as.numeric(X_body %*% s_body)
    mu[mu < 1e-6] <- 1e-6
    mu[idx_y_zero]=1
    
    W <- Diagonal(x = Y/(mu^2))
    grad <- t(X_body) %*% (Y/mu-1)
    H <- t(X_body) %*% W %*% X_body- lambda * rep(1, p)
    
    H <- (H + t(H)) / 2
    tryCatch({ chol(H) }, error = function(cond) {
      eig <- eigen(H, symmetric = TRUE)
      min_eig <- min(eig$values)
      H <<- H + (-min_eig + 0.001) * diag(ncol(H))
    })
    
    Dmat <- H
    Dmat_pd <- as.matrix( nearPD(as.matrix(H), corr = FALSE)$mat )
    
    dvec <- as.numeric(Dmat_pd %*% s_body + grad)
    
    Amat <- diag(length(s_body))
    bvec <- rep(0, length(s_body))
    
    sol <- tryCatch(
      solve.QP(Dmat_pd, dvec, Amat, bvec),
      error = function(e) NULL
    )
    if (is.null(sol)) {
      warning("Iteration ", iter, ": QP failed, fallback but keep going")
      
      s_body_new <- s_body + runif(length(s_body), -1e-8, 1e-8)
      s_body_new[s_body_new < 0] <- 0

      s_old <- s_body_new
      next
    } else {
      s_body_new <- sol$solution
      s_body_new[s_body_new < 1e-7] <- 0
    }
    
    s_new <- s_body_new
    
    mu_new <- as.numeric(X %*% s_new)
    mu_new <- pmax(mu_new, 1e-6)
    ll_new <- sum(Y * log(mu_new) - mu_new)

    if (abs(ll_new - ll_old) < 1e-3) {
      message("Converged by log-likelihood at iter ", iter,
              " (Δll = ", round(ll_new-ll_old, 6), ")")
      s_old = s_new
      break
    }
    ll_old <- ll_new
    s_old <- s_new
    print(iter)
  }
  
  return(s_old)
}


##############################################################################
# 3) update of p (detection probability as a function of distance)
##############################################################################
make_M <- function(dist_round, s, d_unique, delta,T) {
  K <- length(d_unique)
  J <- nrow(dist_round)
  M <- matrix(0, J, K)
  for (k in seq_len(K)) {
    mask   <- (dist_round == d_unique[k]) & (dist_round <= delta)
    M[, k] <- mask %*% s
  }
  do.call(rbind, replicate(T, M, simplify = FALSE)) 
}

#Monotonic and non-negative constraint matrix for p (length K)
make_constraints <- function(K, tol = 1e-6) {
  A1 <- matrix(0, K-1, K)
  for (i in 1:(K-1)) {
    A1[i, i]   <-  1
    A1[i, i+1] <- -1
  }
  b1  <- rep(-tol, K-1)
  Ain <- rbind(A1, diag(K))
  bin <- c(b1, rep(0, K))
  list(Amat = t(Ain), bvec = bin)
}


update_p_fast <- function(Y, dist_round, d_unique, s_vec, delta, T, cons) {
  #(warning) Y must be a (TxJ)x1 vector.
  M     <- make_M(dist_round, s_vec, d_unique, delta,T) #(J x T) x K
  D0    <- crossprod(M) #M'M
  ridge <- max(1e-6, mean(diag(D0)) * 1e-8)
  Dmat  <- (D0 + t(D0)) / 2 + diag(ridge, ncol(D0))
  dvec  <- as.numeric(crossprod(M, Y))
  sol   <- solve.QP(Dmat, dvec, cons$Amat, cons$bvec, meq = 0)$solution
  p2    <- pmin(pmax(sol, 0), 1)
  p3    <- cummin(p2)
  p3 / p3[1]
}

##############################################################################
# Final Update: Alternating Optimization of m and p
##############################################################################

alternating_estimation <- function(Y, dist_mat, d_unique, delta, lambda0, T = 1, max_iter = 30) {

  Y_vec <- as.vector(Y)
  dist_round <- round(dist_mat, 2)       # n × grid_dim
  K          <- length(d_unique)
  cons       <- make_constraints(K)
  k_mat      <- matrix(
    match(as.vector(dist_round), d_unique),
    nrow = nrow(dist_round), 
    ncol = ncol(dist_round)
  )
  
  # initial values
  grid_dim <- ncol(dist_mat)
  p=seq(1,0,length.out=K)

  s=rep(0.0000001, grid_dim)

  X0_vals    <- p[k_mat]*lambda0
  X0      <- matrix(X0_vals, 
                   nrow = nrow(dist_round), 
                   ncol = grid_dim)
  X0[dist_round > delta] <- 0
  X_full0 <- do.call(rbind, replicate(T, X0, simplify = FALSE))
  
  mu0   <- X_full0 %*% s
  mu0   <- pmax(mu0, 1e-6)
  ll_old <- sum(Y_vec * log(mu0) - mu0)
  
  Nhat_old = sum(s)
  
  history  <- list()
  
  for (r in seq_len(max_iter)) {
    # update : p
    p_old = p
    p <- update_p_fast(Y_vec, dist_round, d_unique, s,delta,T, cons)
    
    X_vals <- p[k_mat]*lambda0                 # length = n * grid_dim
    X      <- matrix(X_vals, 
                     nrow = nrow(dist_round), 
                     ncol = grid_dim)
    X[dist_round > delta] <- 0
    X_full <- do.call(rbind, replicate(T, X, simplify = FALSE))
    
    # update : m
    s_old = s
    s_new <- update_s(X_full, Y_vec, s_old)
    Nhat_new = sum(s_new)
    
    rel_diff <- abs(Nhat_new - Nhat_old) / (Nhat_old + 1e-6)
    if (rel_diff > 100000) {
      warning("Nhat spiked at iter ", r, ": rollback to previous s.")
      s <- s_old
      Nhat <- sum(s_old)
    } else {
      s <- s_new
      Nhat <- Nhat_new
      Nhat_old = Nhat_new
    }
    
    mu_hat <- as.numeric(X %*% s)
    mu_hat[mu_hat < 1e-6] <- 1e-6
    ll_new = sum(Y_vec * log(mu_hat) - mu_hat)

    pearson_gof <- sum((Y - mu_hat)^2 / mu_hat)
    
    cat("=== iter:", r, " N_hat:", Nhat, " GOF:", round(pearson_gof, 4), "\n")
    history[[r]] <- list(s = s, p = p, N = Nhat, GOF = pearson_gof, ll=ll_new)
    
    if (abs(ll_new - ll_old) < 1e-3 && abs(Nhat-Nhat_old)<1e-2) {
      message("Converged by log-likelihood at iter ", r,
              " (Δll = ", round(ll_new - ll_old, 6),
              ", ΔN = ", round(Nhat - Nhat_old, 6), ")")
      
      break
    }
    ll_old <- ll_new
   }
  
  list(s = s, p = p, history = history)
}



