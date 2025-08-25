rm(list = ls())


source('data/goral_data.R') # 분석대상 데이터 불러오기
source('functions/sampling_functions/e2dist1.R')

################################### Functions ################################### ################################### ################################### 

# 시뮬레이션 데이터 생성함수
generate.data <- function(N, sigma, lam0, thres, t){
  xrange <- goral$xlim
  yrange <-goral$ylim
  X <- goral$X
  S <- cbind(runif(N, xrange[1], xrange[2]), runif(N, yrange[1], yrange[2]))
  dist <- e2dist1(X,S) # 각 카메라와 activity center 사이의 거리
  lam <- lam0*exp(-(dist*dist)/(2*sigma*sigma))
  lam <- (lam * (dist < thres))
  lamv <- rowSums(lam)
  n.generated <- t(sapply(lamv, function(l) rpois(t, lambda = l)))
  test.data <- list(n=n.generated, X=X, J = nrow(n), K = ncol(n), xlim = xrange, ylim = yrange)
  return(test.data)
}

# 파라미터 조합
params.list <- list(c(100,0.23, 0.1), c(100, 0.23, 0.2), c(100, 0.23, 0.3), c(100, 0.23, 0.4),
                 c(200,0.23, 0.1), c(200, 0.23, 0.2), c(200, 0.23, 0.3), c(200, 0.23, 0.4),
                 c(300,0.23, 0.1), c(300, 0.23, 0.2), c(300, 0.23, 0.3), c(300, 0.23, 0.4))


################################### 병렬처리 ################################### ################################### ################################### 
# iteration과 burnin 설정
n_iter <- 1000
burn_in <- 100
rep <- 3

library(parallel)
source('functions/sampling_functions/rjmcmc.R')

for (params in params.list[1]){ # for some parameter setting
  N <- params[1]
  sigma <- params[2]
  lam0 <- params[3]
  
  test.data <- generate.data(N, sigma, lam0, thres = 1, t = 88)
  
  task_list <- list(
    list(fun = function() rjmcmc(n =  test.data$n, X = test.data$X, N_prior =  c('negbin', 1,0.002), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = test.data$xlim, ylims = test.data$ylim,
                                                        burn_in = burn_in, log = 'simulation.log')),
    list(fun = function() rjmcmc(n =  test.data$n, X = test.data$X, N_prior =  c('negbin', 1,0.002), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = test.data$xlim, ylims = test.data$ylim,
                                                        burn_in = burn_in, log = 'simulation.log')),
    list(fun = function() rjmcmc(n =  test.data$n, X = test.data$X, N_prior =  c('uniform', 1, 300), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = test.data$xlim, ylims = test.data$ylim,
                                                        burn_in = burn_in, log = 'simulation.log'))
    )

  
task_list <- unlist(lapply(task_list, function(x) rep(list(x),rep)), recursive = FALSE)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('test.data','rjmcmc', 'e2dist1', 'n_iter', 'burn_in'))

# CHECK RESULTS
savefile <- paste0('results/data/simulation/', 'N',  N, 'sigma', sigma, 'lam', lam0, '_', format(Sys.time(), "%d-%H-%M.RData"))

results <- parLapply(cl, task_list, function(task) {
  tryCatch({
    task$fun()
  }, error = function(e) {
    message("Error in task: ", conditionMessage(e))
    return(NA)  # 또는 NULL
  })
})

save(results, file = savefile)
stopCluster(cl)
}

################################### 결과 확인 ################################### ################################### ################################### 

# traceplot / density 확인 

## !!! 모종의 이유로 일부 chain이 정상적으로 출력되지 않은 경우 : gelman-rubin 통계량 계산과 요약통계량 산출에서 관심대상 체인만 가지고 분석해주세요.

source('functions/utility/plotting.R')

draw_traceplot(results, n_iter, burn_in, out = 'traceplot.png')
draw_density(results, n_iter, burn_in, out = 'density.png')


# gelman-rubin 통계량을 통한 수렴진단.
library(coda)
# ex) 1,2,3 번째 chain이 동일한 셋팅 하에서 생성된 체인이라면, 이 세개의 체인을 통해 gelman-rubin 통계량을 계산합니다.

ind <- c(1:3)
chain_list <- list()
for (i in ind){
  chain <- list(mcmc(results[[i]]$out[burn_in:n_iter,]))
  chain_list <- append(chain_list, chain)
}

chains <- mcmc.list(chain_list)
gelman_result <- gelman.diag(chains)
print(gelman_result)

# 요약 통계량 산출 

source('functions/utility/summarize.R')
summarize(results, n_iter, burn_in, rep)
