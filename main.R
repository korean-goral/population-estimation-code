rm(list = ls())

# 분석대상 데이터
source('data/goral_data.R') # 분석대상 데이터 불러오기
source('functions/sampling_functions/e2dist1.R')

################################### 병렬처리 ################################### ################################### ################################### 
library(parallel)
source('functions/sampling_functions/rjmcmc.R')

# iteration과 burnin 설정
n_iter <- 1000
burn_in <- 100
rep <- 3

task_list <- list(
  list(fun = function() rjmcmc(n =  goral$n, X = goral$X, N_prior =  c('negbin', 1,0.002), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim,
                               burn_in = burn_in, log = 'realdata.log')),
  list(fun = function() rjmcmc(n =  goral$n, X = goral$X, N_prior =  c('negbin', 1,0.002), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim,
                               burn_in = burn_in, log = 'realdata.log')),
  list(fun = function() rjmcmc(n =  goral$n, X = goral$X, N_prior =  c('uniform', 1, 300), sigma_prior =  c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim,
                               burn_in = burn_in, log = 'realdata.log'))
)

task_list <- unlist(lapply(task_list, function(x) rep(list(x),rep)), recursive = FALSE)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('goral', 'rjmcmc', 'e2dist1','n_iter','burn_in'))

# CHECK RESULTS
savefile <- paste0('results/data/real/',format(Sys.time(), "%d-%H-%M.RData"))

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
