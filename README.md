### 야생동물 개체군 밀도 추정 통계 모델 코드 
bayesian 폴더는 Reversible Jump MCMC 방법을 이용한 베이지안 모형의 코드이고 
<br> poisson 폴더는 거리별 탐지 확률을 비모수적으로 추정하는 포아송 회귀 모형의 코드입니다.

## Bayesian method Introduction - RJMCMC 

작성자 : 정우진

아래는 RJMCMC 코드 활용법에 대한 간단한 설명과 예시입니다.

### 파일 구조

- data : 살제 산양데이터와 전처리 코드
- functions : 알고리즘 구현 파일들
  - sampling functions/rjmcmc.R : RJMCMC 알고리즘
  - utility : 기타 기능코드 들(gelman-rubin test, 시각화 등)
- main.R : 산양데이터 분석 코드
- result : 결과 출력이 저장된 파일
- simulation : 가상실험 관련 코드

산양 데이터 분석을 재현 할 때는 main.R을 통하여 실행하면 됩니다.

### main.R 사용법

``` r
library(parallel) 
source('functions/sampling_functions/rjmcmc.R')

# iteration과 burnin 설정

n_iter <- 1000 
burn_in <- 100
rep <- 3

task_list <- list( list(fun = function() rjmcmc(n = goral$n, X = goral$X, N_prior = c('negbin', 1,0.002), sigma_prior = c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim, burn_in = burn_in, log = 'realdata.log')), 

list(fun = function() rjmcmc(n = goral$n, X = goral$X, N_prior = c('negbin', 1,0.002), sigma_prior = c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim, burn_in = burn_in, log = 'realdata.log')),

list(fun = function() rjmcmc(n = goral$n, X = goral$X, N_prior = c('uniform', 1, 300), sigma_prior = c('gamma', 6, 1/27),
lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim, burn_in = burn_in, log = 'realdata.log')) )

task_list <- unlist(lapply(task_list, function(x) rep(list(x),rep)), recursive = FALSE)

cl <- makeCluster(20)
clusterExport(cl, varlist = c('goral', 'rjmcmc', 'e2dist1','n_iter','burn_in'))

# CHECK RESULTS

savefile <- paste0('results/data/real/',format(Sys.time(), "%d-%H-%M.RData"))

results <- parLapply(cl, task_list, function(task) { 
  tryCatch({
  task$fun() 
    }, error = function(e) {
      message("Error in task: ",
      conditionMessage(e))
      return(NA)}) 
      }
      )

save(results, file = savefile) 
stopCluster(cl) 
```

위 코드는 RJMCMC 함수를 다수 cpu를 이용하여 여러개의 chain을 여러개
돌리고자 병렬로 구현한 것이다.

만약 병렬처리가 불가한 상황이라면,

``` r
result <- rjmcmc(n = goral$n, X = goral$X, N_prior = c('negbin', 1,0.002), sigma_prior = c('gamma', 6, 1/27), lam0_prior = c(1.5, 1.5), niters = n_iter, xlims = goral$xlim, ylims = goral$ylim, burn_in = burn_in, log = 'realdata.log')
```

위와같이 단순 함수만 사용하여 실행하면 된다

위의 병렬코드에서,

- n_iter : 돌리고자 하는 총 iteration 수

- burn_in : 사후 번인 표본 수

- rep : 동일한 체인을 몇 개 돌릴지 결정 \<- gelman-rubin 통계량을
  이용하여 수렴여부를 판단하기 위해 필요하다.

``` r
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
```

main.r의 다음부분은 chain 결과를 해석하기 위한 것이다.

- draw_traceplot : chain의 traceplot을 그려주며, 결과물은
  results/plots/out.png로 저장된다. 여기서 out 파일명은 사용자가
  지정한다.

- draw_traceplot : chain의 사후표본 밀도함수를 그려주며, 결과물은
  results/plots/out.png로 저장된다. 역시 out 파일명은 사용자가 지정한다.

- Gelman-rubin 통계량을 통한 수렴여부 확인을 위해서, ind를 동일한
  셋팅에서 나온 체인들로 설정해준다. 만약 result의 1,2,3 번째 체인이
  동일한 설정에서 나왔다면, ind \<- c(1,2,3)으로 두고 나머지 코드를
  실행하면 된다.

  - 즉, 총 9개 체인이 있는데 , 123/456/789가 각각 동일한 설정이라면, ind
    \<- c(1,2,3) / c(4,5,6)/ c(7,8,9)로 총 세번을 수행해서 결과값 세개를
    얻으면 된다.

- summarize 함수는 결과값 요약통계량을 보여주며, 동일한 환경에서 나온
  체인들의 통계량은 평균을 내서 출력한다

!! 알수 없는 이유로 몇 개의 chain이 NA 를 반환하는 경우가 종종 있다. 이
경우, 요약통계량 산출 시 functions/utility/summarize 함수에서 ind를 위의
gelman-rubin 통계량을 구할 때와 같이 동일한 환경에서 나온 체인으로
설정하여 여러번 결과값을 출력하도록 한다.

### rjmcmc 함수 설명

- rjmcmc 함수는 다음과 같은 인자들을 가지고 있다.

``` r
  rjmcmc(n, X, N_prior, sigma_prior, niters, xlims, ylims, thres = 1, burn_in = 10000, init_delta = 10, init_tune=c(0.1, 0.1, 2), N_updates = 15, log, stochastic = FALSE, lam0_prior = FALSE, check_point = FALSE, monitorS=FALSE)
```

- n : 카메라 트랩 데이터

- X : 카메라 트랩 위치

- N_prior : N에 대한 사전분포

  - ex c(‘negbin’, 1, 0.001) : N 사전분포를 negbin(1, 0.001)로 설정

  - c(’uniform, 1, 500) : N 사전분포를 uniform(1, 500)으로 설정

  - N 사전분포는 uniform과 음이항분포 두 종류 사용 가능하다.

- sigma_prior : sigma에 대한 사전분포

  - ex c(‘gamma’, 6, 1/27) : N 사전분포를 gamma(6, 1/27)로 설정

  - 이 때 1/27은 scale 파라미터로, 위 분포의 평균은 6/27임에 주의

  - c(’uniform, 1, 500) : sigma 사전분포를 uniform(1, 500)으로 설정

  - sigma 사전분포는 uniform과 감마분포 두 종류 사용 가능하다.

- niters: iteration 수

- xlims : Activity center S 샘플링 x범위

<!-- -->

- ylims : Activity center S 샘플링 y범위

- thres : 동물이 camera에 포착될 수 있는 activity center와 camera 거리
  한계

  - thres = 1 : 1km.

- burn_in ; 번인 수

- \<기타 tuning parameters\>

  - init_delta, init_tune : proposal distribution 초기 파라미터를
    지정함.

  - N_updates : 한번의 iteration에서 N을 몇번 업데이트 할지 지정

  - log : 로그파일 위치 지정 : results/logs/log.log에 에서 확인 가능

