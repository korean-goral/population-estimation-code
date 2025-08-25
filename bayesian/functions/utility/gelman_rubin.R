# For gelman-rubin test
burn_in <- 10000

# target indexes
ind <- c(1:3)
chain_list <- list()
for (i in ind){
  chain <- list(mcmc(results[[i]]$out[burn_in:n_iter,]))
  chain_list <- append(chain_list, chain)
  }
chains <- mcmc.list(chain_list)
gelman_result <- gelman.diag(chains)
print(gelman_result)

#gelman_plots
savegraph <- paste0('utility/gelman_plot', format(Sys.time(), "%d-%H-%M-%S"),'.png')
CairoPNG(savegraph, width=1200, height=300 * length(results))
gelman.plot(chains)
dev.off()
