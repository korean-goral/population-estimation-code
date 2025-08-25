#plot

library(Cairo)

draw_traceplot <- function(data, n_iter, burn_in, out){

  savegraph <- paste0('results/plots/',out)
  CairoPNG(savegraph, width=1200, height=300 * length(results))
  par(mfrow=c(length(results), 3)) 
  for (i in 1:length(results)){
    if (is.null(results[[i]]) || anyNA(results[[i]])) {
      for (j in 1:3) {
        plot.new()
        text(0.5, 0.5, paste("No data for result", i), cex=1.5)
      }
    } else {
    mc_samples <- results[[i]]$out
    plot(mc_samples[burn_in:n_iter,1], type="l", ylab="sigma")
    plot(mc_samples[burn_in:n_iter,2], type="l", ylab="lam0")
    plot(mc_samples[burn_in:n_iter,3], type="l", ylab="N")
  }}
  dev.off()
}

#density plot

draw_density <- function(data, n_iter, burn_in, out){
  savegraph <- paste0('results/plots/', out)
  
  CairoPNG(savegraph, width=1200, height=300 * length(results))
  par(mfrow=c(length(results), 3))
  
  for (i in 1:length(results)){
    if (is.null(results[[i]]) || anyNA(results[[i]])) {
      for (j in 1:3) {
        plot.new()
        text(0.5, 0.5, paste("No data for result", i), cex=1.5)
      }
    } else {
      plot(density(results[[i]]$out[burn_in:n_iter,1]), xlab="sigma", main="")
      plot(density(results[[i]]$out[burn_in:n_iter,2]), xlab="lam0", main="")
      hist(results[[i]]$out[burn_in:n_iter,3], xlab="N", freq=FALSE, main="")
    }
  }
  dev.off()
}