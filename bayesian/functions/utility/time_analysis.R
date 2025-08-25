
ind <- c(10:18)
times <- c()
for (i in ind) {
  times <- c(times, results[[i]]$time[3])
}

mean(times)
