rm(list=ls())
setwd("/Users/")
##############################################################################
# 0) Run the functions inside "function_np.R"
##############################################################################

##############################################################################
# 1) library
##############################################################################
library(quadprog)
library(Matrix)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(Rlab)

##############################################################################
# 2) Loading and Standardizing Camera Coordinates (camera -> cam_std)
##############################################################################
camera <- read.csv("TDF_ISRB.csv")
raw_data = read.csv("G_ISRB_v4_day88.csv")

T_all <- ncol(raw_data) - 2
n <- nrow(raw_data)

plot(camera[,2], camera[,3], main = "original UTM", 
     xlab = "UTM_E", ylab = "UTM_N", pch = 19, col = "blue")

min(camera[,2])
min(camera[,3])

# Standardizing
cam_std <- data.frame(
  UTM_E_std = (camera[,2] - 451400) / 500,
  UTM_N_std = (camera[,3] - 4224000) / 500
)

plot(cam_std$UTM_E_std, cam_std$UTM_N_std, main = "",
     xlab = "UTM_E (std)", ylab = "UTM_N (std)", pch = 19, col = "darkgreen")


##############################################################################
# 3) grid setting
##############################################################################
grid_x <- c()
for (i in 1:27) {
  grid_x <- append(grid_x, rep(seq(-0.25, 6.25, length = 27)[i], 19))
}
grid_y <- rep(seq(-0.25, 4.25, length = 19), 27)

grid_center <- data.frame(x = grid_x, y = grid_y)
grid_dim <- nrow(grid_center)  # 35 * 31 = 1085

# visualizing grid and camera (grid_center: black point, cam_extension: red circle)
ggplot() +
  geom_point(data = grid_center,
             aes(x = x, y = y),
             shape = 16, size = 0.5, color = "black") +
  geom_point(data = cam_std,
             aes(x = UTM_E_std, y = UTM_N_std),
             shape = 4, size = 1.5, color = "red") +
  labs(title = "Grid Centers + Cameras",
       x = "X / UTM_E_std", y = "Y / UTM_N_std") +
  theme_minimal()


cam_xy = cam_std
grid_xy   <- grid_center

cam_xy      <- as.matrix(cam_xy)
grid_center <- as.matrix(grid_center)
J           <- nrow(cam_xy)
D2       <- as.matrix(dist(rbind(cam_xy, grid_center)))
dist_mat <- D2[1:J, (J+1):(J + nrow(grid_center))]

dist_round <- round(dist_mat, 2)
d_unique   <- sort(unique(as.vector(dist_round)))

delta_val <- 6 * sqrt(2)

data <- rbind(as.matrix(raw_data[,c(-1,-2)]))
y_data <- as.matrix(data)

lambda0_hat  =quantile((rowSums(y_data>0))/T_all,0.9)

##############################################################################
# 4) model
##############################################################################

res <- alternating_estimation(
    Y        = y_data,
    dist_mat = dist_mat,
    d_unique = d_unique,
    delta    = delta_val,
    lambda0  = lambda0_hat,
    T        = T_all,
    max_iter = 50
  )

res_mean=sum(res$s)

##############################################################################
# 5) bootstrap
##############################################################################
  
B <- 100
Y_base <- y_data
T_all <- ncol(Y_base)

N_boot <- numeric(B)
res_list = vector("list", B)
Y_boot_list = vector("list", B)
lambda0_list = numeric(B)

set.seed(123)
for (b in 1:B) {
    cat("Bootstrap sample", b, "\n")
    
    block_size <- 2
    n_blocks <- ceiling(T_all / block_size)
    max_start <- T_all - block_size + 1
    starts <- sample(seq_len(max_start), size = n_blocks, replace = TRUE)
    b_sample <- unlist(lapply(starts, function(s) seq(s, length.out = block_size)))
    b_sample <- b_sample[1:T_all]
    Y_boot <- Y_base[ , b_sample]

    Y_boot_list[[b]] <- Y_boot
    lambda0_hat  = quantile((rowSums(Y_boot>0))/T_all,0.9)
    lambda0_list[b] = lambda0_hat

    res <- alternating_estimation(
      Y        = Y_boot,
      dist_mat = dist_mat,
      d_unique = d_unique,
      delta    = 8 * sqrt(2),
      lambda0  = lambda0_hat,
      T        = T_all,
      max_iter = 50
    )

    N_boot[b] <- sum(res$s)
    res_list[[b]] <- res
}

cat("\n===== 부트스트랩 결과 요약 =====\n")
cat("N mean:", mean(N_boot), "\n")
cat("N sd:", sd(N_boot), "\n")

y_boot <- sapply(Y_boot_list, sum)

#calculating bootstrap CI
bc_ci <- function(t0, t_boot, alpha=0.05){
  B <- length(t_boot)
  p  <- mean(t_boot < t0) + 0.5*mean(t_boot == t0)
  z0 <- qnorm(p)
  zL <- qnorm(alpha/2); zU <- qnorm(1 - alpha/2)
  aL <- pnorm(2*z0 + zL); aU <- pnorm(2*z0 + zU)
  aL <- pmin(pmax(aL, 1/B), 1 - 1/B)
  aU <- pmin(pmax(aU, 1/B), 1 - 1/B)
  quantile(t_boot, probs = c(aL, aU), type = 6)
}

bc_ci=bc_ci(t0 =res_mean, t_boot = N_boot)

##############################################################################
# 6) Visualization
##############################################################################

# (1) visualizing m
best_s_hat=res$s
s_df <- data.frame(
  x = grid_center[,1],
  y = grid_center[,2],
  s = res$s
)

ggplot(s_df, aes(x = x, y = y)) +
  geom_point(aes(size = ifelse(round(s, 0) > 0, s, NA)), 
             shape = 16, color = "black", alpha = 0.7, na.rm = TRUE) +
  scale_size_continuous(range = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 4, by = 1)) +
  scale_x_continuous(breaks = seq(0, 6, by = 1)) +
  labs(
    x = "UTM_N (std)", y = "UTM_E (std)",
    size = expression(bold(m))
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )


# (2) visualizing p

lambda0 <- quantile((rowSums(y_data>0))/T_all,0.9)
df_lambda_p <- data.frame(
  d = d_unique,
  plambda = res$p * lambda0
)

ggplot(df_lambda_p, aes(x = d, y = plambda)) +
  geom_line(color = "steelblue", size = 1.2) +
  labs(
    x = "distance",
    y = expression(hat(lambda)[0] %.%hat(p)(d)),
    title = ""
  ) +
  scale_x_continuous(limits = c(0, 3)) +
  theme_minimal(base_size = 14)


# (3) visualization of Bootstrap CI
df_N <- data.frame(N_hat = N_boot)

ggplot(df_N, aes(x = N_hat)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", boundary = 0) +
  geom_vline(aes(xintercept = res_mean), color = "blue", linetype = "solid", size = 1.2) +
  geom_vline(xintercept = bc_ci,
             color = "red", linetype = "dotted", size = 1) +
  
  labs(
    title = "",
    x = expression(hat(N)),
    y = "frequency"
  ) +
  theme_minimal(base_size = 20)


# (4) visualization of Bootstrap (p)
lambda0_list <- sapply(Y_boot_list, function(Y_boot) {
  quantile((rowSums(Y_boot > 0)) / T_all, 0.9)
})

p_hat_boot_mat <- sapply(res_list, function(res) res$p)
dim(p_hat_boot_mat)

p_scaled_mat <- sweep(p_hat_boot_mat, MARGIN = 2, STATS = lambda0_list, FUN = "*")
dim(p_scaled_mat)

d_vals <- d_unique

p_mean <- rowMeans(p_scaled_mat)
p_lower <- apply(p_scaled_mat, 1, function(x) mean(x) - 1.96 * sd(x))
p_upper <- apply(p_scaled_mat, 1, function(x) mean(x) + 1.96 * sd(x))

df_p <- data.frame(
  d = d_vals,
  p_mean = p_mean,
  lower = pmax(p_lower,0),
  upper = pmin(p_upper,1)
)

ggplot(df_p, aes(x = d, y = p_mean)) +
  geom_line(color = "steelblue", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "skyblue") +
  labs(
    x = "distance",
    y = expression(hat(lambda)[0]%.%hat(p)(d)),
    title = ""
  ) +
  scale_x_continuous(limits = c(0, 8)) + 
  theme_minimal(base_size = 20)


###########################################################end.