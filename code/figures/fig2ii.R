# Fig 2ii

set.seed(1)
setwd("C:/Users/Adam/Desktop/bgo/code")

library(MASS)
library(tidyverse)
library(reshape2)
library(lhs)
library(cowplot)
library(tikzDevice)

x.star <- seq(0, 1, by = 0.001) # Sequence of points we would like to predict

f <- function(x) {
  sin(5*x) * sin(13*x) * exp(x-0.5)
} # Define objective function

f.star <- f(x.star) # True values

# Radial basis function with adjustable parameters
rbf <- function(x1, x2, par) {
  # Define parameters
  noise <- par[1]
  length <- par[2]
  
  return(noise^2 * exp(-0.5 * ((x1 - x2)/length)^2))
} 

# Prior plot

# Covariance matrix
k.sxsx <- outer(x.star, x.star, rbf, par = c(1, 0.15))

f.star.bar <- 0 # Set the predictive mean to zero

n <- 3
values <- matrix(rep(0, n*length(x.star)), ncol = n)
for (i in 1:n) {
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), k.sxsx)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id="x")

prior_plot <- ggplot() +
  geom_line(data = values, aes(x = x, y = value, group = variable, col = variable)) + 
  # The true function
  geom_line(data = NULL, aes(x = x.star, y = f.star), col = "grey40", lty = 2) +
  # 95% Confidence interval
  geom_ribbon(data = NULL, aes(x = x.star,  
                               ymin = f.star.bar - 1.96*sqrt(diag(k.sxsx)),
                               ymax = f.star.bar + 1.96*sqrt(diag(k.sxsx))),
              fill = "grey80", alpha = 0.3) +
  scale_y_continuous(lim = c(-2.2, 2.2)) +
  labs(y = "$y$") +
  guides(col = FALSE) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Posterior plot

sigma.n <- 0.1
seed.n <- 5
seed <- optimumLHS(n = seed.n, k = 1)
obs <- data.frame(x = seed,
                  y = f(seed) + rnorm(1, sd = sigma.n)
)

seed2 <- optimumLHS(n = seed.n, k = 1)
obs2 <- data.frame(x = seed2,
                  y = f(seed2) + rnorm(1, sd = sigma.n)
)

posterior_plot <- function(obs) {
  
  x <- obs$x
  y <- obs$y
  # Calculate the covariance matrices
  param <- c(1, 0.15)

  k.xx <- outer(x, x, rbf, par = param)
  k.xxs <- outer(x, x.star, rbf, par = param)
  k.xsx <- outer(x.star, x, rbf, par = param)
  k.xsxs <- outer(x.star, x.star, rbf, par = param)

  # Conditioning on the data
  f.star.bar <- k.xsx %*% solve(k.xx + sigma.n^2*diag(1, ncol(k.xx))) %*% y
  cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx + sigma.n^2*diag(1, ncol(k.xx))) %*% k.xxs

  post.n <- 3
  values <- matrix(rep(0, post.n*length(x.star)), ncol = post.n)
  for (i in 1:post.n) {
    values[,i] <- mvrnorm(1, f.star.bar, cov.f.star)
  }
  values <- cbind(x = x.star, as.data.frame(values))
  values <- melt(values, id="x")

  plot_out <- ggplot() +
  # Samples from the posterior
  geom_line(data = values, aes(x = x, y = value, group = variable, col = variable)) + 
  # The true function
  geom_line(data = NULL, aes(x = x.star, y = f.star), col = "grey40", lty = 2) +
  # 95% Confidence interval
  geom_ribbon(data = NULL, aes(x = x.star,  
                               ymin = f.star.bar - 2*sqrt(diag(cov.f.star)),
                               ymax = f.star.bar + 2*sqrt(diag(cov.f.star))),
              fill = "grey80", alpha = 0.3) +
  # Observations
  geom_point(data = NULL, aes(x = obs$x, y = obs$y)) +
  scale_y_continuous(lim = c(-2.2, 2.2), name = "$y$") +
  labs(x = "$x$") +
  guides(col = FALSE) +
  theme_minimal()
  
  return(plot_out)
}

five <- posterior_plot(obs)
ten <- posterior_plot(rbind(obs, obs2))

tikz("fig2iia.tex", width = 5.5, height = 1.9)
prior_plot
dev.off()

tikz("fig2iib.tex", width = 5.5, height = 1.9)
five
dev.off()

tikz("fig2iic.tex", width = 5.5, height = 2.2)
ten
dev.off()
