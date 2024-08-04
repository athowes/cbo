# Fig 2iii

set.seed(1)
setwd("C:/Users/Adam/Desktop/bgo/code")

library(MASS)
library(tidyverse)
library(reshape2)
library(lhs)
library(cowplot)
library(tikzDevice)
library(GPfit)

x.star <- seq(0, 1, by = 0.001) # Sequence of points we would like to predict

f <- function(x) {
  - sin(2*x) * sin(5*x) * exp(x-0.5)
} # Define objective function

f.star <- f(x.star) # True values

# Radial basis function with adjustable parameters
rbf <- function(x1, x2, par) {
  # Define parameters
  noise <- par[1]
  length <- par[2]
  
  return(noise^2 * exp(-0.5 * ((x1 - x2)/length)^2))
} 

seed.n <- 6
seed <- optimumLHS(n = seed.n, k = 1)

obs <- data.frame(x = seed,
                  y = f(seed) + rnorm(1, mean = 0, sd = 0.2)
)

x <- obs$x
y <- obs$y

GP <- GP_fit(x,y)

posterior_plot <- function(obs, param) {
  
  x <- obs$x
  y <- obs$y
  
  # Calculate the covariance matrices

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
    # Mean
    geom_line(data = NULL, aes(x = x.star, y = f.star.bar)) +
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

a <- posterior_plot(obs, param = c(0.2, 0.8))
b <- posterior_plot(obs, param = c(1, 0.1))
c <- posterior_plot(obs, param = c(sqrt(0.6290116), sqrt(1/2*0.67245)))

round(sqrt(0.6290116), 3)
round(sqrt(1/2*0.67245), 3)

a <- a + labs(subtitle = "(a) $(\\sigma_f, l) = (0.2, 0.8)$")
b <- b + labs(subtitle = "(b) $(\\sigma_f, l) = (1, 0.1)$")
c <- c + labs(subtitle = "(c) $(\\sigma_f, l) = (0.793, 0.58)$")

plot_grid(a, b, c, ncol = 1, rel_heights = c(10, 10, 11))

tikz("fig2iii.tex", width = 5.5, height = 6)

dev.off()
