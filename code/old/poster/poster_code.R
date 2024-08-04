# Poster, Bayesian Global Optimization code
# Adam Howes

# Libraries ---------------------------------------------------------------

library(MASS)
library(tidyverse)
library(reshape2)
library(lhs)

# Main --------------------------------------------------------------------

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

# Objective function plot -------------------------------------------------

# A plot of the objective function
ggplot(data = NULL, aes(x = x.star, y = f.star)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(-2.2,2.2))

# Prior funciton plot --------------------------------------------------------------

# Covariance matrix
k.sxsx <- outer(x.star, x.star, rbf, par = c(1, 0.15))

f.star.bar <- 0 # Set the predictive mean to zero

n <- 4
values <- matrix(rep(0, n*length(x.star)), ncol = n)
for (i in 1:n) {
  values[,i] <- mvrnorm(1, rep(0, length(x.star)), k.sxsx)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id="x")
  
ggplot() +
  geom_line(data = values, aes(x = x, y = value, group = variable, col = variable), size = 2) + 
  scale_color_brewer(palette="Oranges") +
  # The true function
  geom_line(data = NULL, aes(x = x.star, y = f.star), size = 2, col = "grey40", lty = 2) +
  # 95% Confidence interval
  geom_ribbon(data = NULL, aes(x = x.star,  
                               ymin = f.star.bar - 1.96*sqrt(diag(k.sxsx)),
                               ymax = f.star.bar + 1.96*sqrt(diag(k.sxsx))),
              fill = "grey80", alpha = 0.3) +
  scale_y_continuous(lim = c(-2.2, 2.2), name=NULL) +
  xlab("input, x") +
  guides(col = FALSE) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())

# Get one data ------------------------------------------------------------

x <- 0.6141048
y <- f(x)

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

ggplot() +
  # Samples from the posterior
  geom_line(data = values, aes(x = x, y = value, group = variable, col = variable), size = 2) + 
  # The true function
  geom_line(data = NULL, aes(x = x.star, y = f.star), size = 2, col = "grey60", lty = 2) +
  # 95% Confidence interval
  geom_ribbon(data = NULL, aes(x = x.star,  
                               ymin = f.star.bar - 2*sqrt(diag(cov.f.star)),
                               ymax = f.star.bar + 2*sqrt(diag(cov.f.star))),
              fill = "grey80", alpha = 0.3) +
  # Observations
  geom_point(data = NULL, aes(x = x, y = y), size = 4) +
  scale_y_continuous(lim = c(-2.2, 2.2), name=NULL) +
  guides(col = FALSE) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Add data ----------------------------------------------------------------

sigma.n <- 0.05
seed.n <- 3
seed <- optimumLHS(n = seed.n, k = 1)
obs <- data.frame(x = seed,
                  y = f(seed) + rnorm(1, sd = sigma.n)
)
x <- obs$x
y <- obs$y

y

# Calculate the covariance matrices
param <- c(1, 0.15)

k.xx <- outer(x, x, rbf, par = param)
k.xxs <- outer(x, x.star, rbf, par = param)
k.xsx <- outer(x.star, x, rbf, par = param)
k.xsxs <- outer(x.star, x.star, rbf, par = param)

# Conditioning on the data
f.star.bar <- k.xsx %*% solve(k.xx + sigma.n^2*diag(1, ncol(k.xx))) %*% y
cov.f.star <- k.xsxs - k.xsx %*% solve(k.xx + sigma.n^2*diag(1, ncol(k.xx))) %*% k.xxs

post.n <- 4
values <- matrix(rep(0, post.n*length(x.star)), ncol = post.n)
for (i in 1:post.n) {
  values[,i] <- mvrnorm(1, f.star.bar, cov.f.star)
}
values <- cbind(x = x.star, as.data.frame(values))
values <- melt(values, id="x")

ggplot() +
  # Samples from the posterior
  geom_line(data = values, aes(x = x, y = value, group = variable, col = variable), size = 2) + 
  scale_color_brewer(palette="Oranges") +
  # The true function
  geom_line(data = NULL, aes(x = x.star, y = f.star), size = 2, col = "grey40", lty = 2) +
  # 95% Confidence interval
  geom_ribbon(data = NULL, aes(x = x.star,  
                               ymin = f.star.bar - 2*sqrt(diag(cov.f.star)),
                               ymax = f.star.bar + 2*sqrt(diag(cov.f.star))),
              fill = "grey80", alpha = 0.3) +
  # Observations
  geom_point(data = NULL, aes(x = obs$x, y = obs$y), size = 4) +
  # Acquisition
  geom_ribbon(data = NULL, aes(x = x.star, ymin = -2, ymax = ci2), 
              size = 2, alpha = 0.3, fill = "#0070CC") +
  # Acquistion max
  geom_point(data = NULL, aes(x = 0.323, y = -2), shape = 17, size = 6, col = 'red') +
  scale_y_continuous(lim = c(-2.2, 2.2), name=NULL) +
  guides(col = FALSE) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())


# Acquisition function ----------------------------------------------------

y.max <- max(y)
sd.f.star <- sqrt(diag(cov.f.star))
std <- (f.star.bar - y.max)/sd.f.star

# Probability of improvment
pi <- pnorm(std)
pi2 <- pi/max(pi)

# Expected improvement
ei <- (f.star.bar - y.max) * pnorm(std) + sd.f.star * dnorm(std)
ei2 <- ei/max(ei)

# Confidence interval
ci <- f.star.bar + 1.95* sd.f.star
ci2 <- ci/max(ci) - 2
which(ci2 < -2)
ci2[687:715,] <- -2

ggplot() +
  geom_line(data = NULL, aes(x = x.star, y = ei2), linetype = "dotted", size = 2) +
  geom_line(data = NULL, aes(x = x.star, y = ci2), size = 2) +
  geom_line(data = NULL, aes(x = x.star, y = pi2), col = "green", size = 2) +
  scale_y_continuous(limits = c(-1,1)) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
