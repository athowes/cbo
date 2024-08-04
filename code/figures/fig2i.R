# Fig 2i

set.seed(1)
setwd("C:/Users/Adam/Desktop/bgo/code")

library(MASS)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(RandomFieldsUtils)
library(tikzDevice)

rbf <- function(x1, x2, par) {
  # Define parameters
  noise <- par[1]
  length <- par[2]
  
  return(noise^2 * exp(-0.5 * ((x1 - x2)/length)^2))
} # Radial basis function with adjustable parameters

mat <- function(x1, x2, par) {
  # Define parameters
  noise <- par[1]
  length <- par[2]
  nu <- par[3]
  
  x <- (x1 - x2) / length^2
  return(noise^2 * matern(x, nu))
} # Matern kernel with adjustable parameters

gp.demo <- function(kernel, par) {
  
  # Sequence of points we would like to predict
  x.star <- seq(0, 10, by = 0.01)
  
  # Covariance matrix
  k.sxsx <- outer(x.star, x.star, kernel, par = par)
  
  # Take samples
  n <- 3
  values <- matrix(rep(0, n*length(x.star)), ncol = n)
  for (i in 1:n) {
    values[,i] <- mvrnorm(1, rep(0, length(x.star)), k.sxsx)
  }
  values <- cbind(x = x.star, as.data.frame(values))
  values <- melt(values, id="x")

  plot <- ggplot() +
    geom_line(data = values, aes(x = x, y = value, group = variable, col = variable)) + 
    scale_y_continuous(lim = c(-4, 4), name=NULL) +
    guides(col = FALSE) +
    theme_minimal() +
    theme(axis.title.x=element_blank())
  
  return(plot)
}

length5 <- gp.demo(kernel = rbf, par = c(1, 5)) + labs(subtitle = "(a) SE with $l = 5$")
length1 <- gp.demo(kernel = rbf, par = c(1, 1)) + labs(subtitle = "(c) SE with $l = 1$")
length0.2 <- gp.demo(kernel = rbf, par = c(1, 0.2)) + labs(subtitle = "(e) SE with $l = 0.2$")

nu5over2 <- gp.demo(kernel = mat, par = c(1, 1, 5/2)) + labs(subtitle = "(b) Matérn with $\\nu = 5/2$")
nu3over2 <- gp.demo(kernel = mat, par = c(1, 1, 3/2)) + labs(subtitle = "(d) Matérn with $\\nu = 3/2$")
nu1over2 <- gp.demo(kernel = mat, par = c(1, 1, 1/2)) + labs(subtitle = "(f) Matérn with $\\nu = 1/2$")

tikz("fig2i.tex", width = 6, height = 6)
grid.arrange(length5, nu5over2, length1,  nu3over2, length0.2, nu1over2, ncol = 2, nrow = 3)
dev.off()
