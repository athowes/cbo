# Fig 1ii

set.seed(1)
setwd("C:/Users/Adam/Desktop/bgo/code")

library(cowplot)
library(GPfit)
library(lhs)
library(tidyverse)
library(tikzDevice)

f <- function(x) {
  z <- x*10
  z*sin(z) + z*cos(2*z)
} # Define objective function

gp.data <- function(f, obs) {
  
  x <- obs$x
  y <- obs$y
  y.min <- min(y)
  
  GPmodel <- GP_fit(x,y)
  
  x.star <- seq(0, 1, by = 0.0001) # Sequence of points to predict
  y.pred <- predict.GP(GPmodel, x.star) # GP predictions
  
  # Extract information
  df <- as.data.frame(y.pred$complete_data) 
  names(df) <- c("x.star", "y.pred", "mse")
  df <- mutate(df, 
               y.upper = y.pred + 2*sqrt(mse), 
               y.lower = y.pred - 2*sqrt(mse),
               ei = (y.min - y.pred) * pnorm((y.min - y.pred)/sqrt(mse)) +
                 sqrt(mse) * dnorm((y.min - y.pred)/sqrt(mse)))
  return(df)
}

# 5 points ----------------------------------------------------------------

# Design points
des.x <- c(0.4554499588, 0.0001036626, 0.7021216746, 0.2028095817, 0.8129379553)
obs <- data.frame(x = des.x,
                  y <- f(des.x))
x <- obs$x
y <- obs$y

df <- gp.data(f, obs)

# Find acquisition maximum
acq.max <- df[which(df$ei == max(df$ei)),]

x.star <- seq(0, 1, by = 0.0001) # Sequence of points to predict
f.star <- f(x.star)   # True values

# Plot objective
obj.plot <- ggplot() +
  # True function
  geom_line(data = NULL, aes(x = x.star, y = f.star), col = "grey60", lty = 2) +
  # Prediction mean
  geom_line(data = df, aes(x = x.star, y = y.pred)) +
  # Confidence bounds
  geom_ribbon(data = df, aes(x = x.star, ymin = y.lower, ymax = y.upper), fill = "#0070CC", alpha = 0.3) +
  # Observed data
  geom_point(data = NULL, aes(x = x, y = y)) +
  labs(y = "$y$") + 
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot acquisition
acq.plot <- ggplot() +
  # Expected improvement
  geom_area(data = df, aes(x = x.star, y = ei), fill = '#5AB953', alpha = 0.5) +
  # Maximium point
  geom_point(data = acq.max, aes(x = x.star, y = 0), shape = 17, size = 3, col = 'red') +
  labs(x = "$x$", y = "EI") +
  scale_y_continuous(breaks = c(0,0.6), limits = c(0,0.6)) +
  theme_minimal()
  
plot5 <- plot_grid(obj.plot, acq.plot, align = "v", nrow = 2, rel_heights = c(2, 1))

tikz(file = "fig1iia.tex", width = 5.5, height = 2.5)

plot5

dev.off()


# 6 points ----------------------------------------------------------------

x <- c(des.x, acq.max$x.star)
obs <- data.frame(x = x,
                  y = f(x))
y.min <- min(y)

df <- gp.data(f, obs)

# Find acquisition maximum
acq.max <- df[which(df$ei == max(df$ei)),]

x <- obs$x
y <- obs$y

# Plot objective
obj.plot <- ggplot() +
  # True function
  geom_line(data = NULL, aes(x = x.star, y = f.star), col = "grey60", lty = 2) +
  # Prediction mean
  geom_line(data = df, aes(x = x.star, y = y.pred)) +
  # Confidence bounds
  geom_ribbon(data = df, aes(x = x.star, ymin = y.lower, ymax = y.upper), fill = "#0070CC", alpha = 0.3) +
  # Observed data
  geom_point(data = NULL, aes(x = x[1:5], y = y[1:5])) +
  # Added point 
  geom_point(data = NULL, aes(x = x[6], y = y[6]), col = 'red') +
  labs(y = "$y$") + 
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot acquisition
acq.plot <- ggplot() +
  # Expected improvement
  geom_area(data = df, aes(x = x.star, y = ei), fill = '#5AB953', alpha = 0.5) +
  # Maximium point
  geom_point(data = acq.max, aes(x = x.star, y = 0), shape = 17, size = 3, col = 'red') +
  labs(x = "$x$", y = "EI") +
  scale_y_continuous(breaks = c(0,0.5), limits = c(0,0.5)) +
  theme_minimal()

plot6 <- plot_grid(obj.plot, acq.plot, align = "v", nrow = 2, rel_heights = c(2, 1))

tikz(file = "fig1iib.tex", width = 5.5, height = 2.5)
plot6
dev.off()

