# Fig 1i

set.seed(1)
setwd("C:/Users/Adam/Desktop/bgo/code")

library(mlrMBO)
library(plot3D)
library(tikzDevice)
library(tidyverse)
library(viridis)

fun = makeSingleObjectiveFunction(
  fn = function(x) {
    
    x1 <- x[1]
    x2 <- x[2]
    
    term1 <- (4-2.1*x1^2+0.3*x1^4) * x1^2
    term2 <- x1*x2
    term3 <- (-4+4*x2^2) * x2^2
    
    y <- term1 + term2 + term3 + 0.6*x2
    return(y)
  },
  par.set = makeNumericParamSet(id = "x", len = 2, lower = -1, upper = 1)
) # Camel function

# Grid search
grid.des = generateGridDesign(par.set = getParamSet(fun), resolution = 5)
grid.des$y = apply(grid.des, 1, fun)
grid.max <- grid.des[which.min(grid.des$y),]

# Random search
random.des = generateRandomDesign(par.set = getParamSet(fun), n = 25L)
random.des$y = apply(random.des, 1, fun)
random.max <- random.des[which.min(random.des$y),]

# Latin-hypercube search
lhs.des = generateDesign(par.set = getParamSet(fun), fun = lhs::randomLHS, n = 25L)
lhs.des$y = apply(lhs.des, 1, fun)
lhs.max <- lhs.des[which.min(lhs.des$y),]

tikz(file = "fig1i.tex", width = 5.5, height = 5)
plot <- autoplot(fun, render.levels = TRUE, render.contours = FALSE) +
  scale_fill_viridis(direction = -1) +
  geom_point(data = grid.des, shape = 21, fill = "white") +
  geom_point(data = grid.max, shape = 21, fill = "grey70") +
  geom_point(data = random.des, shape = 24, fill = "white") +
  geom_point(data = random.max, shape = 24, fill = "grey70") +
  geom_point(data = lhs.des, shape = 22, fill = "white") +
  geom_point(data = lhs.max, shape = 22, fill = "grey70") +
  labs( x = "$x_1$", y = "$x_2$", fill = "$y$") +
  theme_minimal()
plot
dev.off()
