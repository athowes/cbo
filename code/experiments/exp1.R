# Experiment 1: Unconstrained Optimization
# EI vs LHS

library(tidyverse)
library(reshape2)
library(mlrMBO)
library(laGP)
library(lhs)
library(tikzDevice)
library(viridis)
library(Rfast)
library(cowplot)

setwd("C:/Users/Adam/Desktop/bgo/code")

branin <- function(x) {
  
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  
  x1 <- x[,1]
  x2 <- x[,2]
  
  x1bar <- 15*x1 - 5
  x2bar <- 15 * x2
  
  term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
  term2 <- (10 - 10/(8*pi)) * cos(x1bar)
  
  y <- (term1^2 + term2 - 44.81) / 51.95
  return(y)
}

# Expected Improvement 

eps <- sqrt(.Machine$double.eps)

EI <- function(gpi, x, fmin, pred=predGPsep) {
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  p <- pred(gpi, x, lite=TRUE)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d/sigma
  ei <- d*pnorm(dn) + sigma*dnorm(dn)
  return(ei)
}

obj.EI <- function(x, fmin, gpi) {
  - EI(gpi, x, fmin)
}

EI.search <- function(X, y, gpi, multi.start = 5) {
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(X[m,], nrow=1)
  if(multi.start > 1) 
    start <- rbind(start, randomLHS(multi.start-1, ncol(X)))
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    if(EI(gpi, start[i,], fmin) <= eps)
    { out <- list(value=-Inf); next }
    out <- optim(start[i,], obj.EI, method="L-BFGS-B", 
                 lower=0, upper=1, gpi=gpi, fmin=fmin)
    xnew[i,] <- c(out$par, -out$value)
  }
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "x1", "x2", "val")
  solns <- solns[(solns$val > sqrt(.Machine$double.eps)),]
  return(solns)
}

optim.EI <- function(f, ninit, stop) {
  X <- randomLHS(ninit, 2); y <- f(X)
  gpi <- newGPsep(X, y, d=0.1, g=1e-7, dK=TRUE)
  da <- darg(list(mle=TRUE, max=0.5), X)
  mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)$msg
  maxei <- c()
  for(i in (ninit+1):stop) {
    solns <- EI.search(X, y, gpi)
    m <- which.max(solns$val)
    maxei <- c(maxei, solns$val[m])
    xnew <- as.matrix(solns[m,3:4])
    ynew <- f(xnew)
    updateGPsep(gpi, matrix(xnew, nrow=1), ynew)
    mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
    X <- rbind(X, xnew); y <- c(y, ynew)
  }
  deleteGPsep(gpi)
  
  return(list(X=X, y=y, maxei=maxei))
}

# fig5i

set.seed(3)

p.branin <- makeSingleObjectiveFunction(
  fn = function(xx) {
    x1 <- xx[1]
    x2 <- xx[2]
    
    x1bar <- 15*x1 - 5
    x2bar <- 15 * x2
    
    term1 <- x2bar - 5.1*x1bar^2/(4*pi^2) + 5*x1bar/pi - 6
    term2 <- (10 - 10/(8*pi)) * cos(x1bar)
    
    y <- (term1^2 + term2 - 44.81) / 51.95
    return(y)
  },
  par.set = makeParamSet(
    makeNumericParam("x1", lower = 0, upper = 1),
    makeNumericParam("x2", lower = 0, upper = 1)
  )
) # Branin function (rescaled)

sample <- optim.EI(branin, 5, 20)
which.min(sample$y)
round(sample$y[19], 3)
round(sample$X[19,], 3)

rownames(sample$X) <- c()
sample.loc <- as.data.frame(sample$X) %>% mutate(id = row_number())

tikz(file = "fig5i.tex", width = 5.5, height = 5)

autoplot(p.branin, render.levels = TRUE, render.contours = FALSE) + theme_minimal() + 
  labs(x = "$x_1$", y = "$x_2$", fill = "$y$") +
  geom_text(data = sample.loc, aes(label = id), nudge_x = 0.02, nudge_y = 0.02, col = "grey20", size = 3) +
  geom_point(data = sample.loc[1:19,], shape = 21, fill = "white") +
  geom_point(data = sample.loc[20, ], shape = 21, fill = "grey70") +
  scale_fill_viridis(direction = -1)

dev.off()

# fig5ii

reps <- 50

prog.lhs <- matrix(NA, nrow = 20, ncol = reps)
for(r in 1:reps) {
  lhs.des <- generateDesign(par.set = getParamSet(p.branin), fun = lhs::randomLHS, n = 20L)
  lhs.des$y <- apply(lhs.des, 1, p.branin)
  df <- as.data.frame(lhs.des)
  df <- mutate(df, best = cummin(y))
  prog.lhs[,r] <- df$best
}

avg.prog.lhs <- cbind(rowMeans(prog.lhs), rowVars(prog.lhs))
avg.prog.lhs <- as.data.frame(avg.prog.lhs) %>% mutate(id = row_number(), variable = "LHS")
names(avg.prog.lhs) <- c("avg", "var", "id", "variable")

prog.ei <- matrix(NA, nrow = 20, ncol = reps)
for(r in 1:reps) {
  os <- optim.EI(branin, 5, 20)
  df <- cbind(os$X, y = os$y)
  rownames(df) <- c()
  df <- as.data.frame(df)
  df <- mutate(df, best = cummin(y))
  prog.ei[,r] <- df$best
}

avg.prog.ei <- cbind(rowMeans(prog.ei), rowVars(prog.ei))
avg.prog.ei <- as.data.frame(avg.prog.ei) %>% mutate(id = row_number(), variable = "EI")
names(avg.prog.ei) <- c("avg", "var", "id", "variable")

avg.prog <- rbind(avg.prog.ei, avg.prog.lhs)

lineplot <- ggplot(avg.prog, aes(x = id, y = avg, col = variable)) +
  geom_line() +
  labs(x = "Number of samples", y = "Best evaluation") +
  coord_cartesian(ylim = c(-1.2, 0.1)) +
  geom_hline(yintercept = -1.047, linetype="dashed") +
  labs(col = "") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.81)) +
  scale_color_brewer(palette = "Set1")

prog.lhs.tall <- as.data.frame(prog.lhs) %>% mutate(id = row_number())
prog.lhs.tall <- melt(prog.lhs.tall, id = "id")

prog.ei.tall <- as.data.frame(prog.ei) %>% mutate(id = row_number())
prog.ei.tall <- melt(prog.ei.tall, id = "id")

final.best <- as.data.frame(
  cbind(
  filter(prog.ei.tall, id == 20)[,3],
  filter(prog.lhs.tall, id == 20)[,3] 
  )
)

names(final.best) <- c("EI", "LHS")
round(final.best, 3)
final.best <- melt(final.best)

sum(round(final.best$EI, 3) == -1.047)
sum(round(final.best$LHS, 3) == -1.047)

boxplot <- ggplot(final.best, aes(x = variable, y = value, col = variable)) +
  geom_boxplot(alpha = 0.7, fill = "grey100") +
  theme_minimal() +
  labs(x = "", y = "Best evaluation", subtitle = "(b)") +
  guides(col = FALSE)

tikz(file = "fig5ii.tex", width = 5, height = 2.5)
lineplot
dev.off()
