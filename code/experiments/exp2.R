# Experiment 2: Constrained Optimization
# EI vs LHS

library(tidyverse)
library(mlrMBO)
library(laGP)
library(lhs)
library(tikzDevice)
library(tidyverse)
library(viridis)
library(ggpolypath)
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

disk <- function(x) {
  
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  
  x1 <- x[,1]
  x2 <- x[,2]
  
  y <- 2/9 - (x1 - 0.5)^2 - (x2 - 0.5)^2
  return(y)
}

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

p.disk <- makeSingleObjectiveFunction(
  fn = function(xx) {
    x1 <- xx[1]
    x2 <- xx[2]
    
    y <- 2/9 - (x1 - 0.5)^2 - (x2 - 0.5)^2
    return(y)
  },
  par.set = makeParamSet(
    makeNumericParam("x1", lower = 0, upper = 1),
    makeNumericParam("x2", lower = 0, upper = 1)
  )
) # Constraint function

# EIC

eps <- sqrt(.Machine$double.eps)

EIC <- function(gpi, gpc, x, fmin, pred=predGPsep) {
  # Gives the unconstrained expected improvement statistic
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  p <- pred(gpi, x, lite=TRUE)
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d/sigma
  ei <- d*pnorm(dn) + sigma*dnorm(dn)
  
  # Now calculate the probability of satisfying the constraints
  pc <- pred(gpc, x, lite=TRUE)
  sigmac <- sqrt(pc$s2)
  pc <- pnorm(pc$mean/sigmac)
  
  EIC <- ei*pc
  return(EIC)
}

obj.EIC <- function(x, fmin, gpi, gpc) {
  - EIC(gpi, gpc, x, fmin)
}

EIC.search <- function(X, y, gpi, gpc, multi.start = 5) {
  m <- which.min(y)
  fmin <- y[m] # Calculate the incumbent value
  start <- rbind(randomLHS(multi.start, ncol(X)))
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    out <- optim(start[i,], obj.EIC, method="L-BFGS-B", 
                 lower=0, upper=1, gpi=gpi, gpc=gpc, fmin=fmin)
    xnew[i,] <- c(out$par, -out$value)
  }

  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "x1", "x2", "val")
  return(solns)
}

optim.EIC <- function(f, c, ninit, stop) {
  X <- randomLHS(ninit, 2); y <- f(X); yc <- c(X)
  
  gpi <- newGPsep(X, y, d=0.2, g=1e-7, dK=TRUE)
  da <- darg(list(mle=TRUE, max=0.5), X)
  mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)$msg
  
  gpc <- newGPsep(X, yc, d=0.2, g=1e-7, dK=TRUE)
  dac <- darg(list(mle=TRUE, max=0.5), X)
  mleGPsep(gpc, param="d", tmin=dac$min, tmax=dac$max, ab=dac$ab)$msg
  
  maxEIC <- numeric(0)
  for(i in (ninit+1):stop) {
    solns <- EIC.search(X, y, gpi, gpc)
    m <- which.max(solns$val)
    maxEIC <- append(maxEIC, solns$val[m])
    
    xnew <- as.matrix(solns[m,3:4])
    ynew <- f(xnew)
    ycnew <- c(xnew)
    
    updateGPsep(gpi, matrix(xnew, nrow=1), ynew)
    mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
    
    updateGPsep(gpc, matrix(xnew, nrow=1), ycnew)
    mlec <- mleGPsep(gpc, param="d", tmin=dac$min, tmax=dac$max, ab=dac$ab)
    
    X <- rbind(X, xnew); y <- append(y, ynew); yc <- append(yc, ycnew)
  }
  
  return(list(X=X, y=y, yc=yc, maxEIC=maxEIC))
}

# IECI

ECI <- function(gpi = gpi, gpc = gpc, x, xref, fmin=fmin, pred=predGPsep) {
  # Expected conditional improvement statistic
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  if(is.null(nrow(xref))) xref <- matrix(xref, nrow=1)
  
  # Calculate the probability of satisfying the constraints at the reference location
  pc <- pred(gpc, xref, lite=TRUE)
  sigmac <- sqrt(pc$s2)
  pc <- pnorm(pc$mean/sigmac)
  
  # Calculate the expected conditional improvement at xref given that xcand is to be added to the design
  # Weighted by the probabillity of satisfying the constraints at the reference location
  eci <- ieciGPsep(gpi, x, fmin, Xref = xref, w = pc)
  
  return(eci)
}

grid <- expand.grid(seq(0,1,0.1), seq(0,1,0.1))

IECI <- function(xrefgrid, gpi, gpc, x, fmin) {
  ECI_grid <- apply(grid, MARGIN = 1, FUN = ECI, gpi=gpi, gpc=gpc, x=x, fmin=fmin)
  return(sum(ECI_grid))
}

IECI.search <- function(X, y, xrefgrid, gpi, gpc, multi.start = 5) {
  m <- which.min(y)
  fmin <- y[m] # Calculate the incumbent value
  start <- rbind(randomLHS(multi.start, ncol(X)))
  xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
  for(i in 1:nrow(start)) {
    out <- optim(start[i,], IECI, method="L-BFGS-B", 
                 lower=0, upper=1, 
                 xrefgrid=xrefgrid, gpi=gpi, gpc=gpc, fmin=fmin)
    xnew[i,] <- c(out$par, -out$value)
  }
  
  solns <- data.frame(cbind(start, xnew))
  names(solns) <- c("s1", "s2", "x1", "x2", "val")
  return(solns)
}

optim.IECI <- function(f, c, xrefgrid, ninit, stop) {
  X <- randomLHS(ninit, 2); y <- f(X); yc <- c(X)
  
  gpi <- newGPsep(X, y, d=0.2, g=1e-7, dK=TRUE)
  da <- darg(list(mle=TRUE, max=0.5), X)
  mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)$msg
  
  gpc <- newGPsep(X, yc, d=0.2, g=1e-7, dK=TRUE)
  dac <- darg(list(mle=TRUE, max=0.5), X)
  mleGPsep(gpc, param="d", tmin=dac$min, tmax=dac$max, ab=dac$ab)$msg
  
  maxeici <- numeric(0)
  for(i in (ninit+1):stop) {
    solns <- EICI.search(X, y, xrefgrid, gpi, gpc)
    m <- which.max(solns$val)
    maxeici <- append(maxeici, solns$val[m])
    
    xnew <- as.matrix(solns[m,3:4])
    ynew <- f(xnew)
    ycnew <- c(xnew)
    
    updateGPsep(gpi, matrix(xnew, nrow=1), ynew)
    mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
    
    updateGPsep(gpc, matrix(xnew, nrow=1), ycnew)
    mlec <- mleGPsep(gpc, param="d", tmin=dac$min, tmax=dac$max, ab=dac$ab)
    
    X <- rbind(X, xnew); y <- append(y, ynew); yc <- append(yc, ycnew)
  }
  
  return(list(X=X, y=y, yc=yc, maxeici=maxeici))
}

test <- optim.IECI(f = branin, c = disk, xrefgrid = grid, ninit = 5, stop = 20)
rownames(test$X) <- c()
test.df <- as.data.frame(cbind(test$X, y = test$y, yc = test$yc))
test.df <- mutate(test.df, 
             feasible = (yc > 0),
             best = ifelse(feasible, y, NA),
             id = row_number())

autoplot(p.branin, render.levels = TRUE, render.contours = FALSE) + theme_minimal() + 
  geom_text(data = test.df, aes(x = x1, y = x2, label = id), nudge_x = 0.02, nudge_y = 0.02, col = "grey20", size = 3) +
  geom_point(data = test.df, aes(x = x1, y = x2), shape = 21, fill = "white") +
  scale_fill_viridis(direction = -1)


# fig5iii and fig5iv

set.seed(19)

os <- optim.EIC(f = branin, c = disk, 5, 20)

df <- as.data.frame(cbind(os$X, y = os$y, yc = os$yc))
df <- mutate(df, 
             feasible = (yc > 0),
             best = ifelse(feasible, y, NA),
             id = row_number())

circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Posterior mean plots
X.data <- df[,c('x1','x2')]
y.data <- df[,c('y')]
yc.data <- df[,c('yc')]

# Objective function
gpi.plot <- newGPsep(X.data, y.data, d=0.2, g=1e-7, dK=TRUE)
da.plot <- darg(list(mle=TRUE, max=0.5), X)
mleGPsep(gpi.plot, param="d", tmin=da.plot$min, tmax=da.plot$max, ab=da.plot$ab)$msg

obj.mean <- function(x, gpi.plot) {
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  predGPsep(gpi.plot, x, lite=TRUE)$mean
}

obj.mean.smoof <- makeSingleObjectiveFunction(
  fn = function(x) {
    return(obj.mean(x, gpi.plot))
  },
  par.set = makeParamSet(
    makeNumericParam("x1", lower = 0, upper = 1),
    makeNumericParam("x2", lower = 0, upper = 1)
  )
)

plot2 <- autoplot(obj.mean.smoof, render.levels = TRUE, render.contours = TRUE) +
  theme_minimal() + 
  labs(x = "", y = "", fill = "$y$", subtitle = "(b) surrogate objective") +
  scale_fill_gradientn(colors = viridis_pal(option = "viridis", direction = -1)(9), limits=c(-1.77, 4.88), 
                       breaks = c(0, 2, 4), labels = c("0 ", "2 ", "4 ")) +
  theme(legend.key.size = unit(0.3, "cm"))

plot1 <- autoplot(p.branin, render.levels = TRUE, render.contours = TRUE) +
  theme_minimal() + 
  labs(x = "", y = "", fill = "$y$", subtitle = "(a) true objective") +
  guides(fill=FALSE) +
  scale_fill_gradientn(colors = rev(viridis_pal()(9)), limits=c(-1.77, 4.88))

# Constraint function
gpc.plot <- newGPsep(X.data, yc.data, d=0.2, g=1e-7, dK=TRUE)
dac.plot <- darg(list(mle=TRUE, max=0.5), X)
mleGPsep(gpc.plot, param="d", tmin=dac.plot$min, tmax=dac.plot$max, ab=dac.plot$ab)$msg

con.mean <- function(x, gpc.plot) {
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  predGPsep(gpc.plot, x, lite=TRUE)$mean
}

con.mean.smoof <- makeSingleObjectiveFunction(
  fn = function(x) {
    return(con.mean(x, gpc.plot))
  },
  par.set = makeParamSet(
    makeNumericParam("x1", lower = 0, upper = 1),
    makeNumericParam("x2", lower = 0, upper = 1)
  )
) # Constraint function

plot4 <- autoplot(con.mean.smoof, render.levels = TRUE, render.contours = TRUE) + theme_minimal() + 
  labs(x = "", y = "", fill = "$y$", subtitle = "(d) surrogate constraint") +
  scale_fill_gradientn(colors = viridis_pal()(9), limits=c(-0.3, 0.3), breaks = c(-0.2, 0, 0.2), labels = c("-.2", "0", ".2")) +
  theme(legend.key.size = unit(0.3, "cm"))

plot3 <- autoplot(p.disk, render.levels = TRUE, render.contours = TRUE) + theme_minimal() + 
  labs(x = "", y = "", fill = "$y$", subtitle = "(c) true constraint") +
  guides(fill=FALSE) +
  scale_fill_gradientn(colors = viridis_pal()(9), limits=c(-0.3, 0.3))

circle <- circleFun(c(0.5, 0.5), 2*sqrt(2)/3, npoints = 100)
square <- data.frame(x = c(-0.01,-0.01,1.01,1.01,-0.01), y = c(-0.01,1.01,1.01,-0.01,-0.01))
shape <- data.frame(c(rbind(circle, square)))
shape <- mutate(shape, id = c(rep(1L, 100), rep(2L, 5)))

which.min(df$y)
df[15,]

fig5iii <- autoplot(p.branin, render.levels = TRUE, render.contours = FALSE) + theme_minimal() + 
  geom_polypath(data = shape, aes(x, y), fill = "grey30", alpha = 0.5, rule = "evenodd") +
  labs(x = "$x_1$", y = "$x_2$", fill = "$y$") +
  geom_text(data = df, aes(x = x1, y = x2, label = id), nudge_x = 0.02, nudge_y = 0.02, col = "grey20", size = 3) +
  geom_point(data = df[-15,], aes(x = x1, y = x2), shape = 21, fill = "white") +
  geom_point(data = df[15,], aes(x = x1, y = x2), shape = 21, fill = "grey70") +
  scale_fill_viridis(direction = -1)

tikz(file = "fig5iii.tex", width = 5.5, height = 5)
fig5iii
dev.off()

bottom_row <- plot_grid(plot3, plot4, nrow = 1, rel_widths = c(1, 1.25))
top_row <- plot_grid(plot1, plot2, nrow = 1, rel_widths = c(1, 1.25))

tikz(file = "fig5iva.tex", width = 5, height = 2)
top_row
dev.off()

tikz(file = "fig5ivb.tex", width = 5, height = 2)
bottom_row
dev.off()

# fig5v and vi

# Expected Improvement Constrained

reps <- 50
prog.EIC <- matrix(NA, nrow = 20, ncol = reps)

for(r in 1:reps) {
  os <- optim.EIC(f = branin, c = disk, 5, 20)
  df <- as.data.frame(cbind(os$X, y = os$y, yc = os$yc))
  df <- mutate(df, 
               feasible = (yc > 0),
               yfeasible = ifelse(feasible, y, 2),
               best = cummin(yfeasible)
  )
  prog.EIC[,r] <- df$best
}

avg.prog.EIC <- cbind(rowMeans(prog.EIC), rowVars(prog.EIC))
avg.prog.EIC <- as.data.frame(avg.prog.EIC) %>% mutate(id = row_number(), variable = "EIC")
names(avg.prog.EIC) <- c("avg", "var", "id", "variable")

# Latin hyper-cube

reps <- 50
prog.lhs <- matrix(NA, nrow = 20, ncol = reps)

for(r in 1:reps) {
  lhs.des <- generateDesign(par.set = getParamSet(p.branin), fun = lhs::randomLHS, n = 20L)
  lhs.des$y <- apply(lhs.des, 1, p.branin)
  lhs.des$yc <- apply(lhs.des, 1, p.disk)
  df <- as.data.frame(lhs.des)
  df <- mutate(df, 
               feasible = (yc > 0),
               yfeasible = ifelse(feasible, y, 2),
               best = cummin(yfeasible)
  )
  prog.lhs[,r] <- df$best
}

avg.prog.lhs <- cbind(rowMeans(prog.lhs), rowVars(prog.lhs))
avg.prog.lhs <- as.data.frame(avg.prog.lhs) %>% mutate(id = row_number(), variable = "LHS")
names(avg.prog.lhs) <- c("avg", "var", "id", "variable")

# IECI

run_experiment <- function(reps, grid, seed) {
  set.seed(seed)
  prog.ieci <- matrix(NA, nrow = 20, ncol = reps)
  for(r in 1:reps) {
    os <- optim.IECI(f = branin, c = disk, xrefgrid = grid, 5, 20)
    df <- as.data.frame(cbind(os$X, y = os$y, yc = os$yc))
    df <- mutate(df, 
                 feasible = (yc > 0),
                 yfeasible = ifelse(feasible, y, 2),
                 best = cummin(yfeasible)
    )
    prog.ieci[,r] <- df$best
  }
  return(prog.ieci)
}

prog.ieci1 <- run_experiment(reps = 10, grid = expand.grid(seq(0,1,0.05), seq(0,1,0.05)), seed = 1)
prog.ieci2 <- run_experiment(reps = 10, grid = expand.grid(seq(0,1,0.05), seq(0,1,0.05)), seed = 2)
prog.ieci3 <- run_experiment(reps = 10, grid = expand.grid(seq(0,1,0.05), seq(0,1,0.05)), seed = 3)
prog.ieci4 <- run_experiment(reps = 10, grid = expand.grid(seq(0,1,0.05), seq(0,1,0.05)), seed = 4)
prog.ieci5 <- run_experiment(reps = 10, grid = expand.grid(seq(0,1,0.05), seq(0,1,0.05)), seed = 5)

prog.ieci <- cbind(prog.ieci1, prog.ieci2, prog.ieci3, prog.ieci4, prog.ieci5)

avg.prog.ieci <- cbind(rowMeans(prog.ieci), rowVars(prog.ieci))
avg.prog.ieci <- as.data.frame(avg.prog.ieci) %>% mutate(id = row_number(), variable = "IECI")
names(avg.prog.ieci) <- c("avg", "var", "id", "variable")

# Lineplot

avg.prog <- rbind(avg.prog.EIC, avg.prog.ieci, avg.prog.lhs)

library(RColorBrewer)
mycols <- brewer.pal(4, "Set1")[-1]

lineplot <- ggplot(new.avg.prog, aes(x = id, y = avg, col = variable)) +
  geom_line() +
  labs(x = "Number of samples", y = "Best feasible evaluation") +
  coord_cartesian(ylim = c(-1.2, 0.1)) +
  geom_hline(yintercept = -1.047, linetype="dashed") +
  labs(col = "") +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.81)) +
  scale_color_manual(values = rev(mycols))

tikz(file = "fig5v.tex", width = 5, height = 2.5)

lineplot

dev.off()

# Boxplot

prog.lhs.tall <- as.data.frame(prog.lhs) %>% mutate(id = row_number())
prog.lhs.tall <- melt(prog.lhs.tall, id = "id")

prog.EIC.tall <- as.data.frame(prog.EIC) %>% mutate(id = row_number())
prog.EIC.tall <- melt(prog.EIC.tall, id = "id")

prog.ieci.tall <- as.data.frame(prog.ieci) %>% mutate(id = row_number())
prog.ieci.tall <- melt(prog.ieci.tall, id = "id")

final.best <- as.data.frame(
  cbind(
    filter(prog.lhs.tall, id == 20)[,3],
    filter(prog.ieci.tall, id == 20)[,3],
    filter(prog.EIC.tall, id == 20)[,3]
  )
)

names(final.best) <- c("LHS", "IECI", "EIC")
final.best <- melt(final.best)

boxplot <- ggplot(final.best, aes(x = variable, y = value, col = variable)) +
  geom_boxplot(alpha = 0.7, fill = "grey95") +
  theme_minimal() +
  labs(x = "", y = "Best feasible evaluation") +
  scale_color_manual(values = mycols) +
  guides(col = FALSE)

tikz(file = "fig5vi.tex", width = 5, height = 3)

boxplot

dev.off()