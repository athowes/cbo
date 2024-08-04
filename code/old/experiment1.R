# Experiment 1: Branin function

library(tidyverse)
library(reshape2)
library(mlrMBO)
library(DiceKriging)

ps <- makeParamSet(
  makeNumericParam("x1", lower = 0, upper = 1),
  makeNumericParam("x2", lower = 0, upper = 1)
)

branin <- makeSingleObjectiveFunction(
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
  par.set = ps
) # Branin function (rescaled)

autoplot(branin, render.levels = TRUE, render.contours = FALSE) + theme_minimal() + labs(x = "$x_1$", y = "$x_2$", fill = "$y$")

# Bayesian Optimization ---------------------------------------------------

ctrl <- makeMBOControl()
ctrl <- setMBOControlInfill(ctrl, crit = crit.ei)
ctrl <- setMBOControlTermination(ctrl, max.evals = 20L)
sur.lrn <- makeLearner("regr.km", predict.type = "se", covtype = "matern5_2", config = list(show.learner.output = TRUE))

reps <- 100
prog <- matrix(NA, nrow = 20, ncol = reps)
for(r in 1:reps) {
  des <- generateDesign(n = 5L, par.set = getParamSet(branin), fun = lhs::randomLHS)
  res <- mbo(fun = branin, design = des, learner = sur.lrn, control = ctrl, show.info = TRUE)
  opdf <- as.data.frame(res$opt.path)
  opdf <- mutate(opdf, best = cummin(y))
  prog[,r] <- opdf$best
}

prog_bo <- as.data.frame(prog) %>% mutate(id = row_number())
prog_bo <- melt(prog_bo, id = "id")

avg_prog_bo <- rowMeans(prog)
avg_prog_bo <- as.data.frame(avg_prog_bo) %>% mutate(id = row_number())

ggplot() +
  geom_line(data = prog_bo, aes(x = id, y = value, group = variable), alpha = 0.1) +
  geom_line(data = avg_prog_bo, aes(x = id, y = avg_prog_bo), size = 1) +
  theme_minimal()



  

