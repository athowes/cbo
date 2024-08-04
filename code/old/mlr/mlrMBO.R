library(mlrMBO) # Machine learning in R, Model Based Optimization
library(ggplot2)

obj.fun <- makeAlpine02Function(1)
autoplot(obj.fun) + theme_minimal()

surr.km = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", control = list(trace = FALSE))

des = generateDesign(n = 3, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)

control = makeMBOControl()
control = setMBOControlTermination(control, iters = 4)
control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())

run = exampleRun(obj.fun, learner = surr.km, design = des, control = control, show.info = FALSE)

plotExampleRun(run, iters = c(1L, 2L, 3L, 4L), pause = FALSE) + theme_minimal()
