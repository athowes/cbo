library(mlrMBO) # Machine learning in R, Model Based Optimization
library(ggplot2)
set.seed(1)

# Stepwise Bayesian Optimization with mlrMBO ------------------------------

fun <- function(x) {
  x^2 + sin(2 * pi * x) * cos(0.3 * pi * x)
}

ps <- makeParamSet(
  makeNumericParam("x", lower = -3, upper = 3)
)

des <- generateDesign(n = 3, par.set = ps)
des$y <- apply(des, 1, fun)

ctrl <- makeMBOControl()
ctrl <- setMBOControlInfill(ctrl, crit = crit.ei)

opt.state <- initSMBO(
  par.set = ps, 
  design = des, 
  control = ctrl, 
  minimize = TRUE, 
  noisy = FALSE)

plot(opt.state)
prop <- proposePoints(opt.state)

y <- fun(prop$prop.points$x)

updateSMBO(opt.state, x = prop$prop.points, y = y)
plot(opt.state)

replicate(3, {
prop <- proposePoints(opt.state)
y <- fun(prop$prop.points$x)
updateSMBO(opt.state, x = prop$prop.points, y = y)
})

plot(opt.state)


# Human-in-the-loop MBO ---------------------------------------------------

# Define the search space to: over reals in [-1,2] in the first dim, and over integers in [-2,3] in the second dim
ps = makeParamSet(
  makeNumericParam("q", lower = -1, upper = 2),
  makeIntegerParam("v", lower = -2, upper = 3)
)

# Latin Hypercube Sample design, 7 points
des <- generateDesign(n = 7, par.set = ps)

# Function evaluations at these points
des$y = c(1.20, 0.97, 0.91, 3.15, 0.58, 1.12, 0.50)

ctrl <- makeMBOControl()
ctrl <- setMBOControlInfill(ctrl, crit = crit.ei)

opt.state = initSMBO(par.set = ps, design = des, control = ctrl, minimize = TRUE, noisy = FALSE)

plot(opt.state)

# What is the expected improvement suggestion?
proposePoints(opt.state)

# Don't have to stick to what it says, try this point instead
x <- data.frame(q = 1.7, v = 1)

updateSMBO(opt.state, x = x, y = 2.19)

plot(opt.state)

res <- finalizeSMBO(opt.state)
res$x
res$y

f <- function(q, v) 1 + sin(q*5) + 0.1 * (q^2 + v^2)
for (i in 1:10) {
  prop <- proposePoints(opt.state)
  x <- dfRowsToList(df = prop$prop.points, par.set = ps)
  y <- do.call(f, x[[1]])
  updateSMBO(opt.state, x = prop$prop.points, y = y)
}

plot(opt.state)
res <- finalizeSMBO(opt.state)
res$x
res$y

# mlrMBO tutorial ---------------------------------------------------------

fun = makeSingleObjectiveFunction(
  name = "SineMixture",
  fn = function(x) sin(x[1])*cos(x[2])/2 + 0.04 * sum(x^2),
  par.set = makeNumericParamSet(id = "x", len = 2, lower = -5, upper = 5)
)

library(plot3D)
plot3D(fun, contour = TRUE, lightning = TRUE)

# In this simple example we construct the control object with the defaults:
ctrl = makeMBOControl()
# For this numeric optimization we are going to use the Expected Improvement as infill criterion:
ctrl = setMBOControlInfill(ctrl, crit = crit.ei)
# We will allow for exactly 25 evaluations of the objective function:
ctrl = setMBOControlTermination(ctrl, max.evals = 25L)

des = generateDesign(n = 8L, par.set = getParamSet(fun), fun = lhs::randomLHS)
autoplot(fun, render.levels = TRUE) + geom_point(data = des)

sur.lrn = makeLearner("regr.km", predict.type = "se", config = list(show.learner.output = FALSE))

res = mbo(fun = fun, design = des, learner = sur.lrn, control = ctrl, show.info = TRUE)

opdf = as.data.frame(res$opt.path)
autoplot(fun, render.levels = TRUE, render.contours = FALSE) + geom_text(data = opdf, aes(label = dob))

res$x
res$y
