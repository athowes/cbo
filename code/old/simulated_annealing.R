# Simulated Annealing

georgios <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  x1^2 + 50*(1+x1^2)^2 * (x2-sin(x1))^2
}

simulated_annealing <- function(func, s0, niter = 10, step = 0.1) {
  
  # Initialize
  ## s stands for state
  ## f stands for function value
  ## b stands for best
  ## c stands for current
  ## n stands for neighbor
  s_b <- s_c <- s_n <- s0
  f_b <- f_c <- f_n <- func(s0)
  message("It\tBest\tCurrent\tNeigh\tTemp")
  message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
  
  for (k in 1:niter) {
    Temp <- (1- step)^k
    # consider a random neighbor
    s_n <- rnorm(2, s_c, 0.1)
    f_n <- func(s_n)
    # update current state
    if (f_n < f_c | runif(1, 0, 1) < exp((f_c - f_n) / Temp)) {
      s_c <- s_n 
      f_c <- f_n
    }
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n
    }
    message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", k, f_b, f_c, f_n, Temp))
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b))
}     

simulated_annealing(georgios, s0 = c(1,1), niter = 1000, step = 0.01)

models <- vector("list", 10)

for (i in 1:30) {
  models[[i]] <- simulated_annealing(georgios, s0 = rnorm(2, 0, 1), niter = 1000, step = 0.01)
}

vals <- c(rep(0, 30))

for (i in 1:30) vals[i] <- models[[i]]$best_value

boxplot(vals)

models[[1]]
