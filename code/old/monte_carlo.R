# Monte Carlo tests 
# Adam Howes

# A simple Monte Carlo example --------------------------------------------

rm(list = ls())

draws <- seq(from = 100, to = 10000, by = 100)
area <- c()

for (i in 1:length(draws)) {
  int.store <- rnorm(draws[i], mean = 0, sd = 1)
  area[i] <- sum(int.store > -1.96 & int.store < 1.96)/length(int.store)
}

true <- integrate(f = dnorm, lower = -1.96, upper = 1.96, mean = 0, sd = 1)$value

plot(area~draws)
abline(h = true, lty = 2)


# Accept-Reject sampling, tent example --------------------------------------------------

tent <- function(x, mode = 0, width = 2){
  
  max.height <- 1/width
  lo.limit <- mode - width
  hi.limit <- mode + width
  out <- c()
  
  for (i in 1:length(x)){
    
    if(x[i] > lo.limit & x[i] < hi.limit){
      out[i] <- min(max.height*((x[i]-lo.limit)/width), max.height*((hi.limit - x[i])/width))
    }
    else {
      out[i] <- 0
    }
  }
  
  return(out)
}

points <- seq(from = -3, to = 3, by = 0.1)

den <- tent(points)

plot(den, type = "l", ylim = c(0, 0.6))

prop.den <- dnorm(points, mean = 0, sd = 1)
lines(1.3*prop.den, lty = 2)

reps <- 40000
proposal <- rnorm(reps, mean = 0, sd = 1)
decision <- tent(proposal)/(1.3*dnorm(proposal, mean = 0, sd = 1))
keep <- runif(reps) <= decision

plot(density(proposal[keep], ylim = c(0, 0.6)))

# Metropolis-Hastings, random walk proposal -------------------------------

rm(list = ls())

cauchy <- function(x, x0 = 0, gamma = 1) {
  1/(pi*gamma*(1+((x-x0)/gamma)^2))
}

plot(cauchy, xlim = c(-10, 10))

reps <- 50000

chain <- c(0)

for (i in 1:reps){
  proposal <- chain[i] + runif(1, min = -1, max = 1)
  accept <- runif(1) < cauchy(proposal)/cauchy(chain[i])
  chain[i+1] <- ifelse(accept, proposal, chain[i])
}

plot(density(chain[10000:50000]), ylim = c(0, 0.4))

den <- cauchy(seq(from = -10, to = 10, by = 0.1), x = 0, gamma = 1)
lines(den~seq(from = -10, to = 10, by = 0.1), lty = 2, col = 2)

plot(chain, type = "l")
      
# Metropolis-Hastings, harder case ----------------------------------------

rm(list = ls())

target <- function(x) {
  exp(-x[1]^2 - 50*(1+x[1]^2)^2*(x[2]-sin(x[1]))^2)
}

reps <- 500
chain <- matrix(data = NA, ncol = 2, nrow = reps + 1)
chain[1,] <- c(0,0)

for (i in 1:reps){
  proposal <- chain[i,] + runif(2, min = -1, max = 1)
  accept <- runif(1) < target(proposal)/target(chain[i,])
  chain[i+1,] <- ifelse(accept, proposal, chain[i,])
}

par(mfrow = c(1,1))                                               
plot(chain[,1], type = "l")  
plot(chain[,2], type = "l")