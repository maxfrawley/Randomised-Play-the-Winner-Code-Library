# Non centrality parameter function and plot for fixed allocation target
phi <- function(pA, pB){
  qA <- 1 - pA
  qB <- 1 - pB
  rho <- qB/(qA + qB)
  phi <- (pA - pB)^2/((pA*qA)/rho + (pB*qB)/(1-rho))
  return(phi)
}
pB_values <- c(0.1, 0.3, 0.5, 0.7)
pA_values <- seq(0, 0.95, by = 0.01)
phi_values <- sapply(pA_values[pA_values >= pB_values[1]], function(x) phi(x, pB_values[1]))
plot(pA_values[pA_values >= pB_values[1]], phi_values, type = "l", lty = 1, col = 1, ylim = range(0, phi_values), 
     xlab = "pA", ylab = "Noncentrality Parameter")
for (i in 2:length(pB_values)) {
  phi_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) phi(x, pB_values[i]))
  lines(pA_values[pA_values >= pB_values[i]], phi_values, lty = i, col = i, ylim = range(0, phi_values))
}
legend("topleft", legend = paste("pB =", pB_values), lty = 1:length(pB_values), col = 1:length(pB_values))


# Expected proportion of treatment failures function and plot for fixed allocation target
failures <- function(pA, pB){
  qA <- 1 - pA
  qB <- 1 - pB
  rho <- qB/(qA + qB)
  fail <- rho*qA + (1-rho)*qB
  return(fail)
}
pB_values <- c(0.1, 0.3, 0.5, 0.7)
pA_values <- seq(0, 1, by = 0.01)
fail_values <- sapply(pA_values[pA_values >= pB_values[1]], function(x) failures(x, pB_values[1]))
plot(pA_values[pA_values >= pB_values[1]], fail_values, type = "l", lty = 1, col = 1, ylim = range(0, fail_values), 
     xlab = "pA", ylab = "Expected Failure Rate")

for (i in 2:length(pB_values)) {
  fail_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) failures(x, pB_values[i]))
  lines(pA_values[pA_values >= pB_values[i]], fail_values, lty = i, col = i, ylim = range(0, fail_values))
}
legend("topright", legend = paste("pB =", pB_values), lty = 1:length(pB_values), col = 1:length(pB_values))

# Asymptotic Variance of treatment proportions and Cramer-Rao Lower Bound -  Function and Figure
# Asymptotic Variance
asym_variance <- function(pA, pB){
  qA <- 1 - pA
  qB <- 1 - pB
  v <- (qA*qB*(5-2*(qA + qB))) / ((2*(qA+qB)- 1)*(qA + qB)^2)
  return(v)
}
pB_values <- c(0.1, 0.3, 0.4, 0.48)
pA_values <- seq(0, 1, by = 0.01)
var_values <- sapply(pA_values[pA_values >= pB_values[1]], function(x) asym_variance(x, pB_values[1]))
plot(pA_values[pA_values >= pB_values[1]], var_values, type = "l", lty = 1, ylim = range(0, 2.5),  
     xlab = "pA", ylab = "Var(NA(n))")

for (i in 2:length(pB_values)) {
  var_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) asym_variance(x, pB_values[i]))
  lines(pA_values[pA_values >= pB_values[i]], var_values, lty = i)
}
legend("topleft", legend = paste("pB =", pB_values), lty = 1:length(pB_values))

# Cramer-Rao Lower Bound
lower_bound <- function(pA, pB){
  qA <- 1 - pA
  qB <- 1 - pB
  lb <- (qA*qB*(pA+pB))/(qA + qB)^3
  return(lb)
}
pB_values <- c(0.1, 0.3, 0.49)
pA_values <- seq(0, 0.1, by = 0.01)
low_values <- sapply(pA_values[pA_values >= pB_values[1]], function(x) lower_bound(x, pB_values[1]))
plot(pA_values[pA_values >= pB_values[1]], low_values, type = "l", lty = 1, ylim = range(0, 0.4),  
     xlab = "pA", ylab = "Var(NA(n))")

for (i in 2:length(pB_values)) {
  low_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) lower_bound(x, pB_values[i]))
  lines(pA_values[pA_values >= pB_values[i]], low_values, lty = i)
}
legend("topleft", legend = paste("pB =", pB_values), lty = 1:length(pB_values))


# Asymptotic variance vs lower bound
pB_values <- c(0.1, 0.3, 0.48)
pA_values <- seq(0, 0.975, by = 0.01)
par(mfrow = c(1,3))
for (i in 1:length(pB_values)) {
  var_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) asym_variance(x, pB_values[i]))
  lb_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) lower_bound(x, pB_values[i]))
  
  plot(pA_values[pA_values >= pB_values[i]], var_values, type = "l", lty = 1, ylim = range(0, var_values),
       xlab = "pA", ylab = "Var(NA(n))", cex.lab = 1.25, cex.axis = 1.25)
  lines(pA_values[pA_values >= pB_values[i]], lb_values, lty = 2)
  
  legend("topleft", legend = c("Asymptotic Variance", "Lower Bound"), lty = 1:2, cex = 1.2)
  title(main = paste("(", i, ")"), cex = 1.25)
}

par(mfrow=c(1,1))
pB_values <- c(0.1, 0.3, 0.48)
pA_values <- seq(0, 0.975, by = 0.01)
for (i in 1:length(pB_values)) {
  var_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) asym_variance(x, pB_values[i]))
  lb_values <- sapply(pA_values[pA_values >= pB_values[i]], function(x) lower_bound(x, pB_values[i]))
  
  plot(pA_values[pA_values >= pB_values[i]], var_values, type = "l", lty = 1, ylim = range(0, var_values),
       xlab = "pA", ylab = "Var(NA(n))", cex.lab = 1.25, cex.axis = 1.25)
  lines(pA_values[pA_values >= pB_values[i]], lb_values, lty = 2)
  
  legend("topleft", legend = c("Asymptotic Variance", "Lower Bound"), lty = 1:2, cex = 1.2)
}

