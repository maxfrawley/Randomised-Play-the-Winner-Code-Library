# Normal RPW(alpha, beta) Simulation
RPW <- function(n, alpha, beta, pA, pB){
  nA <- numeric(n+1) # Allocate space for delayed response
  nB <- numeric(n+1)
  treat <- numeric(n)
  resp <- numeric(n)
  
  nA[1] <- nB[1] <- alpha 
  
  successes_A <- 0
  successes_B <- 0 
  
  for (i in 1:n) {
    # Treatment Allocation and Response
    treat[i] <- sample(c("A", "B"), 1, prob = c(nA[i], nB[i]))  
    if(treat[i] == "A"){
      resp[i] <- sample(c(1,0), 1, prob = c(pA, 1 - pA))
    }
    else{
      resp[i] <- sample(c(1,0), 1, prob = c(pB, 1 - pB))
    }
    # Updating delayed
    if(treat[i] == "A"){
      if(resp[i] == 1){
        nA[i+1] <- nA[i] + beta
        nB[i+1] <- nB[i]
        successes_A <- successes_A + 1
      } else{
        nB[i+1] <- nB[i] + beta
        nA[i+1] <- nA[i]
      }
    }  
    else{
      if(resp[i] == 1){
        nB[i+1] <- nB[i] + beta
        nA[i+1] <- nA[i]
        successes_B <- successes_B + 1
      } else{
        nA[i+1] <- nA[i] + beta
        nB[i+1] <- nB[i]
      }
    }
  }
  NAn <- sum(treat == "A")
  NBn <- sum(treat == "B")
  pA_hat <- if(NAn > 0){successes_A/NAn} else{0}
  pB_hat <- if(NBn > 0){successes_B/NBn} else{0}
  return(cbind(NAn, pA_hat, NBn, pB_hat))
}
# Average Power Function
average_power <- function(k, n, pA, pB){
  reject <- 0
  type_I <- numeric(k)
  for (i in 1:k){
    RPW_result <- RPW(n, 1, 1, pA, pB)
    nA <- RPW_result[1]
    nB <- RPW_result[3]
    pA_hat <- RPW_result[2] 
    pB_hat <- RPW_result[4]
    
    # Calculate test stat
    RPW_test_stat <- (pA_hat - pB_hat)/(sqrt((pA_hat*(1-pA_hat))/nA + (pB_hat*(1-pB_hat))/nB))
    # Check test stat is not NaN
    if (is.nan(RPW_test_stat)) {
      RPW_test_stat <- 0
    }
    z_critical <- qnorm(1 - 0.05)
    
    # Test Hypothesis
    if (RPW_test_stat > z_critical) {
      reject <- reject + 1
    }
  }
  average_power <- reject/k
  return(average_power)
}

# Power against sample size for fixed pB - Plot
pA_values <- c(0.2, 0.3, 0.4, 0.5)
N_values <- seq(10, 500, by = 10)
power_matrix <- matrix(0, nrow = length(N_values), ncol = length(pA_values))
for (i in 1:length(pA_values)) {
  for (j in 1:length(N_values)) {
    power_matrix[j, i] <- average_power(500, N_values[j], pA_values[i], 0.1)
  }
}
plot(N_values, power_matrix[, 1], pch = 15, col = 1 , xlab = "Sample Size", ylab = "Power", ylim = range(0, power_matrix), cex = 0.75)
for (i in 2:length(pA_values)) {
  points(N_values, power_matrix[, i], pch = 15 + i, col = i, cex = 0.75)
}
legend("bottomright", legend = paste("pA =", pA_values), pch = 15:(15+length(pA_values)), col = 1:length(pA_values))

# Required sample size for RPW
required_sample_size <- function(target_power, pA, pB, k, max_n) {
  current_n <- 10  # starting sample size
  while (current_n <= max_n) {
    avg_power <- average_power(k, current_n, pA, pB)
    if (avg_power >= target_power) {
      return(current_n)
    }
    current_n <- current_n + 1  # increase sample size by 10 in each iteration
  }
  return(NA)  # return NA if no sample size meets the target power within the specified range
}

target_power <- 0.8
pBval <- 0.7
pAvals <- c(0.9)
required_sample_sizes <- numeric(length(pAvals))
for (i in 1:length(pAvals)) {
  required_sample_sizes[i] <- required_sample_size(target_power, pAvals[i], pBval, 1000, 500)
}
data.frame(pA = pAvals, required_sample_size = required_sample_sizes)


# Required Sample Size for Equal Allocation
sample_size_equal <- function(pA, pB, gamma = 0.05, beta = 0.2) {
  z_alpha = qnorm(1-gamma)
  z_beta = qnorm(1 - beta)
  n = (2 * (z_alpha + z_beta)^2 * (pA*(1-pA) + pB*(1-pB))) / (pA - pB)^2
  return(ceiling(n))
}

pBval2 <- 0.1
pAvals2 <- c(0.3)
required_sample_sizes_equal <- numeric(length(pAvals2))
for (i in 1:length(pAvals2)) {
  required_sample_sizes_equal[i] <- sample_size_equal(pAvals2[i], pBval2)
}
data.frame(pA = pAvals2, required_sample_size = required_sample_sizes_equal)


# Power Calculations for Fluoxetine Trial 
# Delayed RPW(alpha, beta) Simulation
RPW_delayed_deterministic <- function(d, n, alpha, beta, pA, pB){
  nA <- numeric(n+d+1) # Allocate space for delayed response
  nB <- numeric(n+d+1)
  
  # For loop for nA[i] and nB[i] up until nA[d] and nB[d]
  for(i in 1:(d+1)){
    nA[i] <- alpha
    nB[i] <- alpha
  }
  treat <- numeric(n)
  resp <- numeric(n)
  
  successes_A <- 0
  successes_B <- 0 
  
  # First d Treatments and Responses
  for(i in 1:d){
    treat[i] <- sample(c("A", "B"), 1, prob = c(nA[i], nB[i]))
    if(treat[i] == "A"){
      resp[i] <- sample(c(1,0), 1, prob = c(pA, 1 - pA))
    }               
    else{
      resp[i] <- sample(c(1,0), 1, prob = c(pB, 1 - pB))
    }
  }
  # Treatments and Responses d+1 to n
  for (i in (d+1):n) {
    # Updating delayed
    if(treat[i-d] == "A"){
      if(resp[i-d] == 1){
        nA[i+1] <- nA[i] + beta
        nB[i+1] <- nB[i]
        successes_A <- successes_A + 1
      } else{
        nB[i+1] <- nB[i] + beta
        nA[i+1] <- nA[i]
      }
    }  
    else{
      if(resp[i-d] == 1){
        nB[i+1] <- nB[i] + beta
        nA[i+1] <- nA[i]
        successes_B <- successes_B + 1
      } else{
        nA[i+1] <- nA[i] + beta
        nB[i+1] <- nB[i]
      }
    }
    # Treatment Allocation and Response
    treat[i] <- sample(c("A", "B"), 1, prob = c(nA[i], nB[i]))  
    if(treat[i] == "A"){
      resp[i] <- sample(c(1,0), 1, prob = c(pA, 1 - pA))
    }
    else{
      resp[i] <- sample(c(1,0), 1, prob = c(pB, 1 - pB))
    }
  }
  NAn <- sum(treat == "A")
  NBn <- sum(treat == "B")
  pA_hat <- if(NAn > 0){successes_A/NAn} else{0}
  pB_hat <- if(NBn > 0){successes_B/NBn} else{0}
  return(cbind(NAn, pA_hat, NBn, pB_hat))
}

# Average Power Function for RPW with Delayed Response
average_power_delayed <- function(k, n, pA, pB){
  reject <- 0
  type_I <- numeric(k)
  for (i in 1:k){
    RPW_result <- RPW_delayed_deterministic(6, n, 1, 1, pA, pB)
    nA <- RPW_result[1]
    nB <- RPW_result[3]
    pA_hat <- RPW_result[2] 
    pB_hat <- RPW_result[4]
    
    # Calculate test stat
    RPW_test_stat <- (pA_hat - pB_hat)/(sqrt((pA_hat*(1-pA_hat))/nA + (pB_hat*(1-pB_hat))/nB))
    # Check test stat is not NaN
    if (is.nan(RPW_test_stat)) {
      RPW_test_stat <- 0
    }
    z_critical <- qnorm(1 - 0.05)
    
    # Test Hypothesis
    if (RPW_test_stat > z_critical) {
      reject <- reject + 1
    }
  }
  average_power <- reject/k
  return(average_power)
}

# Empirical Power Calulcations for Fluoxetine Trial
pA_values <- c(0.6, 0.7)
pB_values <- c(0.2, 0.3)
average_power_matrix <- matrix(0, nrow = length(pA_values)*length(pB_values), ncol = 3)
average_power_matrix[,1] <- rep(pA_values, each = length(pB_values))
average_power_matrix[,2] <- rep(pB_values, length(pA_values))

for (i in 1:length(pA_values)) {
  for (j in 1:length(pB_values)) {
    index <- (i - 1) * length(pB_values) + j
    average_power_matrix[index, 3] <- average_power_delayed(10000, 60, pA_values[i], pB_values[j])
  }
} 
average_power_matrix

# Required Sample Size for RPW with Delayed Responses
required_sample_size_delayed <- function(target_power, pA, pB, k, max_n) {
  current_n <- 10  # starting sample size
  while (current_n <= max_n) {
    avg_power <- average_power_delayed(k, current_n, pA, pB)
    if (avg_power >= target_power) {
      return(current_n)
    }
    current_n <- current_n + 1 # increase sample size by 10 in each iteration
  }
  return(NA)  # return NA if no sample size meets the target power within the specified range
}

# Required Sample Sizes for Fluoxetine Trial
target_power <- 0.8
pA_values <- c(0.6, 0.7, 0.8)
pB_values <- c(0.2, 0.3)
required_sample_sizes <- matrix(0, nrow = length(pA_values)*length(pB_values), ncol = 3)
required_sample_sizes[,1] <- rep(pA_values, each = length(pB_values))
required_sample_sizes[,2] <- rep(pB_values, length(pA_values))

for (i in 1:length(pA_values)) {
  for (j in 1:length(pB_values)) {
    index <- (i - 1) * length(pB_values) + j
    required_sample_sizes[index, 3] <- required_sample_size_delayed(target_power, pA_values[i], pB_values[j], 1000, 500)
  }
} 
required_sample_sizes


# Required Sample Sizes for ECMO Trial
target_power <- 0.8
pBval <- 0.2
pAvals <- c(0.55, 0.6, 0.65, 0.70)
required_sample_sizes <- numeric(length(pAvals))
for (i in 1:length(pAvals)) {
  required_sample_sizes[i] <- required_sample_size_delayed(target_power, pAvals[i], pBval, 1000, 500)
}
data.frame(pA = pAvals, required_sample_size = required_sample_sizes)
