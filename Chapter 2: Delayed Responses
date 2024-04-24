# Normal RPW(alpha, beta) Simulation with no Delayed Response
RPW <- function(n, alpha, beta, pA, pB){
  nA <- numeric(n+1) # Allocate space for delayed response
  nB <- numeric(n+1)
  
  treat <- numeric(n)
  resp <- numeric(n)
  
  nA[1] <- nB[1] <- alpha # Set inital urn composition
  
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
  return(cbind(NAn, NAn/n, NBn, NBn/n))
}

# RPW for Deterministic Delayed Response for (d time units) Simulation
RPW_delayed_deterministic <- function(d, n, alpha, beta, pA, pB){
  nA <- numeric(n+d+1) # Allocate space for delayed response
  nB <- numeric(n+d+1)
  
  # For loop for nA[i] and nB[i] up until nA[d] and nB[d] - Intial urn composition
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
  return(cbind(NAn, NAn/n, NBn, NBn/n))
}

# Simulations for mean allocation to treatment A for delay up to 15 units
mean_allocation <- function(n_simulation, n, alpha, beta, pA, pB){
  RPW_allocation <- matrix(0, nrow = 16, ncol = n_simulation)
  rownames(RPW_allocation) <- c(paste(0:15))
  for (i in 1:n_simulation) {
    RPW_allocation[1, i] <- RPW(n, alpha, beta, pA, pB)[2]
    for(j in 2:16){
      RPW_allocation[j,i] <- RPW_delayed_deterministic(j, n, alpha, beta, pA, pB)[2]
    }
  }
  mean_values <- rowMeans(RPW_allocation)
  return(mean_values)
}

# Table for Delayed Response for pA and fixed pB
pA_values <- c(0.3, 0.5, 0.7, 0.9)
delayed_matrix <- matrix(0, nrow = length(pA_values), ncol = 16)
for (i in 1:length(pA_values)) {
    delayed_matrix[i, ] <- mean_allocation(5000, 50, 1, 1, pA_values[i], 0.1)
}
colnames(delayed_matrix) <- c(0:15)
rownames(delayed_matrix) <- pA_values
round(delayed_matrix[,c(1, 3, 6, 11, 16)], 3)

# Figure for Delayed Response
plot(0:15, delayed_matrix[1,], pch = 15, col = 1 , xlab = "Time Units of Delay", ylab = "Proportion of Patients Allocated to Treatment A", ylim = range(delayed_matrix), cex = 0.75)
for (i in 2:length(pA_values)) {
  points(0:15, delayed_matrix[i,], pch = 15 + i, col = i, cex = 0.75)
}
legend("topright", legend = paste("pA =", pA_values), pch = 15:(15+length(pA_values)), col = 1:length(pA_values))


# Table for Delayed Response in Fluoxetine Trial
mean_allocation <- function(n_simulation, n, alpha, beta, pA, pB){
  RPW_allocation <- numeric(n_simulation)
  RPW_delayed_allocation <- numeric(n_simulation)
  for (i in 1:n_simulation) {
    RPW_allocation[i] <- RPW(n, alpha, beta, pA, pB)[2]
    RPW_delayed_allocation[i] <- RPW_delayed_deterministic(6, n, alpha, beta, pA, pB)[2]
  }
  return(c(mean(RPW_allocation), mean(RPW_delayed_allocation)))
}

pA_values <- c(0.6, 0.7, 0.8)
pB_values <- c(0.1, 0.2, 0.3)
delayed_matrix <- matrix(0, nrow = length(pA_values), ncol = length(pB_values))
for (i in 1:length(pA_values)) {
  for (j in 1:length(pB_values)){
    delayed_matrix[i, j] <- mean_allocation(5000, 50, 1, 1, pA_values[i], pB_values[j])[2]
  }
}
colnames(delayed_matrix) <- pB_values
rownames(delayed_matrix) <- pA_values
round(delayed_matrix,3)
