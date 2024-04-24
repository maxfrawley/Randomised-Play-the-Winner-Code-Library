# RPW(alpha, beta) Simulation
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

# CI Functions for CIs for pA and pB
# Bootstrap Confidence Interval for pA and pB
bootstrap_ci1 <- function(observed_p, observed_N, B, gamma) {
  pA <- observed_p[1]
  pB <- observed_p[2]

  NAn <- observed_N[1]
  NBn <- observed_N[2]
  n <- NAn + NBn
  
  K <- length(observed_p)
  
  bootstrap_p_matrix <- matrix(NA, nrow = B, ncol = K)
  
  # Bootstrap procedure
  for (b in 1:B) {
    # Simulate treatment assignments and responses
    simulated_responses <- RPW(n, 1, 1, pA, pB) # Setting alpha = beta = 1
     
    # Compute bootstrap estimates
    bootstrap_p_matrix[b, 1] <- simulated_responses[2]
    bootstrap_p_matrix[b, 2] <- simulated_responses[4]
  }
  
  # Order bootstrap estimates for each treatment
  ordered_bootstrap_p <- t(apply(bootstrap_p_matrix, 2, sort))
  
  # Calculate bootstrap confidence intervals
  lower_bound <- ordered_bootstrap_p[, ceiling(B*gamma/2)]
  upper_bound <- ordered_bootstrap_p[, floor(B*(1 - gamma/2))]
  
  # Return confidence intervals
  return(cbind(observed_p, lower_bound, upper_bound))
}
# Coverage Probabilities and Average Lengths for Bootstrap CI 1
simulate_bootstrap_ci1 <- function(x, observed_p, observed_N, B, gamma) {
  coverage_probs <- matrix(NA, nrow = x, ncol = 2)
  lengths <- matrix(NA, nrow = x, ncol = 2)
  for (i in 1:x) {
    simulated_data <- RPW(observed_N[1] + observed_N[2], 1, 1, observed_p[1], observed_p[2])
    simulated_p <- c(simulated_data[2], simulated_data[4])
    simulated_N <- c(simulated_data[1], simulated_data[3])
    result <- bootstrap_ci1(simulated_p, simulated_N, B, gamma)
    lower_bound <- result[, 2]
    upper_bound <- result[, 3]
    coverage_probs[i,1] <- lower_bound[1] <= observed_p[1] & observed_p[1] <= upper_bound[1]
    coverage_probs[i,2] <- lower_bound[2] <= observed_p[2] & observed_p[2] <= upper_bound[2]
    lengths[i,] <- upper_bound - lower_bound
  }
  coverage_prob <- colMeans(coverage_probs)
  avg_lengths <- colMeans(lengths)
  return(cbind(observed_p, coverage_prob, avg_lengths))
}

# Large-sample Confidence Interval for pA and pB
large_sample_ci1 <- function(observed_p, observed_N, gamma) {
  pi <- observed_p
  Ni <- observed_N
  qi <- 1 - pi
  se <- sqrt((pi*qi)/Ni)
  z_value <- qnorm(1-gamma/2)
  lower_bound <- pi - z_value*se
  upper_bound <- pi + z_value*se
  return(cbind(observed_p, lower_bound, upper_bound))
}
# Coverage Probabilities and Average Lengths for Large-sample CI 1
simulate_large_sample_ci1 <- function(x, observed_p, observed_N, gamma) {
  coverage_probs <- matrix(NA, nrow = x, ncol = 2)
  lengths <- matrix(NA, nrow = x, ncol = 2)
  for (i in 1:x) {
    simulated_data <- RPW(observed_N[1] + observed_N[2], 1, 1, observed_p[1], observed_p[2])
    simulated_p <- c(simulated_data[2], simulated_data[4])
    simulated_N <- c(simulated_data[1], simulated_data[3])
    result <- large_sample_ci1(simulated_p, simulated_N, gamma)
    lower_bound <- result[, 2]
    upper_bound <- result[, 3]
    coverage_probs[i,1] <- (lower_bound[1] <= observed_p[1]) & (observed_p[1] <= upper_bound[1])
    coverage_probs[i,2] <- (lower_bound[2] <= observed_p[2]) & (observed_p[2] <= upper_bound[2])
    lengths[i,] <- upper_bound - lower_bound
  }
  coverage_prob <- colMeans(coverage_probs)
  avg_lengths <- colMeans(lengths)
  return(cbind(observed_p, coverage_prob, avg_lengths))
}


# CI Functions for CIs for delta = pA - pB
# Bootstrap Confidence Interval for delta = pA - pB
bootstrap_ci2 <- function(observed_p, observed_N, B, gamma) {
  pA_hat <- observed_p[1]
  pB_hat <- observed_p[2]
  delta_hat <- pA_hat - pB_hat
  nA <- observed_N[1]
  nB <- observed_N[2]
  n <- nA + nB

  bootstrap_p <- numeric(B)

  for (b in 1:B) {
    simulated_responses <- RPW(n, 1, 1, pA_hat, pB_hat)
    bootstrap_p[b] <- simulated_responses[2] - simulated_responses[4]
  }
  
  # Order bootstrap estimates for each treatment
  ordered_bootstrap_p <- sort(bootstrap_p)
  
  # Calculate bootstrap confidence intervals
  lower_bound <- ordered_bootstrap_p[ceiling(B*gamma/2)]
  upper_bound <- ordered_bootstrap_p[floor(B*(1 - gamma/2))]
  
  # Return confidence intervals
  return(cbind(delta_hat, lower_bound, upper_bound))
}
# Coverage Probabilities and Average Lengths for Bootstrap CI 2
simulate_bootstrap_ci2 <- function(x, observed_p, observed_N, B, gamma) {
  coverage_probs <- numeric(x)
  lengths <- numeric(x)
  delta_hat <- observed_p[1] - observed_p[2]
  for (i in 1:x) {
    simulated_data <- RPW(observed_N[1] + observed_N[2], 1, 1, observed_p[1], observed_p[2])
    simulated_p <- c(simulated_data[2], simulated_data[4])
    simulated_N <- c(simulated_data[1], simulated_data[3])
    result <- bootstrap_ci2(simulated_p, simulated_N, B, gamma)
    lower_bound <- result[2]
    upper_bound <- result[3]
    coverage_probs[i] <- (lower_bound <= delta_hat) & (delta_hat <= upper_bound)
    lengths[i] <- upper_bound - lower_bound
  }
  coverage_prob <- mean(coverage_probs)
  avg_lengths <- mean(lengths)
  return(cbind(delta_hat, coverage_prob, avg_lengths))
}

# Large-sample Confidence Interval for delta = pA - pB
large_sample_ci2 <- function(observed_p, observed_N, gamma) {
  pA_hat <- observed_p[1]
  pB_hat <- observed_p[2]
  delta_hat <- pA_hat - pB_hat
  nA <- observed_N[1]
  nB <- observed_N[2]
  
  se <- sqrt(pA_hat*(1-pA_hat)/nA + pB_hat*(1-pB_hat)/nB)
  z_value <- qnorm(1-gamma/2)
  
  lower_bound <- delta_hat - z_value*se
  upper_bound <- delta_hat + z_value*se
  return(cbind(delta_hat, lower_bound, upper_bound))
}
# Coverage Probabilities and Average Lengths for Large-sample CI 2
simulate_large_sample_ci2 <- function(x, observed_p, observed_N, gamma) {
  coverage_probs <- numeric(x)
  lengths <- numeric(x)
  delta_hat <- observed_p[1] - observed_p[2]
  for (i in 1:x) {
    simulated_data <- RPW(observed_N[1] + observed_N[2], 1, 1, observed_p[1], observed_p[2])
    simulated_p <- c(simulated_data[2], simulated_data[4])
    simulated_N <- c(simulated_data[1], simulated_data[3])
    result <- large_sample_ci2(simulated_p, simulated_N, gamma)
    lower_bound <- result[2]
    upper_bound <- result[3]
    coverage_probs[i] <- (lower_bound <= delta_hat) & (delta_hat <= upper_bound)
    lengths[i] <- upper_bound - lower_bound
  }
  coverage_probs <- na.omit(coverage_probs)
  lengths <- na.omit(lengths)
  
  coverage_prob <- mean(coverage_probs)
  avg_lengths <- mean(lengths)
  return(cbind(delta_hat, coverage_prob, avg_lengths))
}


# Singular simulated RPW trial
RPW_result <- RPW(100, 1, 1, pA, pB) # Input required values of pA and pB
observed_p <- c(RPW_result[2], RPW_result[4])
observed_N <- c(RPW_result[1], RPW_result[3])

# Calculated Values for pA and pB CIs
simulate_bootstrap_ci1(1000, observed_p, observed_N, 1000, 0.05)
simulate_large_sample_ci1(10000, observed_p, observed_N, 0.05)

# Calculated Values for delta = pA - pB CIs
simulate_bootstrap_ci2(1000, observed_p, observed_N, 1000, 0.05)
simulate_large_sample_ci2(10000, observed_p, observed_N, 0.05)


# Fluoxetine Trial
# Normal Stratum
# Observed Trial Outcomes
observed_p <-  c(0.619, 0.476)
observed_N <- c(21, 21)

# Bootstrap CIs for pA, pB and delta
bootstrap_ci1(observed_p, observed_N, 10000, 0.05)
bootstrap_ci2(observed_p, observed_N, 10000, 0.05)

# Large-sample CIs for pA, pB and delta
large_sample_ci1(observed_p, observed_N, 0.05)
large_sample_ci2(observed_p, observed_N, 0.05)

# Shortened Stratum
# Observed Trial Outcomes
observed_p <-  c(0.600, 0.333)
observed_N <- c(20, 21)

# Bootstrap CIs for pA, pB and delta
bootstrap_ci1(observed_p, observed_N, 10000, 0.05)
bootstrap_ci2(observed_p, observed_N, 10000, 0.05)

# Large-sample CIs for pA, pB and delta
large_sample_ci1(observed_p, observed_N, 0.05)
large_sample_ci2(observed_p, observed_N, 0.05)
