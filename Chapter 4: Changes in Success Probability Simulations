# RPW(alpha_A, alpha_B, beta) Simulation - Initial Urn Composition of different values
RPW2 <- function(n, alpha_A, alpha_B, beta, pA, pB){
  nA <- numeric(n+1)
  nB <- numeric(n+1)
  treat <- numeric(n)
  resp <- numeric(n)
  
  nA[1] <- alpha_A
  nB[1] <- alpha_B
  
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
  return(list(
    'nA' = nA,
    'NAn' = NAn,
    'pA_hat' = successes_A/NAn,
    'nB' = nB,
    'NBn' = NBn,
    'pB_hat' = successes_B/NBn))
}

# Success change simulation
simulation_success_change <- function(n, n1, alpha, beta, pA1, pB1, pA2, pB2){
  # Simulation for first n1 trial
  RPW_1 <- RPW2(n1, alpha, alpha, beta, pA1, pB1)
  
  nA1 <- RPW_1$nA
  NAn1 <- RPW_1$NAn
  nB1 <- RPW_1$nB
  NBn1 <- RPW_1$NBn
  
  alpha_A <- nA1[n1]
  alpha_B <- nB1[n1]
  
  # Simulation for subsequent n2 = n-n1 trials
  RPW_2 <- RPW2(n-n1, alpha_A, alpha_B, beta, pA2, pB2)
  
  nA2 <- RPW_2$nA
  NAn2 <- RPW_2$NAn
  nB2 <- RPW_2$nB
  NBn2 <- RPW_2$NBn
  
  # Total results
  nA <- c(nA1, nA2)
  nB <- c(nB1, nB2)
  NAn <- NAn1 + NAn2
  NBn <- NBn1 + NBn2
  
  return(list(
    'nA' = nA,
    'NAn1' = NAn1,
    'NAn2' = NAn2,
    'NAn' = NAn,
    'NAn/n' = NAn/n,
    'nB' = nB,
    'NBn1' = NBn1,
    'NBn2' = NBn2, 
    'NBn' = NBn,
    'NBn/n' = NBn/n
    ))
}

# Mean allocation of success change 
changes_mean_allocation <- function(n_simulation, n, n1, alpha, beta, pA1, pB1, pA2, pB2){
  RPW_changes_allocation <- numeric(n_simulation)
  for (i in 1:n_simulation) {
    RPW_allocation <- E_NA(n, alpha, beta, pA1, pB1)
    RPW_changes_allocation[i] <- as.numeric(simulation_success_change(n, n1, alpha, beta, pA1, pB1, pA2, pB2)[5])
  }
  return(c(RPW_allocation, mean(RPW_changes_allocation)))
}

# Table of N_A(n)/n for Scenario 2
sim1 <- simulation_success_change(500, 125, 5, 1, 0.7, 0.2, 0.45, 0.45)
nA1 <- sim1$nA
nB1 <- sim1$nB

sim2 <- simulation_success_change(500, 250, 5, 1, 0.7, 0.2, 0.45, 0.45)
nA2 <- sim2$nA
nB2 <- sim2$nB

sim3 <- simulation_success_change(500, 375, 5, 1, 0.7, 0.2, 0.45, 0.45)
nA3 <- sim3$nA
nB3 <- sim3$nB

par(mfrow = c(1,3))

plot(1:length(nA1), nA1, type = "l", col = "blue", xlab = "Number of Treatments", ylab = "Number of Treatment Balls", main = "(a): n = 125",
     cex.lab=1.5, cex.axis=1.5, ylim = range(0, max(nA1, nB1)))
points(length(nA1), nA1[length(nA1)], pch = 16, col = "blue")
text(length(nA1), nA1[length(nA1)], nA1[length(nA1)], pos = 2, offset = 0.5, cex = 1.5)
points(nB1, type = "l", col = "green")
points(length(nB1), nB1[length(nB1)], pch = 16, col = "green")
text(length(nB1), nB1[length(nB1)], nB1[length(nB1)], pos = 2, offset = 0.5, cex = 1.5)
abline(v = 125, col = "red", lty = 2)
legend("topleft", legend = c("Treatment A", "Treatment B"),
       col = c("blue", "green"), lty = 1, cex = 1.2)

plot(1:length(nA2), nA2, type = "l", col = "blue", xlab = "Number of Treatments", ylab = "Number of Treatment Balls", main = "(b): n = 250",
     cex.lab=1.5, cex.axis=1.5, ylim = range(0, max(nA2, nB2)))
points(length(nA2), nA2[length(nA2)], pch = 16, col = "blue")
text(length(nA2), nA2[length(nA2)], nA2[length(nA2)], pos = 2, offset = 0.5, cex = 1.5)
points(nB2, type = "l", col = "green")
points(length(nB2), nB2[length(nB2)], pch = 16, col = "green")
text(length(nB2), nB2[length(nB2)], nB2[length(nB2)], pos = 2, offset = 0.5, cex = 1.5)
abline(v = 250, col = "red", lty = 2)
legend("topleft", legend = c("Treatment A", "Treatment B"),
       col = c("blue", "green"), lty = 1, cex = 1.2)

plot(1:length(nA3), nA3, type = "l", col = "blue", xlab = "Number of Treatments", ylab = "Number of Treatment Balls", main = "(c): n = 375",
     cex.lab=1.5, cex.axis=1.5, ylim = range(0, max(nA3, nB3)))
points(length(nA3), nA3[length(nA3)], pch = 16, col = "blue")
text(length(nA3), nA3[length(nA3)], nA3[length(nA3)], pos = 2, offset = 0.5, cex = 1.5)
points(nB3, type = "l", col = "green")
points(length(nB3), nB3[length(nB3)], pch = 16, col = "green")
text(length(nB3), nB3[length(nB3)], nB3[length(nB3)], pos = 2, offset = 0.5, cex = 1.5)
abline(v = 375, col = "red", lty = 2)
legend("topleft", legend = c("Treatment A", "Treatment B"),
       col = c("blue", "green"), lty = 1, cex = 1.2)
sim1 <- sim1[-1]
sim1 <- sim1[-5]
sim1 <- rbind(sim1)

sim2 <- sim2[-1]
sim2 <- sim2[-5]
sim2 <- rbind(sim2)

sim3 <- sim3[-1]
sim3 <- sim3[-5]
sim3 <- rbind(sim3)

(success_prob_change_results <- rbind(sim1, sim2, sim3))



# Treatment Proportion against Intervention Point Plot
par(mfrow=c(1,1))
n <- 100
N_values = seq(5, 100, by = 5)

success_change1 <- numeric(length(N_values))
success_change1[length(N_values)] <- E_NA(n, 1, 1, 0.7, 0.2)
success_change2 <- numeric(length(N_values))
success_change2[length(N_values)] <- E_NA(n, 1, 1, 0.7, 0.2)
success_change3 <- numeric(length(N_values))
success_change3[length(N_values)] <- E_NA(n, 1, 1, 0.7, 0.2)

for (i in 1:(length(N_values)-1)){
  success_change1[i] <- changes_mean_allocation(1000, n, N_values[i], 1, 1, 0.7, 0.2, 0.6, 0.3)[2]
  success_change2[i] <- changes_mean_allocation(1000, n, N_values[i], 1, 1, 0.7, 0.2, 0.45, 0.45)[2]
  success_change3[i] <- changes_mean_allocation(1000, n, N_values[i], 1, 1, 0.7, 0.2, 0.3, 0.6)[2]
}

library(ggplot2)
df1 <- data.frame(N_values, success_change1)
df2 <- data.frame(N_values, success_change2)
df3 <- data.frame(N_values, success_change3)

colnames(df1) <- c("N_values", "success_change")
colnames(df2) <- c("N_values", "success_change")
colnames(df3) <- c("N_values", "success_change")

df_combined <- rbind(
  transform(df1, Scenario = "Scenario 1"),
  transform(df2, Scenario = "Scenario 2"),
  transform(df3, Scenario = "Scenario 3")
)

# Plot the combined data
ggplot(df_combined, aes(x = N_values, y = success_change, color = Scenario)) +
  geom_smooth(size = 1, se = FALSE) +
  geom_point(size = 2) +
  labs(x = "Intervention Point",
       y = "Treatment Proportion to A") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red"))
