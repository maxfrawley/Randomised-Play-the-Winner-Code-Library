# Expected Expected Value of N_A(n)/n Function
E_NA <- function(n, alpha, beta, p_A, p_B) {
  q_B <- 1 - p_B
  E_NA <- numeric(n)
  for (i in 1:n) {
    prod_term <- 1
    for (k in (i+1):n) {
      prod_term <- prod_term*(1 + (beta*(p_A - q_B))/(2*alpha + beta*(k - 1)))
    }
    E_NA[i] <- ((alpha + beta*(i - 1)*q_B) / (2*alpha + beta*(i - 1)))  *prod_term
  }
  return(sum(E_NA)/n)
}

# Tables for Exact Expected Value of N_A(n)/n
p_values <- seq(0.1, 0.9, by = 0.01)
result_table <- matrix(NA, nrow = length(p_values), ncol = length(p_values))
colnames(result_table) <- p_values
rownames(result_table) <- p_values
for (i in 1:length(p_values)) {
  for (j in 1:length(p_values)) {
    # Calculate expected value
    result_table[i, j] <- E_NA(n = 50, alpha = 1, beta = 1, p_A = p_values[i], p_B = p_values[j])
  }
}
result_df <- as.data.frame(result_table)
print(result_df)

result_table2 <- matrix(NA, nrow = length(p_values), ncol = length(p_values))
colnames(result_table2) <- p_values
rownames(result_table2) <- p_values
for (i in 1:length(p_values)) {
  for (j in 1:length(p_values)) {
    # Calculate expected value
    result_table2[i, j] <- E_NA(n = 50, alpha = 5, beta = 1, p_A = p_values[i], p_B = p_values[j])
  }
}
result_df2 <- as.data.frame(result_table2)
print(result_df2)

# 3D surface plots for 
library(plotly)
plot_ly(z = as.matrix(result_df), x = colnames(result_df), y = rownames(result_df), type = "surface",
        colorscale = "Blues") %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(text = "pB", font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      yaxis = list(
        title = list(text = "pA", font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      zaxis = list(
        title = list(text = "Expected Allocation Proportion", font = list(size = 18)),
        tickfont = list(size = 14)
      )
    )
  )

plot_ly(z = as.matrix(result_df2), x = colnames(result_df2), y = rownames(result_df2), type = "surface",
        colorscale = "Blues") %>%
  layout(
    scene = list(
      xaxis = list(
        title = list(text = "pB", font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      yaxis = list(
        title = list(text = "pA", font = list(size = 18)),
        tickfont = list(size = 14)
      ),
      zaxis = list(
        title = list(text = "Expected Allocation Proportion", font = list(size = 18)),
        tickfont = list(size = 14)
      )
    )
  )
