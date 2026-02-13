
library(dplyr)

#### Inverse logit to calculate theta ####
inverse_logit <- function(mu, alpha, beta) {
  ans <- exp(mu + alpha + beta) / (1 + exp(mu + alpha + beta))
  return(ans)
}

#### Function to remove labels for unknown data ####
sample_known <- function(n, m, data){
  id <- sample(1:n, size = m * n)
  r_m <- rep(NA, n)
  r_m[id] <- 1
  y_mat <- cbind(data, r_m)
  y_mat$r_m <- y_mat$r_m * y_mat$r
  y_mat <- y_mat %>% relocate(r, .after = last_col())
  y_mat <- y_mat[order(y_mat$r_m, decreasing = TRUE),]
  return(y_mat)
}

#### Function to turn some labels incorrect ####

sample_wrong <- function(n, m, s, data){
  m_known <- sum(!is.na(data[,ncol(data)-1]))
  id <- sample(1:m_known, size = s * m_known)
  
  for (i in id) {
    current <- data[i, ncol(data)-1]
    
    if (current == 1) {
      data[i, ncol(data)-1] <- sample(c(2, 3), 1)
    } else if (current == 2) {
      data[i, ncol(data)-1] <- sample(c(1, 3), 1)
    } else if (current == 3) {
      data[i, ncol(data)-1] <- sample(c(1, 2), 1)
    }
  }
  
  return(data)
}

#### Function to obtain data for one of the n ####
get_data <- function(n, seed, theta_mat){
  set.seed(seed)
  
  true_prop_known <- matrix(nrow = 6, ncol = 3)
  
  y_mat <- data.frame(matrix(NA, nrow = n, p))
  y_mat$r<- rep(1:R, length.out = n)
  
  for(i in 1:n){
    r <- y_mat[i, ncol(y_mat)]
    for(j in 1:p){
      y_mat[i,j] <- rbinom(1,1,prob = theta_mat[r,j])
    }
  }
  
  # shuffle the rows 
  # indices <- sample(1:n, n)
  #  y_mat <- y_mat[indices,]
  
  y_m10 <- sample_known(n, 0.1, y_mat)
  y_m30 <- sample_known(n, 0.3, y_mat)
  
  y_m10s10 <- sample_wrong(n, 0.1, 0.1, y_m10)
  y_m10s30 <- sample_wrong(n, 0.1, 0.3, y_m10)
  y_m10s50 <- sample_wrong(n, 0.1, 0.5, y_m10)
  
  y_m30s10 <- sample_wrong(n, 0.3, 0.1, y_m30)
  y_m30s30 <- sample_wrong(n, 0.3, 0.3, y_m30)
  y_m30s50 <- sample_wrong(n, 0.3, 0.5, y_m30)
  
  all_mat <- list(m10s10 = y_m10s10,
                  m10s30 = y_m10s30,
                  m10s50 = y_m10s50,
                  m30s10 = y_m30s10,
                  m30s30 = y_m30s30,
                  m30s50 = y_m30s50)
  
  iter = 1
  for(s in 1:6){
    true_prop_known[iter, ] <- round(as.numeric(table(all_mat[[s]]$r_m))/
                                            sum(!is.na(all_mat[[s]]$r_m)), 4)
    iter = iter + 1
  }
  
  return(append(all_mat, list(true_prop_known)))
}

#### Outward function to get repetitions of datasets ####

repeat_data <- function(iter, n) {
  total_list <- vector(mode = "list", length = iter)
  
  # for(i in 1:iter){
  #   total_list[[i]] <- get_data(n = n, seed = i, theta_mat)
  # }
  
  total_list <- lapply(1:iter, function(i) {
  get_data(n = n, seed = i, theta_mat)
})

  return(total_list)
}

#### Compute theta matrix given a vector ####
theta_func <- function(par_vec, R = 3, p = 5, q = 2){
  # q is just to start after mu
  # R-2 bc R-1 parameters in the vector, and last one is 0 - sum of these
  # same idea with p-2
  mu <- par_vec[1]
  alpha_vec <- c(par_vec[q:(q+R-2)], 0-sum(par_vec[q:(q+R-2)]))
  beta_vec <- c(par_vec[(q+R-1):(q+R-1 + p-2)], 
                0-sum(par_vec[(q+R-1):(q+R-1 + p-2)]))
  
  theta_mat <- outer(alpha_vec, beta_vec, function(a, b) inverse_logit(mu, a, b))
  
  #   theta_mat <- matrix(nrow = R, ncol = p)

  # for(r in 1:R){
  #   for(j in 1:p){
  #     theta_mat[r,j] <- inverse_logit(mu, alpha_vec[r], beta_vec[j])
  #   }
  # }
  return(theta_mat)
}

#### calculate complete data ll part 2

# sum_i sum_r sum_j zir(yij log thetarj + (1-yij) log (1-thetarj))
complete_data_ll_p2 <- function(theta_mat, z_mat, y_data, R){
  n <- nrow(y_data)
  p <- ncol(y_data) - 2
  ll <- 0 
  for(i in 1:n){
    for(j in 1:p){
      for(r in R){
        ll = ll + z_mat[i,r] * (y_data[i,j] * log(theta_mat[r,j]) + 
                                  (1 - y_data[i,j]) * log(1 - theta_mat[r,j])) 
      }
    }
  }
  return(ll)
}

optim_complete_data_ll <- function(par_vec, z_mat, R, p, y_data, theta_func){
  
  current_theta_mat <- theta_func(par_vec, R, p)
  current_theta_mat[current_theta_mat <= 0] <- 0.00001
  ans <- complete_data_ll_p2(current_theta_mat, z_mat, y_data, R)
  return(ans)
}



