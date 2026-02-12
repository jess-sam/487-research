
#### Model Fitting ####

model_fit <- function(data_list, n, p = 5, R = 3, theta_func, pi_start, param_start, ll){
  # number of reps = 100
  #reps <- length(data_list)
  #num_scen <- 6
  reps <- 2
  num_scen <- 2
  
  mu_matrix <- matrix(nrow = reps, ncol = num_scen)
  alpha_matrix_R1 <- matrix(nrow = reps, ncol = num_scen)
  alpha_matrix_R2 <- matrix(nrow = reps, ncol = num_scen)

  beta_matrix_p1 <- matrix(nrow = reps, ncol = num_scen)
  beta_matrix_p2 <- matrix(nrow = reps, ncol = num_scen)
  beta_matrix_p3 <- matrix(nrow = reps, ncol = num_scen)
  beta_matrix_p4 <- matrix(nrow = reps, ncol = num_scen)

  pi1_matrix <- matrix(nrow = reps, ncol = num_scen)
  pi2_matrix <- matrix(nrow = reps, ncol = num_scen)
  pi3_matrix <- matrix(nrow = reps, ncol = num_scen)
  
  acc_matrix <- matrix(nrow = reps, ncol = num_scen)
  time_matrix <- matrix(nrow = reps, ncol = num_scen)

  for(rep in 1:2){
    for(d in 2:3){
      
      start_time <- Sys.time()
      
      # Set up the group membership matrix, z_mat
      y_data <- data_list[[rep]][[d]]
      
      label_swap_id <- matrix(nrow = 1, ncol = R)
      
      n_unknown <- sum(is.na(y_data$r_m))
      n_known <- n - n_unknown
      y_data_unknown <- y_data[is.na(y_data$r_m), -c((ncol(y_data)-1):ncol(y_data))]
      
      q = 1
      z_known <- matrix(NA, n_known, R)
      z_mat <- matrix(NA, n_unknown, R)
      
      for(g in y_data$r_m[!is.na(y_data$r_m)]) {
        z_known[q,g] <- 1
        z_known[q,-g] <- 0
        q = q + 1
      }
      
      # initialise params
      
      pi_vect <- pi_start
      param_vec <- param_start
      theta_mat <- theta_func(param_vec, R, p)
      
      these_pars <- c(param_vec, pi_vect)
      previous_pars <- rep(0, (R+p-1 + R))
      
      # Expectation and maximisation steps
      nloops <- 0
      
      while (max(abs(these_pars-previous_pars)) > 0.01) {
        # E step, use current pars to estimate z:
        numerators <- matrix(NA, n_unknown, R)
        for (i in 1:n_unknown){
          for (r in 1:R) {
            numerators[i,r] <- (
              pi_vect[r]*prod(theta_mat[r,]^y_data_unknown[i,]*
                                (1-theta_mat[r,])^(1-y_data_unknown[i,]))
            )
          }
        }
        
        z_mat <- numerators/apply(numerators,1,sum)
        z_mat <- rbind(z_known, z_mat)
        
        # M step, update parameter estimates by maximisation:
        # Update pi vector:
        pi_vect <- apply(z_mat,2,mean)
        
        # maximise model params
        optim_res <- optim(par = param_vec,
                           fn = optim_complete_data_ll,
                           z_mat = z_mat,
                           R = 3,
                           p = 5,
                           y_data = y_data,
                           theta_func = theta_func,
                           method = "L-BFGS-B",
                           hessian = F,
                           control = list(maxit = 10000, fnscale = -1))
        
        param_vec <- optim_res$par
        
        # update params
        previous_pars <- these_pars
        these_pars <- c(param_vec, pi_vect)
      }
      end_time <- Sys.time()
      time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
      time_matrix[rep, d-1] <- time_taken
      
      
      labels_est <- rep(NA, n_unknown)
      for(i in 1:n_unknown){
        labels_est[i] <- which(z_mat[(i+n_known),] == max(z_mat[i+n_known,]))
      }
      
      acc_matrix[rep, d-1] <- mean(labels_est == y_data$r[-c(1:n_known)])
      #print(labels_est)
      #print(y_data$r[-c(1:n_known)])
        #sum(labels_est == y_data$r[-c(1:n_known)])/n_unknown
      
      mu_matrix[rep, d-1] <- param_vec[1]
      alpha_matrix_R1[rep, d-1] <- param_vec[2]
      alpha_matrix_R2[rep, d-1] <- param_vec[3]
      beta_matrix_p1[rep, d-1] <- param_vec[4]
      beta_matrix_p2[rep, d-1] <- param_vec[5]
      beta_matrix_p3[rep, d-1] <- param_vec[6]
      beta_matrix_p4[rep, d-1] <- param_vec[7]

      pi1_matrix[rep, d-1] <- pi_vect[1]
      pi2_matrix[rep, d-1] <- pi_vect[2]
      pi3_matrix[rep, d-1] <- pi_vect[3]
      
    }
  }
    
  
  return(list(mu_matrix = mu_matrix,
              ar1 = alpha_matrix_R1,
              ar2 = alpha_matrix_R2,
              bp1 = beta_matrix_p1,
              bp2 = beta_matrix_p2,
              bp3 = beta_matrix_p3,
              bp4 = beta_matrix_p4,
              pi1 = pi1_matrix,
              pi2 = pi2_matrix,
              pi3 = pi3_matrix,
              acc = acc_matrix,
              time = time_matrix))
}
