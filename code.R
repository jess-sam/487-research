#### Simulation experiments ####

source("functions.R")
source("model_fit.R")

#### Parameter setup #####

p <- 5
R <- 3

alpha_vec <- c(-2, 0, 2)
beta_vec <- c(-2, -1.5, 0.3, 1, 2.2)
mu <- 0

#### True theta matrix ####
theta_mat <- matrix(NA, nrow = R, ncol = p)

for(r in 1:R){
  for(j in 1:p){
    theta_mat[r,j] <- inverse_logit(mu, alpha_vec[r], beta_vec[j])
  }
}


#### Lists of 100 repetitions of all the n = 300, 1000, 3000 datasets)
all_n300 <- repeat_data(iter = 100, n = 300)
all_n1000 <- repeat_data(iter = 100, n = 1000)
all_n3000 <- repeat_data(iter = 100, n = 3000)



#################################################################


# Set starting values for parameters 

pi_start <- rep(1/R, R)

set.seed(123)
mu_rand <- runif(1, min = mu - 2, max = mu + 2)
alpha_rand <- sapply(alpha_vec[-3], function (x){runif(1, min = x - 2, max = x + 2)})
beta_rand <- sapply(beta_vec[-5], function (x){runif(1, min = x - 2, max = x + 2)})


param_start <- c(mu_rand, alpha_rand, beta_rand)


res <- model_fit(all_n3002, 
                 n = 300, 
                 theta_func = theta_func, 
                 pi_start = pi_start,
                 param_start = param_start,
                 ll = optim_complete_data_ll)

res#res$mu_matrix
#t(res$mu_matrix)

all_n3002[[1]][[7]]
res$
res

#####################


#### EM code ####

#test_data <- get_data(n = 300, seed = 0)
#test_dataset <- test_data$m10s10

#y_data_ex <- all_n300[[1]]$m10s10


#X_mat <- test_dataset[-ncol(test_dataset)]

#X_mat_unknown <- y_data_ex[is.na(y_data_ex$r_m), -ncol(y_data_ex)]

#X_mat_unknown

#n <- nrow(y_data_ex)
#n_unknown <- sum(is.na(y_data_ex$r_m))
#n_known <- n - n_unknown
############################

# create theta matrix with R = 3 rows

#mean_theta <- sum(X_mat)/(n*p) # Average prob. success
#theta_start <- matrix(NA,R,p)
# random variation

#for (r in 1:R) {
#  theta_start[r,] <- runif(p, mean_theta/2,(1+mean_theta)/2)
#}



all_n3002 <- repeat_data(iter = 10, n = 300)
table(all_n3002[[1]][[2]]$r_m) / sum(!is.na(all_n3002[[1]][[2]]$r_m))


length(is.na(all_n3002[[1]][[2]]$r_m))

all_n3002[[3]]
