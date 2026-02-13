
# n = 3000 ----

source("functions.R")
source("model_fit.R")

# Parameter setup ----
set.seed(123)

p <- 5
R <- 3

alpha_vec <- c(-2, 0)
beta_vec <- c(-2, -1.5, 0.3, 1)
mu <- 0

# True theta matrix ----

theta_mat <- theta_func(c(mu, alpha_vec, beta_vec))

# Get datasets ----
all_n1000 <- repeat_data(iter = 100, n = 3000)

# Set starting values for parameters ----
pi_start <- rep(1/R, R)

set.seed(123)
mu_rand <- runif(1, min = mu - 1, max = mu + 1)
alpha_rand <- sapply(alpha_vec[-3], function (x){runif(1, min = x - 1, max = x + 1)})
beta_rand <- sapply(beta_vec[-5], function (x){runif(1, min = x - 1, max = x + 1)})

param_start <- c(mu_rand, alpha_rand, beta_rand)

# Fit Model ----
res3000 <- model_fit(all_n3000, 
                     n = 3000, 
                     theta_func = theta_func, 
                     pi_start = pi_start,
                     param_start = param_start,
                     ll = optim_complete_data_ll)

save(res3000, all_n3000, file = "n3000-res.Rdata")
