
#MTH210 PROJECT
#BY ISHA GOEL


# data input
obs_data <- c(0.907,0.562,0.363,2.016,1.066,2.425,1.035,2.434,0.659,
              1.730,1.225,0.754,1.362,1.276,0.312,1.601,0.911,2.481,
              2.087,2.846,0.993,1.369,1.958,0.908,1.069)


# function to estimate alpha given lambda
alpha_func <- function(data_vals, lam){
  n_val <- length(data_vals)
  a_val <- -n_val / sum(log(1 - exp(-lam * data_vals)))
  return(a_val)
}

# newton raphson for lambda 
nr_solver <- function(data_vals){
  
  tol <- 0.001
  diff <- 100
  max_iter <- 100
  iter <- 0
  
  n_val <- length(data_vals)
  
  # initial guess using median
  lam_curr <- log(2) / median(data_vals)
  a_curr   <- alpha_func(data_vals, lam_curr)
  
  while(diff > tol && iter < max_iter){
    
    # first derivative
     g1 <- n_val/lam_curr - sum(data_vals) +
      (a_curr - 1) * sum((data_vals * exp(-lam_curr * data_vals)) /
                           (1 - exp(-lam_curr * data_vals)))
    
    # second derivative
     g2 <- -n_val/(lam_curr^2) -
      (a_curr - 1) * sum((data_vals^2 * exp(-lam_curr * data_vals)) /
                           ((1 - exp(-lam_curr * data_vals))^2))
    
    # update step
    lam_new <- lam_curr - g1 / g2
    
    # update alpha
    a_new <- alpha_func(data_vals, lam_new)
    
    # convergence check
     diff <- abs(lam_new - lam_curr) + abs(a_new - a_curr)
    
    iter <- iter + 1
    lam_curr <- lam_new
    a_curr   <- a_new
  }
  
  return(list(lambda = lam_curr,
              alpha  = a_curr,
              iter   = iter))
}



# mle estimation
res_main <- nr_solver(obs_data)

lam_hat <- res_main$lambda
a_hat   <- res_main$alpha
it_used <- res_main$iter

cat("Lambda MLE:", lam_hat)
cat("Alpha MLE :", a_hat)
cat("Iterations:", it_used)


# non parametric bootstrap
n_val <- length(obs_data)
B <- 1000

set.seed(123)

lam_np <- numeric(B)
a_np   <- numeric(B)


for(i in 1:B){
  
  # resampling with replacement
  samp <- sample(obs_data, size = n_val, replace = TRUE)
  
  res_np <- nr_solver(samp)
  
  lam_np[i] <- res_np$lambda
  a_np[i]   <- res_np$alpha
}

# remove invalid values
lam_np_clean <- lam_np[is.finite(lam_np)]
a_np_clean   <- a_np[is.finite(a_np)]

# ci using quantiles
lam_ci_np <- quantile(lam_np_clean, c(0.025, 0.975))
a_ci_np   <- quantile(a_np_clean,  c(0.025, 0.975))

cat("Non-Parametric Bootstrap 95% CI")
cat("Lambda CI: [", lam_ci_np[1], ",", lam_ci_np[2], "]")
cat("Alpha CI: [", a_ci_np[1],  ",", a_ci_np[2],  "]")


# parametric bootstrap
set.seed(123)

lam_p <- numeric(B)
a_p   <- numeric(B)

# generator function
gen_sample <- function(n_val, lam, a_val){
  u <- runif(n_val)
  x <- (-1/lam) * log(1 - u^(1/a_val))
  return(x)
}

for(i in 1:B){
  
  # generate from fitted model
  samp_p <- gen_sample(n_val, lam_hat, a_hat)
  
  res_p <- nr_solver(samp_p)
  
  lam_p[i] <- res_p$lambda
  a_p[i]   <- res_p$alpha
}


# remove invalid values
lam_p_clean <- lam_p[is.finite(lam_p)]
a_p_clean   <- a_p[is.finite(a_p)]

# ci
lam_ci_p <- quantile(lam_p_clean, c(0.025, 0.975))
a_ci_p   <- quantile(a_p_clean,  c(0.025, 0.975))


#results
cat("ALL RESULTS")
cat("Parametric Bootstrap 95% CI")
cat("Lambda CI: [", lam_ci_p[1], ",", lam_ci_p[2], "]\n")
cat("Alpha  CI: [", a_ci_p[1],  ",", a_ci_p[2],  "]\n")

cat("Non-Parametric Bootstrap 95% CI")
cat("Lambda CI: [", lam_ci_np[1], ",", lam_ci_np[2], "]")
cat("Alpha CI: [", a_ci_np[1],  ",", a_ci_np[2],  "]")


cat("Lambda MLE:", lam_hat)
cat("Alpha MLE :", a_hat)
cat("Iterations:", it_used)

