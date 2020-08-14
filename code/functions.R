# set of support functions

simulate_covid <- function(allparvals, model, nsim=500){
  sims <- pomp::simulate(model, 
                         params=allparvals, 
                         nsim=nsim, format="data.frame", 
                         include.data=FALSE)
  
}


summarize_contact_tracing <- function(alpha, kappa, beta0, beta1, gamma, sigma, q, theta, Ntot, inivals, nsim=500){
  # function to return the mean of a a large number of contact tracing simulations
  
  parvals <- c(log_beta_0 = log(beta0/Ntot), 
                  log_beta_1 = log(beta1/Ntot), 
                  log_alpha = log(alpha*beta0/Ntot),
                  log_kappa = log(kappa),
                  log_gamma = log(gamma),
                  log_sigma = log(sigma),
                  q = q, 
                  log_theta_cases = log(theta)
  )
  
  sims <- simulate_covid(c(parvals,inivals), model=contact_tracing_model, nsim=nsim)
  
  outbreak_size <- sims %>%
    group_by(.id) %>%
    summarize(size = max(R))
  
  return(mean=mean(outbreak_size$size))
  
}

vCT <- Vectorize(summarize_contact_tracing, c("alpha", "kappa"))   #vectorized version of function to return mean of contact tracing


summarize_quarantine <- function(alpha, kappa, beta0, beta1, gamma, sigma, q, theta, Ntot, inivals, nsim=500){
  # function to return the mean of a a large number of contact tracing simulations
  
  parvals <- c(log_beta_0 = log(beta0/Ntot), 
               log_beta_1 = log(beta1/Ntot), 
               log_alpha = log(alpha*beta0/Ntot),
               log_kappa = log(kappa),
               log_gamma = log(gamma),
               log_sigma = log(sigma),
               q = q, 
               log_theta_cases = log(theta)
  )
  
  sims <- simulate_covid(c(parvals,inivals), model=quarantine_model, nsim=nsim)
  
  outbreak_size <- sims %>%
    group_by(.id) %>%
    summarize(size = max(R))
  
  return(mean=mean(outbreak_size$size))
  
}

vQ <- Vectorize(summarize_quarantine, c("alpha", "kappa"))   #vectorized version of function to return mean of quarantine

summarize_certification <- function(alpha, kappa, beta0, beta1, gamma, sigma, delta, theta, Ntot, inivals, nsim=500){
  # function to return the mean of a a large number of contact tracing simulations
  
  parvals <- c(log_beta_0 = log(beta0/Ntot), 
               log_beta_1 = log(beta1/Ntot), 
               log_alpha = log(alpha),
               log_kappa = log(kappa),
               log_gamma = log(gamma),
               log_sigma = log(sigma),
               log_delta = log(delta),
               log_theta_cases = log(theta)
  )
  
  sims <- simulate_covid(c(parvals,inivals), model=certification_model, nsim=nsim)
  
  outbreak_size <- sims %>%
    group_by(.id) %>%
    mutate(R = Rc+Ru) %>%
    summarize(size = max(R))
  
  return(mean=mean(outbreak_size$size))
  
}

vCert <- Vectorize(summarize_certification, c("alpha", "kappa"))   #vectorized version of function to return mean of contact tracing