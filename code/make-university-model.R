# make-university-model.R
#
# This script generates a pomp object for an SLAIR model of COVID 19 
# Running this script saves an RDS object: university-model.RDS
# The pomp object can be used for simulating trajectories and fitting

# # Clear the decks ---------------------------------------------------------

# rm(list = ls(all.names = TRUE))


# Load libraries ----------------------------------------------------------
# library(dplyr)
# library(pomp)
# library(here) #to simplify loading/saving into different folders

############################################################################
# Code to define  process model -------------------------------------------
############################################################################

#1 step function of model
pomp_step <- Csnippet(
    "
  // C indexes at 0, for R users making things 1 bigger and start 
  // with index 1, i.e. leave trans[0] empty, same for rate[0]
  
  double rate[7];
  double trans[7];

  double foi;  // total force of infection
  double f;  // force of infection (symptomatic)
  double g;  // force of infection (asymptomatic)

  // Forces of Infection (foi):
  foi = exp(log_beta)*(I + exp(log_bL)*L + exp(log_bA)*A);   // total force of infection
  f = (1-a)*foi;                                            // foi from S to L
  g = a*foi;                                                // foi from S to A
  
  // Compute the transition rates
  // Transitions from S
  rate[1] = f;                                        //infection, movement from S to L
  rate[2] = g;                                        //infection, movement from S to A

  // Transitions from L
  rate[3] = exp(log_sigma);                           //movement from L to I
  rate[4] = exp(log_xi);                              //movement from L to R

  // Transitions from I
  rate[5] = exp(log_gammaI);                          //recovery from I to R

  // Transitions from A
  rate[6] = exp(log_gammaA)+exp(log_xi);              //movement from A to R

  // Compute the state transitions
  reulermultinom(2, S, &rate[1], dt, &trans[1]); //move from S to L (latent infection) or A (asympotmatic infection)
  reulermultinom(2, L, &rate[3], dt, &trans[3]); //move from L to I (symtom onset) or R (removal due to testing)
  reulermultinom(1, I, &rate[5], dt, &trans[5]); //move from I to R
  reulermultinom(1, A, &rate[6], dt, &trans[6]); //move from A to R (natural recovery + removal due to testing) 

  // Apply transitions to state variables
  S += - trans[1] - trans[2];
  L += trans[1] - trans[3] - trans[4];
  I += trans[3] - trans[5];
  A += trans[2] - trans[6];
  R += trans[4] + trans[5] + trans[6];
  C_report += trans[5];
  C_surveil += trans[4] + trans[6];
  "
)

# C snippet for initial condition specification ---------------------------

rinit <- Csnippet(
  "
  S = nearbyint(S_0);
  L = nearbyint(L_0);
  I = nearbyint(I_0);
  A = nearbyint(A_0);
  R = nearbyint(R_0);
  C_report = 0;
  C_surveil = 0;
  "
)


############################################################################
# Code to define estimation components of model 
############################################################################

# Define likelihood function ----------------------------------------------

# Details on dnbinom_mu:
# dnbinom_mu is negative binomial parameterized by mu (see ?rnbinom)
# Given data (e.g. x = cases) and model-predicted value (e.g. size = C_new) and mean (e.g theta1), estimate probability.
# the value 1 in the last slot corresponds to log = 1.

#R code: dbinom(x=reports, size=H, prob=rho, log=log), C code: lik = dbinom(reports,H,rho,give_log)

dmeas <- Csnippet(
  "
  double d1; // log likelihood for C_symp
  double d2; // log likelihood for C_asymp
  double theta1; // dispersion parameter for C_symp
  double theta2; // dispersion parameter for C_asymp
  theta1 = exp(log_theta_cases_reported);
  theta2 = exp(log_theta_cases_surveilled);
  
  if(ISNA(reports)) {
    d1 = 0;  // loglik is 0 if no observations
  } else {
    d1 = dnbinom_mu(reports, theta1, C_report, 1); // negative binomial parameterized by mean/mu, which is 
  }
  
  if(ISNA(surveil)) {
    d2 = 0;  // loglik is 0 if no observations
  } else {
    d2 = dnbinom_mu(surveil, theta2, C_surveil, 1); // negative binomial parameterized by mean/mu, which is 
  }
  
  
  lik = d1 + d2;  // sum the individual likelihoods
  lik = (give_log) ? lik : exp(lik);  // return loglik or exp(lik)
  "
)


# Define process simulator for observations  ------------------------------
# given the model states (C_report, C_surveil)
# produce simulated data/obervations based on the specified distribution
# rnbinom_mu is negative binomial parameterized by mu (see ?rnbinom)
# Given model-predicted value (e.g. size = C_report) and mean (e.g theta1), 
# generate a single random value for expected observed cases (i.e. n=1).
# size/C_report is number of trials, the outcome is expected number of successes, given the specified mean
# Some examples from pomp website
#R code: reports = rbinom(n=1, size=H, prob=rho), C code: reports = rbinom(H,rho)

rmeas <- Csnippet(
  "
  double theta1;
  double theta2;
  theta1 = exp(log_theta_cases_surveilled);
  theta2 = exp(log_theta_cases_reported);
  reports = rnbinom_mu(theta1, C_report);  // for forecasting. 
  surveil = rnbinom_mu(theta2, C_surveil);  // for forecasting. 
  "
)


############################################################################
# Code to define variables, parameters and parameter transformations
############################################################################

# State variables to track ------------------------------------------------
varnames <- c("S", "L", "I", "A", "R", "C_report", "C_surveil")


# Parameters --------------------------------------------------------------
# Parameter and variable names
model_pars <- c("log_beta",   #basic transmissibility
                "log_bL",     # relative transmissibility of latent individuals
                "log_bA",     # relative transmissibility of asymptomatic individuals
                "a",          # fraction asymptomatic
                "log_xi",     # surveilance test rate
                "log_sigma",  # reciprocal of latent period
                "log_gammaI", # reciprocal of infectious period
                "log_gammaA"  # reciprocal of asymptomatic incubation period
)

measure_pars <- c("log_theta_cases_surveilled", "log_theta_cases_reported")

# Initial conditions of state variables are also parameters
ini_pars <- c("S_0", "L_0", "I_0", "A_0", "R_0")

#combine names of all parameters into a single vector
#note that values for parameters are specified in a separate script and can be loaded and assigned as needed
parnames <- c(model_pars,measure_pars,ini_pars)


#######################################################################
# Load cleaned data ---------------------------------------------------------
#######################################################################
pseudo_data <- data.frame(
  Date = seq.Date(from = as.Date("2020-08-20"), to = as.Date("2020-12-31"), by = "day"))
pseudo_data$time=seq(1,dim(pseudo_data)[1])

# # Data, if we get it, will be defined here
# filename = here('data/clean-CT-data.RDS') #data from covidtracking.com
# dat <- readRDS(filename)
# pomp_data <- dat %>%
#   dplyr::filter(state_full == "Georgia") %>%
#   dplyr::select(Date, cases) %>%
#   dplyr::arrange(Date) %>%
#   right_join(pseudo_data, by = "Date") %>%
#   dplyr::select(-hold) %>%
#   mutate(time = 1:n()) %>%
#   dplyr::select(time, cases)

# Define the pomp model object --------------------------------------------
university_model <- pomp(
  data = pseudo_data, 
  times = "time",
  t0 = 0,
  dmeasure = dmeas,
  rmeasure = rmeas,
  rinit = rinit,
  rprocess = euler(step.fun = pomp_step, delta.t = 1/20),
  statenames = varnames,
  paramnames = parnames, 
  obsnames = c("reports", "surveil"),
  accumvars = c("C_report", "C_surveil") 
)

# Save the pomp object ----------------------------------------------------
filename = here('output/university-model.RDS')
saveRDS(university_model, filename)
