###############################
### Likelihood functions
### Arnaud Godin
### Spring 2018
### McGill University
###############################

binom_ll <- function (N, Obs, Mod) {
  #Computes the log-likelihood of X
  # N is the sample size
  # Obs is the observed prevalence in the sample
  # Mod is the predicted prevalence according to the model

  LL <- N * Obs * log(Mod) + N * (1 - Obs) * log(1 - Mod)
  return(LL)
}