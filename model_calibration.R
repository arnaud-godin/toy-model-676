###############################
### Model calibration
### Arnaud Godin
### Spring 2018
### McGill University
###############################

source('run_model.R')
library(lhs)

# ---- Using brute force and latin hypercube sampling ----

# Number of simulations
nSim <- 500

# Sampling using LHS and rescaling the data

sampleLHS <- randomLHS(n = nSim, k = 2)
scaledBeta <- 0.5 + (0.9 - 0.5) * sampleLHS[,1]
scaledSigma <- 0.02 + (0.2 - 0.02) * sampleLHS[,2]
priorLHS <- data.frame(beta = scaledBeta, sigma = scaledSigma)

# The LHS

choice <- function (predEstim, LCI, UCI) {
  ifelse(predEstim >= LCI & 
         predEstim <= UCI, 1, 0)
      }

targetLHS <- NULL
for (i in 1:nSim) {
  pred <- RunModel(beta = priorLHS$beta[i], sigma = priorLHS$sigma[i],
                   data = surv.track)
  predEstim <- pred$prev.estim$model.prev
  #Fit on the target
  target <- sum(mapply(choice, predEstim, surv.track$prev.LCI,
                            surv.track$prev.UCI))
  #Check all points fit withing the bounds
  targetI <- ifelse(target == 11, 1, 0) 
  targetLHS <- c(targetLHS, targetI)
}

sum(targetLHS)
toSelect <- which(targetLHS == 1)
summary(priorLHS$beta[toSelect])
summary(priorLHS$sigma[toSelect])
