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
nSim <- 2500

# Sampling using LHS and rescaling the data

sampleLHS <- randomLHS(n = nSim, k = 2)
sampleLHS[,1] <- 0.5 + (0.9 - 0.5) * sampleLHS[,1]
sampleLHS[,2] <- 0.2 + (0.2 - 0.02) * sampleLHS[,2]
priorLHS <- data.frame(beta = sampleLHS[,1], sigma = sampleLHS[,2])

# The LHS

targetLHS <- NULL
for (i in 1:nSim) {
  pred <- RunModel(beta = priorLHS$beta[i], sigma = priorLHS$sigma[i],
                   data = surv.track)
  predEstim <- pred[1]$model.prev
  for (j in 1:11) {
      targetI <- ifelse(predEstim[j] >= surv.track$prev.LCI[j] &
                        predEstim[j] <= surv.track$prev.UCI[j], 1, 0)
  targetLHS <- c(targetLHS, targetI)
  }
}
sum(targetLHS)
toSelect <- which(targetLHS == 1)
summary(priorLHS$beta[toSelect])
summary(priorLHS$sigma[toSelect])
