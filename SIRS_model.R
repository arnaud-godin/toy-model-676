###############################
### HCV transmission model
### Arnaud Godin
### Spring 2018
### McGill University
###############################

require(stats)
source('force_of_infection.R', chdir = TRUE)

# ---- Model of HCV transmission ----

SISModel <- function (t, pop, parms) {
	# --- The compartments ---

	# --Low risk individuals--

	SLow <- pop[1]		# Susceptible
	ALow <- pop[2]		# Acute 
	CLow <- pop[3]		# Chronic
	TxLow <- pop[4]		# On treatment
	FLow <- pop[5]		# Failed treatment

	# Population size for low risk individuals
	NLow <- SLow + ALow + CLow + TxLow + FLow 

	# The infectious, which imply the assumtion that people undergoing
	# treatment are not infectious because their viral load is low 
	ILow <- ALow + CLow + FLow

	# --High risk individuals--

	SHigh <- pop[6]		#Same disposition as before	
	AHigh <- pop[7]		
	CHigh <- pop[8]		
	TxHigh <- pop[9]
	FHigh <- pop[10]

	# Population size for high risk individuals
	NHigh <- SHigh + AHigh + CHigh + TxHigh + FHigh

	# Infectious in the high risk population
	IHigh <- AHigh + CHigh + FHigh
	
	# --- Sigma parameter of treatment coverage ---

	# The data is the following (2006
	txRate <- data.frame(year = c(3, 6, 12, 18),
                      r = c(0.024, 0.047, 0.05, 0.055))
	estim.sigma <- splinefun(txRate$year, txRate$r, method = 'hyman')

	############## The ODE system with parameter list #####################
	### beta = Probability of effective contact
	### alpha.A = Spont. clearance rate; multiplied by 2 (1/D)
	### gamma = HCV mortality rate (assumed equal for C and F)	
	### alpha.Tx = HCV Treatment efficacy
	### There is a total of 8 parameters to estimate
	#######################################################################

	# --- The overall differential equations model ---

	with(as.list(parms), {
		# The sigma parameter
		sigma <- estim.sigma(t)

		# Computing the force of infection for both groups

		FOI <- FoI(beta, ILow, NLow, IHigh, NHigh)

		# The differential equation system for the low risk
		dS <- - FOI$low * SLow + alpha.A * 2 * ALow + (alpha.Tx) * 4 * TxLow -
      r * SLow + r * NLow
		dA <- FOI$low * SLow - (1 - alpha.A) * 2 * ALow - alpha.A * 2 * ALow -
		 r * ALow
		dC <- (1 - alpha.A) * 2 * ALow - sigma * CLow - gamma * CLow	- r * CLow 
		dTx <- sigma * CLow - (1 - alpha.Tx) * 4 * TxLow + sigma * FLow -
		 alpha.Tx * 4 * TxLow  
		dF <- (1 - alpha.Tx) * 4 * TxLow - sigma * FLow - gamma * FLow - r * FLow

		# The differential equation system for the high risk

		#Returning the results
		result <- c(dS, dA, dC, dTx, dF)
		list(result)
	})

}

