###############################
### HCV transmission model
### Arnaud Godin
### Spring 2018
### McGill University
###############################

require(stats)
source('force_of_infection.R', chdir = TRUE)

# ---- MODIFICATION ----
a <- rnorm(1,1,1)

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

	# Total population size
	N <- NLow + NHigh

	# The total infected
	I <- ILow + IHigh
	
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

		# -- Differential equations system for the low risk --

		dSLow <- - FOI$low * SLow + alpha.A * 2 * ALow + alpha.Tx * 4 * TxLow -
     r * SLow + r * N * (NLow/N)
  
		dALow <- FOI$low * SLow - (1 - alpha.A) * 2 * ALow - r * ALow -
		 alpha.A * 2 * ALow
 
		dCLow <- (1 - alpha.A) * 2 * ALow - sigma * CLow - gamma * CLow	- r * CLow

		dTxLow <- sigma * CLow - (1 - alpha.Tx) * 4 * TxLow + sigma * FLow -
		 alpha.Tx * 4 * TxLow

		dFLow <- (1 - alpha.Tx) * 4 * TxLow - sigma * FLow - gamma * FLow -
		 r * FLow

		# -- Differential equations system for the high risk --

		dSHigh <- - FOI$high * SHigh + alpha.A * 2 * AHigh + alpha.Tx * 4 * TxHigh -
		 r * SHigh + r * N * (1 - (NLow/N))

		dAHigh <- FOI$high * SHigh - (1 - alpha.A) * 2 * AHigh - r * AHigh -
		 alpha.A * 2 * AHigh

		dCHigh <- (1 - alpha.A) * 2 * AHigh - sigma * CHigh - gamma * CHigh -
		 r * CHigh

		dTxHigh <- sigma * CHigh - (1 - alpha.Tx) * 4 * TxHigh + sigma * FHigh -
		 alpha.Tx * 4 * TxHigh

		dFHigh <- (1 - alpha.Tx) * 4 * TxHigh - sigma * FHigh - gamma * FHigh -
		 r * FHigh

		#Returning the results
		result <- c(dS, dA, dC, dTx, dF)
		list(result)
	})

}

