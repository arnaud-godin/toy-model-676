###############################
### Run of the model
### Arnaud Godin
### Spring 2018
### McGill University
###############################

require(ggplot2)
require(deSolve)
source('SIRS_model.R')

# The run of the model is a function by itself to then be used in the
# calibration process

# ---- Processing survey data ----

# Preparing the different data points
year <- c(2003:2013)
surv.prevalence <- c(.643, .699, .707, .703, .735, .750, .759, .753, .707, 
	.696, .720)
surv.N <- c(389, 581, 461, 566, 532, 492, 464, 489, 478, 461, 486)
surv.I <- round(surv.prevalence * surv.N)
surv.S <- surv.N - surv.I

# The data frame
surv.track <- data.frame(year = year, N = surv.N, I = surv.I, S = surv.S,
 prevalence = surv.prevalence)

#The confidence interval for prevalence (binomial)
surv.track$prev.LCI <- sapply(surv.track, function(x,y) binom.test(x = surv.track$I, y = surv.tra))
	surv.track$prev.UCI[i] <- binom.test(surv.track$I[i],
	 surv.track$N[i])$confint[2]
}

# ---- Run of the model ----

RunModel <- function(beta=2, alpha.A=0.15, alpha.Tx=.9, gamma=0.75/100,
 sigma = 0.024, dt=.01, data=surv.track, plot = ''){
	# --- Initialize the model parameters ---
	# There are only acute infections at first

	N <- 10000	# Total number of people
	A0 <- (N*0.643)/2	# Acutely infected (Prevalent cases/2)
	C0 <- (N*0.643)/2	# Chronic cases (Prevalent cases/2)
	Tx0 <- 0
	F0 <- 0	
	S0 <- N-A0-C0 # The susceptibles are total population - acute

	#	Sigma definition in the model

	# The initial population and parameters
	init.pop <- c('S' = S0, 'A' = A0, 'C' = C0, 'Tx' = Tx0, 'F' = F0)
	params <- c(beta = beta, alpha.A = alpha.A,
	 alpha.Tx = alpha.Tx, gamma = gamma, sigma = sigma)

	# Time definition for the model
	time.0 <- 2003
	time.t <- 2018
	vec.time <- seq(0, time.t - time.0, by = dt)

	# --- Run the model ---

	out <- as.data.frame(ode(init.pop, vec.time, SIRSModel, params))

	# Set time to a reasonable scale
	out$time <- time.0 + vec.time
	out$I <- out$A + out$C + out$F
	out$N <- out$S + out$A + out$C + out$Tx + out$F
	out$prev <- out$I / out$N

	# Estimated prevalence in the years of the surveys 
	out.prev <- subset(out$prev, out$time %in% c(2003:2013))
	prev.estim <- data.frame(c(2003:2013), out.prev)

	#Returning graphical representation
	base.graph <- ggplot(data = out, aes(x = time)) +
	theme(panel.grid.major.x = element_blank(),
        title = element_text(family = 'Times', size = 12),
        axis.title.x = element_text(family = 'Times', size = 10),
        axis.title.y = element_text(family = 'Times', size = 10),
        legend.title = element_text(family = 'Times', size = 10)) +
	theme_bw()

	epid.trajectory <- base.graph +
	geom_line(aes(y = S, color = 'Susceptibles')) +
	geom_line(aes(y = I, color = 'Infectious')) +
	geom_line(aes(y = Tx, color = 'On treatment')) +
	labs(x = 'Time (years)', y = 'Number of people')

	epid.prevalence <- base.graph +
	geom_line(aes(y = prev, color = 'Prevalence')) +
	labs(x = 'Time (years)', y = 'Prevalence')

	#Plotting either the trajectory or the prevalence
	if(plot == ''){
		return(prev.estim)
	} else if (plot == 'trajectory'){
		return(list(prev.estim = prev.estim, plot.traject = epid.trajectory))
	} else if (plot == 'prevalence'){
		return(list(prev.estim = prev.estim, plot.prev = epid.prevalence))
	} else {
		return(print('Options are: "", "trajectory", "prevalence"'))
	}
}

