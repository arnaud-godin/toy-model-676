##########################
### Force of infection    
### Arnaud Godin       
### Spring 2018        
### McGill University  
##########################

FoI <- function (beta, ILow, NLow, IHigh, NHigh) {
  # Function for the force of infection
  # For now, it will only allow for assortative mixing, but it will be extended
  # to allow for other types of mixing this is

  # --- Mixing matrix for semi-assortative mixing --- 

  # Assortative mixing (Todo -> command from func.parm to choose mixing type)
  mixLL <- 0.88     # The mixing probability in low risk
  mixHH <- mixLL    # Assumed to be the same as in the low risk 

  Mx <- matrix(data = c(mixLL, 1 - mixLL, mixHH, 1 - mixHH), nrow = 2, ncol = 2)
	
  # --- The forces of infection themselves

  # -- Low risk group --

  lambdaLow <- beta * ((Mx[1,1] * ILow/NLow) + (Mx[1,2] * IHigh/NHigh))

  # -- High risk group --

  lambdaHigh <- 1.5 * beta * ((Mx[2,2] * IHigh/NHigh) + (Mx[2,1] * ILow/NLow))

  # Returning the forces of infection
	return(list(FOILow = lambdaLow, FOIHigh = lambdaHigh))
}