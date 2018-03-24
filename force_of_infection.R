##########################
### Force of infection ###   
### Arnaud Godin       ###
### Spring 2018        ###
### McGill University  ###
##########################

FoI <- function(beta, I, N){
	lambda <- beta * (I/N)
	return(lambda)
}