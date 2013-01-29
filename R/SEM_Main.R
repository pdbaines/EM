
## SEM Code

"DM.estimate" <- function(theta.mle,theta.t,y.obs,update,fixed,verbose=FALSE)
{
	# NOTE: We assume that "update" is a function formatted in the 
	# same form as required for the EM() function i.e.,
	# update <- function(theta.t,y.obs,fixed,verbose){...}

	# Generic code to estimate the DM matrix from an EM-update function.
	p <- length(theta.mle)
	if (length(theta.t) != p){
		stop("length of 'theta.mle' and 'theta.t' must match")
	}
	
	# Return matrix of EM-update map derivatives:
	R <- matrix(NA,nrow=p,ncol=p)

	for (i in 1:p){

		# Update the user on current progress:
		if (verbose){
			cat(ppaste("Computing DM for element ",i,"...\n"))
		}

		# Define theta.t.i:
		theta.t.i <- theta.mle
		theta.t.i[i] <- theta.t[i]

		# Run one-iteration of EM from \theta^{t}(i):
		theta.t1.i <- update(theta.t.i,y.obs,fixed,verbose)

		# Compute R (see eq. 3.3.1 in Meng & Rubin, 1991):
		R[i,] <- (theta.t1.i - theta.mle)/(theta.t[i] - theta.mle[i])

		if (verbose){
			cat("======================================\n")
			cat(ppaste("Computing row ",i,"/",p," of DM matrix...\n"))
			cat(ppaste("theta^{(t)}(",i,"):\n"))
			print(theta.t.i)
			cat(ppaste("theta^{(t+1)}(",i,"):\n"))
			print(theta.t1.i)
			cat("======================================\n")
			cat(ppaste("Finished element ",i,"...\n"))
		}

	}

	# Return all stuff...
	return(R)
}


