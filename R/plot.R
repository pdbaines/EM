
"updateplot" <- function(obj,main="",addline=TRUE)
{
	if (is.null(obj$paths)){
		stop("No paths found")
	}
	if (is.null(obj$iter)){
		stop("Not a valid EM object")
	}
	plot(y=obj$paths$theta[2:t1$iter],x=obj$paths$theta[1:(t1$iter-1)],xlab="theta^{(t)}",ylab="theta^{(t+1)}",main=main)
	if (addline){
		if (obj$iter>2){
			abline(coef(lm(obj$paths$theta[2:obj$iter]~obj$paths$theta[1:(obj$iter-1)])),col="red",lwd=1.5)
		}
	}
}

