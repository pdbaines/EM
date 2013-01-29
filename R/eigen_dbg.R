
########################################
## FUNCTIONS FOR EIGENVALUE DEBUGGING ##
########################################

"eigen.vals" <- function(x)
{
  ret <- try({eigen(x)$values},silent=TRUE)
  if (class(ret)=="try-error"){
    ret <- rep(NA,times=nrow(x))
  } else {
    # Take care of possible imaginary values and ordering:
    ret <- sort(Mod(ret),decreasing=TRUE) ## What about negative eigenvalues????
  }
  return(ret)
}

"clean.DM" <- function(x)
{
  return(x$DM.path[sapply(x$DM.path,function(y){class(y)!="try-error"})])
}

"eigen.plots" <- function(obj.list,alg.obj.names,DM.plots=TRUE,map.plots=TRUE,dm.mfrow=c(1,1),map.mfrow=c(1,1),print.every=Inf)
{
  ####
  ## TODO:
  ## Change input from alg.obj.names to list of objects...
  ####
  
  ####
  ## Notes:
  ## DM.plots requires compute.DM=TRUE for EM algorithms
  ## map.plots only requires keep.paths=TRUE
  ####
  
  ## Functionize this for all algorithms...
  par.names <- c("alpha","beta")

  ret <- vector("list",length(obj.list))
  for (rep.no in 1:length(obj.list)){
    ret[[rep.no]] <- list("DM"=NULL,"DM.eigen"=NULL)
    ret[[rep.no]]$DM <- vector("list",length=length(alg.obj.names))
    ret[[rep.no]]$DM.eigen <- vector("list",length=length(alg.obj.names))
    ret[[rep.no]]$slopes <- vector("list",length=length(alg.obj.names))
    names(ret[[rep.no]]$DM) <- names(alg.obj.names)
    names(ret[[rep.no]]$DM.eigen) <- names(alg.obj.names)
    names(ret[[rep.no]]$slopes) <- names(alg.obj.names)
  }
  
  for (rep.no in 1:length(obj.list)){
  
    for (i in 2:length(alg.obj.names)) # EXCLUDE optim()
    {
  
      tmp.alg <- names(alg.obj.names)[i]
      tmp.res <- ppaste("obj.list[[rep.no]]$'",tmp.alg,"'")
      eval(parse(text=ppaste("tmp.res.optim <- obj.list[[rep.no]]$'",names(alg.obj.names)[1],"'")))
      
      if (DM.plots){
	par(mfrow=dm.mfrow)
	## First: eigenvalue paths...
	eval(parse(text=paste("res.DM.sub <- clean.DM(",tmp.res,")",sep="")))
	if (!is.null(res.DM.sub)){
	  res.DM.eigen <- t(sapply(res.DM.sub,eigen.vals))
	  res.DM <- t(sapply(res.DM.sub,as.numeric))
	  matplot(res.DM.eigen,type="l",
		  main=paste(tmp.alg,": Max(ev) = ",round(median(res.DM.eigen[,1],na.rm=TRUE),4)),
		  xlab="iteration",ylab="eigenvalues(DM)")
	  # Store the DM map and eigenvalues:
	  ret[[rep.no]]$DM[[i]] <- res.DM
	  ret[[rep.no]]$DM.eigen[[i]] <- res.DM.eigen
	} else {
	  # Unable to compute eigenvalues...
	}
      } ## END if (DM.plots){...}

      if (map.plots){
	par(mfrow=map.mfrow)
	## Next: parameter paths...
	eval(parse(text=paste("tmp.nits <- ",tmp.res,"$iterations",sep="")))
	ret[[rep.no]]$slopes[[i]] <- rep(NA,length(par.names))
	names(ret[[rep.no]]$slopes[[i]]) <- par.names
	if (!is.null(tmp.nits) && !is.na(tmp.nits)){
	  for (j in 1:2){
	    ## Approximate the update map:
	    tmp.t  <- eval(parse(text=ppaste(tmp.res,"$paths$theta[1:(tmp.nits-1),j]")))
	    tmp.t1 <- eval(parse(text=ppaste(tmp.res,"$paths$theta[2:(tmp.nits),j]")))
	    sub.its <- c(1:length(tmp.t))[ceiling(length(tmp.t)/2):length(tmp.t)]
	    tmp.lm <- lm(tmp.t1[sub.its]~tmp.t[sub.its])
	    ret[[rep.no]]$slopes[[i]][j] <- coef(tmp.lm)[2]
	    eval(parse(text=paste("plot(x=tmp.t,y=tmp.t1,type='b',xlab='theta^{(t)}',ylab='theta^{(t+1)}',
				  main='",tmp.alg,": ",par.names[j],", slope=",round(coef(tmp.lm)[2],6),"')",sep="")))
	    if (!any(is.na(coef(tmp.lm)))){
	      abline(coef(tmp.lm),col="blue",lty=2,lwd=0.5)
	    } else {
	      print(tmp.lm)
	    }
	    points(x=tmp.res.optim$par[j],y=tmp.res.optim$par[j],pch=17,cex=1.8,col="red")
	    plot(y=tmp.t1[sub.its],x=tmp.t[sub.its],main=ppaste("Zoomed Subset: ",tmp.alg,": ",par.names[j],", slope=",round(coef(tmp.lm)[2],6)),
		xlab="theta^{(t)}",ylab="theta^{(t+1)}",type="b")
	    if (!any(is.na(coef(tmp.lm)))){
	      abline(coef(tmp.lm),col="blue",lty=2,lwd=0.5)
	    }
	  }
	} else {
	  # Not able to plot marginal parameter update map
	}
	
      } ## END if (map.plots){...}
    } ## END for (i in 2:length(alg.obj.names)){...} 
    
    if (is.finite(print.every) && (rep.no%%print.every==0)){
      cat(ppaste("Finished replicate number ",rep.no,"...\n"))
    }
	
  } ## END for (rep.no in 1:length(obj.list)){...}
  
  return(ret)
}

