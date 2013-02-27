My first edit... :)

EM: An R Package
================

Author: Paul D. Baines
Original Date: 01/16/2011
Version Date: 01/30/2013

General tools for the EM Algorithm, including MCEM convergence monitoring.

To install:

R CMD INSTALL EM

To use:

(1) Write an EM-update function of the form:

  my.em.update.func <- function(theta.t,y.obs,fixed,verbose){
    # Do cool stuff...
    theta.t1 <- ... # Find update
    return(theta.t1)
  }

(2) Use the "EM" function to perform the EM algorithm with your update function,
    this includes all built-in convergence monitoring, options to save sample
    paths, store log-likehood paths (if the log-likelihood function is supplied), 
    and other useful stuff.

(2') Use the "MCEM" function if your EM algorithm is a Monte Carlo EM algorithm.
     The MCEM function includes improved parameter estimation and convergence
     monitoring tools specifically tailored to MCEM algorithms.

Enjoy! :)

Git basics:

Local files ===> (git add) ===> Index (git commit) ===> Local repository ===> (git push) ===> Remote repository


Git setup:

git config --global user.name "Paul Baines"
git config --global user.email "pdbaines@ucdavis.edu"

Got to new directory:

git init
git add . 
git commit -m "First comment for project"

Note that the git add includes all subdirectories.

To push to Github, use Windows or Mac client, and select "Push to Github".

Or, 

git remote add origin https://github.com/pdbaines/LogNlogS.git
git push -u origin master

Do local editing, then to commit changes and send to remote repo:

git commit -m "Update after changes
git push

or: git push origin
or: git push origin master

git push 

To update local files to those on Github:

git pull

git pull https://github.com/pdbaines/LogNLogS.git
git pull origin master

git fetch https://github.com/pdbaines/LogNLogS.git
git diff

git fetch origin

Check status:

git status



