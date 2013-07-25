".onLoad" <-
    function (libname, pkgname)
{
    ## figure out this year automatically 
    this.year <- substr(as.character(Sys.Date( )), 1, 4)
    
    ## echo output to screen
    packageStartupMessage("##\n## EM")
    packageStartupMessage("## Copyright (C) 2013-", this.year,
                          " Paul D. Baines\n##\n", sep="")

    ## Load any compiled code:
    #library.dynam(pkgname, pkgname, lib.loc=libname)
}






