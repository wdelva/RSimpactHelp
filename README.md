# The RSimpactHelper package

<!-- Created by Wim Delva, 21 September 2018 -->


This repository contains R functions that facilitate post-simulation analysis in R of simulation output for models developed with SimpactCyan (<https://simpactcyan.readthedocs.io/en/latest/>). 


## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system

  We have only tested this code on personal computers (OS X Version 10.11.6 and Linux Ubuntu Version 16.04), and on the lengau cluster at the Cape Town Centre for High Performance Computing (CHPC) and the golett cluster of the Flemish Supercomputer Centre (VSC).

### Required software

  R version 3.4.4

  Seq-Gen version 1.3.4. <https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4> Simulates viral evolution across a transmission network.

  FastTree version 2.1.10. <http://www.microbesonline.org/fasttree/#Install> Reconstructs a phylogenetic tree from a large alignment dataset.

  To install and load the package, run the following lines of code:
    
  library(devtools)

  install_github("wdelva/RSimpactHelp”, dependencies = TRUE)
  
  library(RSimpactHelper)

   

## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0 License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this work is suitable for commercial use, that licensees can modify the work, that they must release the source alongside with Derivative Work, and that Derivative Work must be released under the same terms.


## CONTACT INFORMATION

Please contact Prof Wim Delva with questions regarding the files provided and their authorized usage:

Wim Delva
Email: <DELVAW@sun.ac.za>

## INFORMATION FOR RSIMPACTHELP DEVELOPERS

### Adding package dependencies

If you want to use a function from a package, it not only needs to be listed in the NAMESPACE file, but you must also list the package as a dependency in the DESCRIPTION file. 

If the package is on CRAN, it should be listed under "Imports" in the DESCRIPTION file. 

If the package is hosted on another repository or website, it gets listed under "Imports", as well as "Remotes" in the DESCRIPTION file. How you list it under "Remotes" is dependent upon where you have to install it from. Here are the directions for listing it under Remotes: https://cran.r-project.org/web/packages/devtools/vignettes/dependencies.html . 
If the packages is from Bioconductor you need to, first list "biocViews" before "Imports", and then include the name of package in "Imports". 
