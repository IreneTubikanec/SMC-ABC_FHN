
#-----------------------------------------------------------------------------------
# Author: Irene Tubikanec
# Date:  2024-05-24
#-----------------------------------------------------------------------------------
# Required packages for SBP SMC-ABC inference in the stochastic FHN model
#-----------------------------------------------------------------------------------

#Rcpp
library(Rcpp)
library(RcppNumerical)
library(devtools)
find_rtools(T)

#Parallelization
library(foreach)
library(doSNOW)
library(doParallel)

#Simulation of a path of the stochastic FHN model (Cpp-code in the SplittingFHN-package)
library(SplittingFHN)

#Multivariate normal distribution
library(mvnfast) 

