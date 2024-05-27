
Author: Irene Tubikanec
Date:   2024-05-24

Description: This code provides an implementation of the structure-based and preserving (SBP) SMC-ABC method for parameter estimation in the stochastic FitzHugh-Nagumo (FHN) model, proposed in the paper:
         
Inference for the stochastic FitzHugh-Nagumo model from real action potential data via approximate Bayesian computation,
by Adeline Samson, Massimiliano Tamborrino and Irene Tubikanec.

In particular, it reproduces the estimation results of the setting T=200, Delta_obs=0.02 when simulated data is observed (and can be easily adapted to other settings and real data experiments).

How does the code work?

1. Install relevant packages (see the R-file required_packages)

2. Install the package SplittingFHN (splitting-simulation of paths of the stochastic FHN model using C++ Code)

3. Run the R-file main_SMC_ABC_FHN (SBP SMC-ABC inference for the stochastic FHN model)
After termination of the algorithm, the kept ABC posterior samples and corresponding weights are stored into the folder ABC_Results.

4. Run the R-file plot_results (visualization of estimation results)
It will create a figure visualizing the marginal ABC posterior densities. The figure also reports the underlying true parameters values. 

Code description:
For a detailed description of the code, please consider the respective files

Licence information:
Please consider the txt-files LICENCE and COPYING.GPL.v3.0




