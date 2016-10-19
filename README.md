# cure-estimator

Simulation code for the manuscript 
"Estimating Causal Effects from a Randomized Clinical Trial when Noncompliance is Measured with Error".
Designed for R 3.2.0

Each R file is documented internally.

sim.R - main simulation code. Written to run simulation on multiple cores
using the parallel package. 

true_mu.R - finds the value of E{Y*(1, 0)} = mu(1, 0) for given simulation parameters. 

functions.R - functions used in sim.R and true_mu.R

