# cure-estimator

Simulation code for the manuscript 
"Estimating Causal Effects from a Randomized Clinical Trial when Noncompliance is Measured with Error".
Designed for R 3.2.0

Each R file also has documentation in the file.

sim.R - main simulation code. Written to run simulation on 5 servers
using 12 ores.

true_mu.R - finds the value of E{Y*(1, 0)} = mu(1, 0) for given simulation parameters. 

functions.R - functions used by all the R files.

example_data.txt - example data set for use in example_analysis.R

example_analysis.R - example analysis of the data in example_data.txt
