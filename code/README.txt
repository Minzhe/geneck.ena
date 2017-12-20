-----------------------
analysis
-----------------------
Contains the code to compare different network construction methods.
	job - folder to generate and store job script to submit to the server
	combine.result.R - combine result from different simulation settings
	compare.net.roc.fun.R - function to compare network construction methods
	compare.ROC.R - wrap up the comparison function to R command line tool
	roc.fun.R - utility function to calculate auc.

-----------------------
lib
-----------------------
Library of functions needed for different algorithms.
	CMI2NI.R - core funciton of cmi2ni method
	ENA.R - core function of ena method
	glasso_sf.R - core function of glasso_sf method
	PCAC_CMI.R - core funciton of pcacmi method

-----------------------
replicate
-----------------------
Code to simulate network and compare performance of different method in replicate and plot boxplot.
	auc.replicate.R - pipeline to simulate network and compare performance in different simulation settings

-----------------------
simulation
-----------------------
Code to simulate network.
	ggm.R - core funciton to simulate ggm network.
	self_defined_network.Rdata - contains the network structure to simulate for ggm
	simulate.one.round.R - wrapped function to simulate network
	simulation.fun.R - core function to simulate network default method
	simulation.ggm.fun.R - wrapped function to simulate network with ggm