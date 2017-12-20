#################################################################
###                    one round simulate                     ###
#################################################################

source("code/simulation/simulation.fun.R")
source("code/simulation/simulation.ggm.fun.R")

simulate.network <- function(base.folder, method, size, n.sample, noise, seed = 0) {
      if (method == "simulation") {
            network.simulation(base.folder = base.folder, size = size, sample = n.sample, noise = noise, seed = seed)
      } else if (method == "ggm") {
            simulate.ggm(base.folder = base.folder, size = size, n.sample = n.sample, sigma = noise, seed = seed)
      }
}
