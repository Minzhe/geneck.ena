#################################################################
###                       simulate.all.R                      ###
#################################################################


###############  simulation  ###############
source("code/simulation/simulation.fun.R")
sizes <- c("tiny", "small", "moderate", "middle", "large", "huge")
n.samples <- c(50, 100, 200, 500)
sigmas <- c(0.2, 0.5, 1)

set.seed(1234)

for (size in sizes) {
      for (n.sample in n.samples) {
            for (sigma in sigmas) {
                  network.simulation(base.folder = "simulation", size = size, sample = n.sample, noise = sigma)
            }
      }
}
                  


##############  simulation ggm  ##################
source("code/simulation/simulation.ggm.fun.R")
sizes <- c("tiny", "small", "moderate", "middle", "large", "huge")
n.samples <- c(50, 100, 200, 500)
sigma2s <- c(0.2, 0.5, 1)
for (size in sizes) {
      for (n.sample in n.samples) {
            for (sigma2 in sigma2s) {
                  simulate.ggm(base.folder = "ggm", size = size, n.sample = n.sample, sigma2 = sigma2)
            }
      }
}
