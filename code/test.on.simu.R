###                         test.on.simu.R                          ###
### =============================================================== ###
# This R script is test ena performance on simulation data.

source("code/ena.R")

sim.huge <- read.csv("data/huge.nSamp500.Sigma1.csv", check.names = FALSE)
est_edge.huge <- network.ena(expr.data = sim.huge, n.perm = 10, sig.quant = 0.99)
truth.huge <- read.csv("data/huge.csv")
est_edge.huge.roc <- pre.edge.roc(est_edge = est_edge.huge, truth_edge = truth.huge)





### ================  function =================== ###
clean.net <- function(net) {
    idx <- net$Source >= net$Target
    net[idx,] <- data.frame(net$Target[idx], net$Source[idx])
    idx <- net$Source == net$Target
    net <- net[!idx,]
    return(net)
}

pre.edge.roc <- function(est_edge, truth_edge) {
    est_edge.roc <- data.frame(est_edge, Truth = 0)
    for (i in 1:nrow(truth_edge)) {
        Source <- truth_edge[i,1]
        Target <- truth_edge[i,2]
        idx <- which(Source == est_edge.roc$Node1 & Target == est_edge.roc$Node2)
        if (sum(idx) != 0) {
            est_edge.roc$Truth[idx] <- 1
            next
        } else {
            idx <- which(Source == est_edge.roc$Node2 & Target == est_edge.roc$Node1)
            est_edge.roc$Truth[idx] <- 1
        }
        
    }
    return(est_edge.roc)
}