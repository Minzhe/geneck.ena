###                         ENA.R                          ###
### ====================================================== ###
# This R script is function to use ensemble network aggregation method to constrcut gene network.


ena.rank <- function(net.list, method) {
    if (!(method %in% c("inverse.sum", "sum", "inverse.prod"))) {
        print("Method should be 'inverse.sum', 'sum' or 'inverse.prod'.")
        return()
    }
    
    n.net <- length(net.list)
    for (i in 1:n.net) {
        net.list[[i]][lower.tri(net.list[[i]], diag = TRUE)] <- NA
        net.list[[i]][upper.tri(net.list[[i]], diag = FALSE)] <- rank(-net.list[[i]][upper.tri(net.list[[i]], diag = FALSE)])
        if (method %in% c("inverse.sum", "inverse.prod")) {
            net.list[[i]] <- 1/net.list[[i]]
        }
    }
    if (method %in% c("inverse.sum", "sum")) {
        net.ena <- Reduce("+", net.list)
        net.ena <- net.ena / n.net
    } else {
        net.ena <- Reduce("*", net.list)
        net.ena < net.ena^(1/n.net)
    }
    if (method %in% c("inverse.sum", "inverse.prod")) net.ena <- 1/net.ena
    return(net.ena)
}


##############  function  ################
perm.net <- function(net.list, method) {
    for (i in 1:length(net.list)) {
        net.list[[i]][upper.tri(net.list[[i]], diag = FALSE)] <- sample(net.list[[i]][upper.tri(net.list[[i]], diag = FALSE)])
    }
    net.perm <- ena.rank(net.list, method)
    return(net.perm)
}

perm.net.n <- function(net.list, method, n.perm) {
    perm.v <- c()
    for (i in 1:n.perm) {
        net.perm <- perm.net(net.list, method)
        perm.v <- c(perm.v, net.perm[upper.tri(net.perm, diag = FALSE)])
    }
    return(perm.v)
}