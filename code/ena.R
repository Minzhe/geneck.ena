###                         ena.R                          ###
### ====================================================== ###
# This R script is function to use ena to constrcut gene network.

suppressMessages(library(corpcor))
suppressMessages(library(GeneNet))
suppressMessages(library(CDLasso))
suppressMessages(library(glasso))
source("code/lib/glasso_SF.R")
source("code/lib/PCA_CMI.R")
source("code/lib/CMI2NI.R")
suppressMessages(library(space))
source("code/lib/ENA.R")
library(reshape2)

network.ena <- function(expr.data, n.perm, sig.quant) {
    if (sig.quant <= 0 | sig.quant >= 1) {
        stop('Input error: parameter sig.quant for ena should be between 0 and 1.')
    }
    p <- ncol(expr.data)
    n <- nrow(expr.data)
    gene.index <- colnames(expr.data)
    expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
    
    est_edge <- list()
    ### GeneNet
    print("Constructing using GeneNet")
    print("=======================================")
    pcor_est <- pcor.shrink(expr.mat)
    test_pcor <- network.test.edges(pcor_est)
    est_edge.GeneNet <- matrix(0, p, p)
    for (i in 1:nrow(test_pcor)) {
        est_edge.GeneNet[test_pcor$node1[i], test_pcor$node2[i]] <- test_pcor$prob[i]
    }
    est_edge[[1]] <- est_edge.GeneNet
    file.remove("Rplots.pdf")
    
    ### ns
    print("Constructing using neighborhood selection")
    print("=======================================")
    est_res.ns <- matrix(0, p, p)
    alpha <- 0.2
    for (k in 1:p) {
        rsp <- expr.mat[, k]
        prd <- t(expr.mat[, -k])
        lam <- sqrt(sum(rsp ^ 2)) * qnorm(alpha / (2 * p ^ 2), lower.tail = F)
        out <- l2.reg(prd, rsp, lambda = lam)
        est_res.ns[k, -k] <- out$estimate
    }
    est_edge[[2]] <- (abs(est_res.ns) + t(abs(est_res.ns))) / 2
    
    ### glasso
    print("Constructing using GLASSO")
    print("=======================================")
    S <- t(expr.mat) %*% expr.mat / n
    out <- glasso(S, rho = 0.6)
    est_edge[[3]] <- abs(out$wi)
    
    ### glasso-sf
    print("Constructing using GLASSO-SF")
    print("=======================================")
    out <- glasso_sf(expr.mat, alpha = 0.3)
    est_edge[[4]] <- abs(out$wi)
    
    ### pcacmi
    print("Constructing using PCACMI")
    print("=======================================")
    out <- pca_cmi(t(expr.mat), 0.03)
    est_edge[[5]] <- abs(out$Gval)
    
    ### cmi2ni
    print("Constructing using CMI2NI")
    print("=======================================")
    out <- cmi2ni(t(expr.mat), 0.03)
    est_edge[[6]] <- abs(out$Gval)
    
    ### space
    print("Constructing using SPACE")
    print("=======================================")
    out <- space.joint(expr.mat, lam1 = 1 * n, iter = 5)
    est_edge[[7]] <- abs(out$ParCor)
    
    ### ena
    print("Ensembling network using ENA")
    print("=======================================")
    est_edge.ena.mat <- ena.rank(net.list = est_edge, method = "inverse.sum")
    
    ### permutate
    print("... Permutating.")
    perm.v <- perm.net.n(net.list = est_edge, method = "inverse.sum", n.perm = n.perm)
    
    ### edge confidence
    print("... Calculating confidence.")
    est_edge.ena <- setNames(melt(est_edge.ena.mat), c("Node1", "Node2", "Value"))
    est_edge.ena <- est_edge.ena[est_edge.ena$Node1 < est_edge.ena$Node2,]
    est_edge.ena[,1] <- gene.index[est_edge.ena[,1]]
    est_edge.ena[,2] <- gene.index[est_edge.ena[,2]]
    
    percentile <- ecdf(perm.v)
    est_edge.ena$Prob <- round(1 - percentile(est_edge.ena$Value), 5)
    
    return(est_edge.ena)
}


