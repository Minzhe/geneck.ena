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
source("code/lib/BayesianGLasso.R")
source("code/lib/ENA.R")
library(reshape2)

network.ena <- function(expr.data, method = "inverse.sum", n.perm, sig.quant) {
      if (sig.quant <= 0 | sig.quant >= 1) {
            stop('Input error: parameter sig.quant for ena should be between 0 and 1.')
      }
      p <- ncol(expr.data)
      n <- nrow(expr.data)
      gene.index <- colnames(expr.data)
      expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
      
      est_edge <- list()
      ### GeneNet
      cat("Constructing using GeneNet ...")
      invisible(capture.output(pcor_est <- pcor.shrink(expr.mat)))
      invisible(capture.output(test_pcor <- network.test.edges(pcor_est)))
      est_edge.GeneNet <- matrix(NA, p, p)
      for (i in 1:nrow(test_pcor)) {
            est_edge.GeneNet[test_pcor$node1[i], test_pcor$node2[i]] <- test_pcor$prob[i]
      }
      est_edge[[1]] <- est_edge.GeneNet
      if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
      # store prob
      est_edge.prob <- test_pcor[test_pcor$node1 < test_pcor$node2, c(2,3,6)]
      colnames(est_edge.prob) <- c("Node1", "Node2", "GeneNet")
      est_edge.prob$GeneNet <- round(est_edge.prob$GeneNet, 5)
      cat(" OK\n")
      
      ### ns
      cat("Constructing using neighborhood selection ...")
      est_res.ns <- matrix(0, p, p)
	  lambda <- 0.3
      alpha <- (1-lambda)^2 * p^2
      for (k in 1:p) {
            rsp <- expr.mat[, k]
            prd <- t(expr.mat[, -k])
            lam <- sqrt(sum(rsp ^ 2)) * qnorm(alpha / (2 * p ^ 2), lower.tail = F)
            out <- l2.reg(prd, rsp, lambda = lam)
            est_res.ns[k, -k] <- out$estimate
      }
      est_edge[[2]] <- (abs(est_res.ns) + t(abs(est_res.ns))) / 2
      # store prob
      est_edge.prob <- cbind(est_edge.prob, NS = mapply(FUN = function(x,y) est_edge[[2]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      ### glasso
      cat("Constructing using GLASSO ...")
      S <- t(expr.mat) %*% expr.mat / n
      out <- glasso(S, rho = 0.1)
      est_edge[[3]] <- abs(out$wi)
      # store prob
      est_edge.prob <- cbind(est_edge.prob, GLasso = mapply(FUN = function(x,y) est_edge[[3]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      ### glasso-sf
      cat("Constructing using GLASSO-SF ...")
      out <- glasso_sf(expr.mat, alpha = 0.1)
      est_edge[[4]] <- abs(out$wi)
      # store prob
      est_edge.prob <- cbind(est_edge.prob, GLasso.sf = mapply(FUN = function(x,y) est_edge[[4]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      ### pcacmi
      cat("Constructing using PCACMI ...")
      invisible(capture.output(out <- pca_cmi(t(expr.mat), 0.03)))
      est_edge[[5]] <- abs(out$Gval)
      # store prob
      est_edge.prob <- cbind(est_edge.prob, PCACMI = mapply(FUN = function(x,y) est_edge[[5]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      # ### cmi2ni
      # cat("\nConstructing using CMI2NI ...")
      # out <- cmi2ni(t(expr.mat), 0.03)
      # est_edge[[6]] <- abs(out$Gval)
      # # store prob
      # est_edge.prob <- cbind(est_edge.prob, cmi2ni = mapply(FUN = function(x,y) est_edge[[6]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      # cat("\nOK\n")
      
      ### space
      cat("Constructing using SPACE ...")
      invisible(capture.output(out <- space.joint(expr.mat, lam1 = 0.2 * n, iter = 5)))
      est_edge[[6]] <- abs(out$ParCor)
      # store prob
      est_edge.prob <- cbind(est_edge.prob, SPACE = mapply(FUN = function(x,y) est_edge[[6]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      ### BayesianGLasso
      cat("Constructing using BayesianGLasso ...")
      a <- 10^(-2); b <- 10^(-6); iter <- 2000; burn <- 1000
      invisible(capture.output(out <- blockGLasso_s(expr.mat, iterations = iter, burnIn = burn, lambdaPriora = a, lambdaPriorb = b, verbose = FALSE)))
      est_edge[[7]] <- abs(out)
      # store prob
      est_edge.prob <- cbind(est_edge.prob, BayesianGLasso = mapply(FUN = function(x,y) est_edge[[7]][x,y], est_edge.prob$Node1, est_edge.prob$Node2))
      cat(" OK\n")
      
      names(est_edge) <- c("GeneNet", "NS", "GLasso", "GLasso-sf", "PCACMI", "SPACE", "BayesianGLasso")
      
      ### ena
      cat("Ensembling network using ENA ...")
      est_edge.ena.mat <- ena.rank(net.list = est_edge, method = method)
      cat(" OK\n")
      
      ### permutate
      cat("Permutating ...")
      perm.v <- perm.net.n(net.list = est_edge, method = method, n.perm = n.perm)
      cat(" OK\n")
      
      ### edge confidence
      cat("Calculating confidence ...")
      est_edge.ena <- setNames(melt(est_edge.ena.mat), c("Node1", "Node2", "Value"))
      est_edge.ena <- est_edge.ena[est_edge.ena$Node1 < est_edge.ena$Node2,]
      percentile <- ecdf(perm.v)
      est_edge.ena$ena <- round(1 - percentile(est_edge.ena$Value), 5)
      
      # store prob
      est_edge.prob <- merge(est_edge.prob, est_edge.ena[,-3], by = c("Node1", "Node2"))
      
      # set gene name
      est_edge.prob$Node1 <- gene.index[est_edge.prob$Node1]
      est_edge.prob$Node2 <- gene.index[est_edge.prob$Node2]
      cat(" OK\n")
      
      return(est_edge.prob)
}

############### add true label  ###################
prep.edge.roc <- function(est_edge, truth_edge) {
      est_edge.roc <- data.frame(est_edge, Truth = rep(0, nrow(est_edge)))
      for (i in 1:nrow(truth_edge)) {
            Source <- truth_edge[i,1]
            Target <- truth_edge[i,2]
            idx <- which((est_edge.roc$Node1 == Source & est_edge.roc$Node2 == Target) | (est_edge.roc$Node2 == Source & est_edge.roc$Node1 == Target))
            est_edge.roc$Truth[idx] <- 1
      }
      return(est_edge.roc)
}