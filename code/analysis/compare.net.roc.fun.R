############################################################
###                compare.net.roc.R                     ###
############################################################

suppressMessages(library(corpcor))
suppressMessages(library(GeneNet))
suppressMessages(library(CDLasso))
suppressMessages(library(glasso))
source("code/lib/glasso_SF.R")
source("code/lib/PCA_CMI.R")
source("code/lib/CMI2NI.R")
suppressMessages(library(space))
source("code/lib/BayesianGLasso.R")
source("code/ena.R")
source("code/analysis/roc.fun.R")


compare.net.roc <- function(base.folder, size, n.sample, sigma, verbose = TRUE, plot = TRUE, save = TRUE, store.temp = TRUE) {
      
      data.path <- paste("data/", base.folder, "/", size, ".nSample", n.sample, ".sigma", sigma, ".csv", sep = "")
      truth.path <- paste("data/", base.folder, "/", size, ".truth.csv", sep = "")
      rdata.path <- paste("data/", base.folder, "/auc/res.", size, ".nSample", n.sample, ".sigma", sigma, ".RData", sep = "")
      cache.path <- paste("data/", base.folder, "/cache/cache.", size, ".nSample", n.sample, ".sigma", sigma, ".txt", sep = "")
      
      ### catch terminal output
      if (store.temp) sink(cache.path)
      cat("Working with ", data.path, " ...\n", sep = "")
      cat("=================================================\n\n")
     
      ########## read data #############
      net.data <- read.data(data.path = data.path, truth.path = truth.path)
      expr.data <- net.data$expr.data; expr.mat <- net.data$expr.mat;
      truth.edges <- net.data$truth.edges; truth.mat <- net.data$truth.mat
      p <- net.data$p; n <- net.data$n
      
      #########  method  ##########
      AUC <- rep(NA, 8); names(AUC) <- c("GeneNet", "NS", "GLasso", "GLasso-sf", "PCACMI", "SPACE", "BayesianGLasso", "ENA")
      
      ### GeneNet
      cat("------------- GeneNet -------------\n")
      invisible(capture.output(pcor_est <- pcor.shrink(expr.mat)))
      invisible(capture.output(test_pcor <- network.test.edges(pcor_est)))
      est_edge.GeneNet <- matrix(NA, p, p)
      for (i in 1:nrow(test_pcor)) {
            est_edge.GeneNet[test_pcor$node1[i], test_pcor$node2[i]] <- test_pcor$prob[i]
      }
      if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
      roc.GeneNet <- cal.ROC(truth = truth.mat, prediction = est_edge.GeneNet, multi = FALSE)
      AUC[1] <- cal.AUC(roc.GeneNet)
      
      ### ns
      cat("------------- NS -------------\n")
      lambda <- c(seq(0.01, 0.1, 0.01), seq(0.12, 0.5, 0.02), seq(0.54, 0.98, 0.04))
      est_edge.ns <- array(NA, dim=c(p, p, length(lambda)))
      for (i in 1:length(lambda)) {
            alpha <- (1-lambda[i])^2 * p^2
            est_res.ns <- matrix(0, p, p)
            for (k in 1:p) {
                  rsp <- expr.mat[, k]
                  prd <- t(expr.mat[, -k])
                  lam <- sqrt(sum(rsp ^ 2)) * qnorm(alpha / (2 * p ^ 2), lower.tail = F)
                  out <- l2.reg(prd, rsp, lambda = lam)
                  est_res.ns[k, -k] <- out$estimate
            }
            est_edge.ns[,,i] <- (abs(est_res.ns) + t(abs(est_res.ns))) / 2
            if (verbose) print.temp.auc(param = c(lambda[i], lam), truth = truth.mat, est_edge = est_edge.ns[,,i])
      }
      roc.ns <- cal.ROC(truth = truth.mat, prediction = est_edge.ns, multi = TRUE)
      AUC[2] <- cal.AUC(roc.ns)
      
      
      ### glasso
      cat("------------- glasso -------------\n")
      S <- t(expr.mat) %*% expr.mat / n
      rho <- c(seq(0.01, 0.1, 0.01), seq(0.12, 0.5, 0.02), seq(0.54, 0.8, 0.04))
      est_edge.glasso <- array(NA, dim=c(p, p, length(rho)))
      for (i in 1:length(rho)) {
            out <- glasso(s = S, rho = rho[i])
            est_edge.glasso[,,i] <- abs(out$wi)
            if (verbose) print.temp.auc(param = rho[i], truth = truth.mat, est_edge = est_edge.glasso[,,i])
      }
      roc.glasso <- cal.ROC(truth = truth.mat, prediction = est_edge.glasso, multi = TRUE)
      AUC[3] <- cal.AUC(roc.glasso)
      
      
      ### glasso-sf
      cat("------------- glassosf -------------\n")
      rho <- c(seq(0.01, 0.1, 0.01), seq(0.12, 0.5, 0.02), seq(0.54, 0.8, 0.04))
      est_edge.glassosf <- array(NA, dim = c(p, p, length(rho)))
      for (i in 1:length(rho)) {
            out <- glasso_sf(expr.mat, alpha = rho[i])
            est_edge.glassosf[,,i] <- abs(out$wi)
            if (verbose) print.temp.auc(param = rho[i], truth = truth.mat, est_edge = est_edge.glassosf[,,i])
      }
      roc.glassosf <- cal.ROC(truth = truth.mat, prediction = est_edge.glassosf, multi = TRUE)
      AUC[4] <- cal.AUC(roc.glassosf)
      
      
      ### pcacmi
      cat("------------- pcacmi -------------\n")
      invisible(capture.output(out <- pca_cmi(t(expr.mat), 0.03)))
      est_edge.pcacmi <- abs(out$Gval)
      roc.pcacmi <- cal.ROC(truth = truth.mat, prediction = est_edge.pcacmi, multi = FALSE)
      AUC[5] <- cal.AUC(roc.pcacmi)
      
      
      ### space
      cat("------------- space -------------\n")
      lambda <- c(seq(0.01, 0.1, 0.01), seq(0.12, 0.5, 0.02), seq(0.54, 0.98, 0.04))
      est_edge.space <- array(NA, dim = c(p, p, length(lambda)))
      for (i in 1:length(lambda)) {
            invisible(capture.output(out <- space.joint(expr.mat, lam1 = lambda[i] * n, iter = 5)))
            est_edge.space[,,i] <- abs(out$ParCor)
            if (verbose) print.temp.auc(param = lambda[i], truth = truth.mat, est_edge = est_edge.space[,,i])
      }
      roc.space <- cal.ROC(truth = truth.mat, prediction = est_edge.space, multi = TRUE)
      AUC[6] <- cal.AUC(roc.space)
      
      
      ### BayesianGLasso
      cat("------------- BayesianGLasso -------------\n")
      a <- 10^(-2); b <- 10^(-6); iter <- 2000; burn <- 1000
      est_edge.bayesglasso <- blockGLasso_s(expr.mat, iterations = iter, burnIn = burn, lambdaPriora = a, lambdaPriorb = b, verbose = FALSE)
      roc.bayesianglasso <- cal.ROC(truth = truth.mat, prediction = est_edge.bayesglasso, multi = FALSE)
      AUC[7] <- cal.AUC(roc.bayesianglasso)
      
      
      ### ena
      cat("------------- ena -------------\n")
      if (verbose) {
            est_edge.ena <- network.ena(expr.data = expr.data, method = "inverse.sum", n.perm = 10, sig.quant = 0.99)
      } else {
            invisible(capture.output(est_edge.ena <- network.ena(expr.data = expr.data, method = "inverse.sum", n.perm = 10, sig.quant = 0.99)))
      }
      est_edge.ena <- prep.edge.roc(est_edge = est_edge.ena, truth_edge = truth.edges)
      roc.ena <- cal.ROC(truth = est_edge.ena$Truth, prediction = est_edge.ena$ena, multi = FALSE)
      AUC[8] <- cal.AUC(roc.ena)
      
      est_edge.all = list(est_edge.GeneNet = est_edge.GeneNet,
                          est_edge.ns = est_edge.ns,
                          est_edge.glasso = est_edge.glasso,
                          est_edge.glassosf = est_edge.glassosf,
                          est_edge.pcacmi = est_edge.pcacmi,
                          est_edge.space = est_edge.space,
                          est_edge.bayesglasso = est_edge.bayesglasso,
                          est_edge.ena = est_edge.ena)
      
      ##########  roc  ###########
      if (plot) {
            cols = c("#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#AA4499")
            plot.roc(roc = roc.GeneNet, main = "ROC", color = cols[1])
            plot.roc(roc = roc.ns, color = cols[2], add = TRUE)
            plot.roc(roc = roc.glasso, color = cols[3], add = TRUE)
            plot.roc(roc = roc.glassosf, color = cols[4], add = TRUE)
            plot.roc(roc = roc.pcacmi, color = cols[5], add = TRUE)
            plot.roc(roc = roc.space, color = cols[6], add = TRUE)
            plot.roc(roc = roc.ena, color = cols[7], add = TRUE)
            plot.roc(roc = roc.bayesianglasso, color = cols[8], add = TRUE)
            legend(0.85, 0.5, legend = c("GeneNet", "NS", "GLasso", "GLasso-sf", "PCACMI", "SPACE", "BayesianGLasso", "ENA"), col = cols, lty = 1, lwd = 3)
      }
      
      ### close sink
      if (store.temp) sink()
      
      res <- list(AUC = AUC,
                  roc.all = list(roc.GeneNet = roc.GeneNet,
                                 roc.ns = roc.ns,
                                 roc.glasso = roc.glasso,
                                 roc.glassosf = roc.glassosf,
                                 roc.pcacmi = roc.pcacmi,
                                 roc.space = roc.space,
                                 roc.ena = roc.ena,
                                 roc.bayesianglasso = roc.bayesianglasso),
                  net.all = est_edge.all)
      
      ### save result
      if (save) save(res, file = rdata.path)
      
      return(res)
}


#################### utility function #######################

read.data <- function(data.path, truth.path) {
      ### truth
      truth.mat <- as.matrix(read.csv(truth.path))
      truth.edges <- get.truth.edges(truth = truth.mat)
      ### data
      expr.data <- read.csv(data.path)
      expr.data <- expr.data[,colnames(truth.mat)]
      p <- ncol(expr.data)
      n <- nrow(expr.data)
      expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
      
      return(list(expr.data = expr.data, expr.mat = expr.mat, p = p, n = n, truth.edges = truth.edges, truth.mat = truth.mat)) 
}

get.truth.edges <- function(truth) {
      row.names(truth) <- colnames(truth)
      truth.edge <- setNames(melt(truth), c("Node1", "Node2", "Value"))
      truth.edge$Value[truth.edge$Value != 0] <- 1
      truth.edge <- truth.edge[truth.edge$Value == 1,]
      return(truth.edge)
}

print.temp.auc <- function(param, truth, est_edge) {
      temp.roc <- cal.ROC(truth = truth, prediction = est_edge, multi = FALSE)
      cat(param, cal.AUC(temp.roc), "\n")
}