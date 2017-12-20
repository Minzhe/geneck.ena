# ##################  plot roc curve  #################
# cols = c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
# plot.roc(roc = roc.GeneNet, main = "ROC", color = cols[1])
# plot.roc(roc = roc.ns, color = cols[2], add = TRUE)
# plot.roc(roc = roc.glasso, color = cols[3], add = TRUE)
# plot.roc(roc = roc.glassosf, color = cols[4], add = TRUE)
# plot.roc(roc = roc.pcacmi, color = cols[5], add = TRUE)
# plot.roc(roc = roc.space, color = cols[6], add = TRUE)
# plot.roc(roc = roc.ena, color = cols[7], add = TRUE)
# legend(0.85, 0.5, legend = c("GeneNet", "NS", "glasso", "glasso-sf", "pcacmi", "space", "ena"), col = cols, lty = 1, lwd = 3)


################  ROC with different parameters  #############
cal.ROC <- function(truth, prediction, multi = FALSE) {
      if (multi) {
            N <- dim(prediction)[3]
            roc <- matrix(NA, nrow = N, ncol = 2)   # FPR TPR
            for (i in 1:N) {
                  tab <- tabulate_error(truth != 0, prediction[,,i] != 0)
                  roc[i,1] <- tab[2,1]/sum(tab[,1]) # false positive rate
                  roc[i,2] <- tab[2,2]/sum(tab[,2]) # true positive
            }
      } else {
            cutoff <- seq(0.01, 0.99, 0.01)
            roc <- matrix(NA, nrow = length(cutoff), ncol = 2)   # FPR TPR
            for (i in 1:length(cutoff)) {
                  if (class(truth) == "matrix" & class(prediction) == "matrix") {
                        tab <- tabulate_error(upper.mat(truth) != 0, abs(upper.mat(prediction)) >= cutoff[i])
                  } else {
                        tab <- tabulate_error(truth != 0, abs(prediction) >= cutoff[i])
                  }
                  roc[i,1] <- tab[2,1]/sum(tab[,1]) # false positive rate
                  roc[i,2] <- tab[2,2]/sum(tab[,2]) # true positive
            }
      }
      roc <- rbind(c(1, 1), roc)
      roc <- rbind(roc, c(0, 0))
      return(roc)
}

tabulate_error <- function(gamma_true, gamma) {
      table <- matrix(0L, 2, 2);
      p <- length(gamma_true);
      for (i in 1:p) {
            table[gamma[i] + 1, gamma_true[i] + 1] <- table[gamma[i] + 1, gamma_true[i] + 1] + 1;
      }
      # table content
      #         true f   true t
      # pred f
      # pred t
      return (table);
}

upper.mat <- function(data.mat) {
      return(data.mat[upper.tri(data.mat)])
}

plot.roc <- function(roc, main, color = "black", add = FALSE) {
      if (!add) {
            plot(roc[, 1], roc[, 2], ylim = c(0, 1), xlim = c(0, 1), type = "l", ylab = "True positive rate", xlab = "False positive rate", main = main, cex.lab = 1.5, cex.axis = 1.5, lwd = 2, col = color)
      } else {
            lines(roc[, 1], roc[, 2], ylim = c(0, 1), xlim = c(0, 1), type = "l", cex.lab = 1.5, cex.axis = 1.5, lwd = 2, col = color)
      }
      
}

###############  AUC  ###################
cal.AUC <- function(roc) {
      area <- 0
      roc <- roc[order(-roc[,1]),]
      for (i in 1:(nrow(roc)-1)) {
            dot1 <- roc[i,]
            dot2 <- roc[i+1,]
            area <- area + trap.area(dot1, dot2)
      }
      return(area)
}

trap.area <- function(dot1, dot2) {
      x1 <- dot1[1]; y1 <- dot1[2]
      x2 <- dot2[1]; y2 <- dot2[2]
      area <- (y1 + y2) * abs(x1 - x2) / 2
      return(area)
}
