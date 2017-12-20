#############################################
###             replicate.R               ###
#############################################

setwd("/project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena")
suppressMessages(source("code/simulation/simulate.network.R"))
suppressMessages(source("code/analysis/compare.net.roc.fun.R"))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

########################  function  ##############################
auc.replicate <- function(rounds, method, size, n.sample, noise, verbose = TRUE, plot = TRUE, seed = 0) {
      auc.rep <- data.frame(matrix(nrow = rounds, ncol = 8, dimnames = list(1:rounds, c("GeneNet", "NS", "GLasso", "GLasso-sf", "PCACMI", "SPACE", "BayesianGLasso", "ENA"))))
      for (i in 1:rounds) {
            cat("i: ", i, "; seed: ", seed, "\n", sep = "")
            simulate.network(base.folder = "temp", method = method, size = size, noise = noise, n.sample = n.sample, seed = seed)
            if(seed > 0) seed <- seed + 1
            res.tmp <- compare.net.roc(base.folder = "temp", size = size, n.sample = n.sample, sigma = noise, verbose = verbose, plot = FALSE, save = FALSE, store.temp = FALSE)
            auc.rep[i,] <- auc.tmp <- round(res.tmp$AUC, 5)
            print(auc.tmp); cat("\n")
      }
      
      res <- list(auc.rep = auc.rep,
                  simu = method,
                  size = size,
                  n.sample = n.sample,
                  noise = noise)
      if (plot) {
            plot.box(param = res, save = TRUE)
      }
      save(res, file = paste("data/replicate/", res$simu, ".", res$size, ".nSample", res$n.sample, ".sigma", res$noise, ".RData", sep = ""))
      return(res)
}


plot.box <- function(param, save = TRUE) {
      auc.data <- melt(param$auc.rep, variable.name = "method", value.name = "AUC")
      if (param$simu == "simulation") {
            setting <- paste("(", param$size, ", ", param$n.sample, " samples, ", "noise ", param$noise, ")", sep = "")
      } else if (param$simu == "ggm") {
            setting <- paste("(GGM, ", param$size, ", ", param$n.sample, " samples, ", "noise ", param$noise, ")", sep = "")
      }
      title <- paste("AUC of ROC curves of different methods\n", setting, sep = "")
      ### box plot
      p <- ggplot(auc.data, aes(x = method, y = AUC, fill = method)) + 
            geom_boxplot() +
            # geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1)
            geom_jitter(shape = 16, position = position_jitter(0.2)) +
            ggtitle(title) +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
      plot(p)
      
      if (save) {
            plot.path <- paste("data/replicate/boxplot/", param$simu, ".", param$size, ".nSample", param$n.sample, ".sigma", param$noise, ".png", sep = "")
            ggsave(plot.path, width = 6, height = 4)
      }
}


########################  parse arguments  ########################
parser <- ArgumentParser(description = "This pipline is to compare gene network construction methods based on their auc of ROC curve in replications, and plot the result in boxplot.")
parser$add_argument("rounds", type = "integer", help = "number of rounds to run")
parser$add_argument("simu", type = "character", help = "method to simulate data, whether simulation or ggm")
parser$add_argument("size", type = "character", help = "size of network, e.g. moderate")
parser$add_argument("samp", type = "integer", help = "number of samples when simulating network data")
parser$add_argument("noise", type = "double", help = "noise when simulating network data")
parser$add_argument("-p", "--plot", action = "store_true", help = "whether to plot or not")
parser$add_argument('-v', "--verbose", action = "store_true", help = "whether verbose printing")
parser$add_argument("-s", "--seed", default = 0, type = "integer", help = "seed to set, positive number (0 means random)")

args <- parser$parse_args()
print(args)
res <- auc.replicate(rounds = args$rounds,
                     method = args$simu,
                     size = args$size,
                     n.sample = args$samp,
                     noise = args$noise,
                     plot = args$plot,
                     verbose = args$verbose,
                     seed = args$seed)