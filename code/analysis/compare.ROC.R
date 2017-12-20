############################################################
###                    Compare.ROC.R                     ###
############################################################

setwd("/project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena")
suppressMessages(library(argparse))
source("code/analysis/compare.net.roc.fun.R")

#############  1. Parse comandline argument  #################
parser <- ArgumentParser(description = "This pipline is to compare gene network construction methods based on their auc of ROC curve.")
parser$add_argument("-p", "--path", type = "character", help = "base data folder")
parser$add_argument("-s", "--size", type = "character", help = "size of network")
parser$add_argument("-v", "--sigma2", type = "double", help = "sigma2 when simulating network data")
parser$add_argument("-n", "--num", type = "integer", help = "number of samples when simulating network data")

args <- parser$parse_args()
folder <- args$path
size <- args$size
sigma2 <- args$sigma2
n <- args$num

res <- compare.ggm.net.roc(base.folder = folder, size = size, n.sample = n, sigma = sigma2, verbose = TRUE, plot = TRUE, save = TRUE, store.temp = TRUE)
