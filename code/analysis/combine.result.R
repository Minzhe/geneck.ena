#####################################################
###                combine result                 ###
#####################################################

res.files <- list.files(path = "data/ggm/auc/")
n <- length(res.files)
auc.table <- data.frame(matrix(nrow = n, ncol = 10, dimnames = list(1:n, c("size", "n_sample", "sigma2", "GeneNet", "NS", "glasso", "glassosf", "pcacmi", "space", "ena"))))
for (i in 1:n) {
      print(i)
      res.file <- res.files[i]
      size <- strsplit(res.file, ".", fixed = TRUE)[[1]][2]
      n.sample <- as.numeric(gsub("nSample", "", strsplit(res.file, ".", fixed = TRUE)[[1]][3]))
      sigma2 <- as.numeric(gsub(".RData", "", strsplit(res.file, "sigma")[[1]][2]))
      load(paste("data/ggm/auc/", res.file, sep = ""))
      auc.table[i,] <- c(size, n.sample, sigma2, round(res$AUC, 5))
}
write.csv(auc.table, file = "data/ggm/auc/summary.csv", row.names = FALSE)
