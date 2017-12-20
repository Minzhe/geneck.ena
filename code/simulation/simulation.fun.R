#################################################################
###    read the edge defination and simulate the regulation   ###
#################################################################

## size should be one of c("tiny","small","moderate","middle","large", "huge")
## sample can be any positive number.
## noise  can be any positive number (for our study, we should try 0.25, 0.5, 0.75, 1, 1.5 and 2)

#####################  function  ###########################
network.simulation <- function(base.folder, size, sample, noise) {
    dat <- read.csv(paste("data/", size, ".csv", sep = ""))
    
    #dim(dat)
    #head(dat)
    
    dat <- dat[dat$Source != dat$Target, ]
    dat[dat$Source > dat$Target, ] <- dat[dat$Source > dat$Target, 2:1]
    dat[dat$Source > dat$Target, ]
    
    Gene <- union(dat$Source, dat$Target)
    nGene <- length(Gene)
    
    print(c(size, dim(dat)[1], nGene, sample, noise))
    
    ######################
    
    N <- dim(dat)[1]
    dat$regulation <- round(sign(rbinom(N, 1, 0.5) - 0.5) * rnorm(N, 0.5, 0.2), 2)
    dat$type <- sign(dat$regulation)
    
    write.csv(dat, paste("data/", base.folder, "/truth ", size, ".csv", sep = ""), row.names = F, quote = F)
    
    #### initiate the matrix expr to store the simulated expression level ####
    
    nSamp <- sample
    expr <- matrix(0, nGene, nSamp)
    rownames(expr) <- Gene
    colnames(expr) <- paste("S", 1:nSamp, sep = "")
    
    sigma <- noise
    
    unknown <- unique(dat$Target)
    known <- setdiff(Gene, unknown)
    
    for (i in known){
        expr[as.character(i), ] <- rbinom(nSamp, 1, 0.5) * 2 + rnorm(nSamp, 0, sigma) - 1
    }
    
    while (length(unknown) > 0) {
        for (i in unknown) {
            # i <- unknown[1]
            source <- as.character(dat$Source[dat$Target == i])
            if (all(source %in% known)) {
                #print(dat[dat$Target == i,]); cat("\n")
                val <- ifelse(rep(length(source) == 1, nSamp),
                              expr[as.character(source), ] * dat[dat$Target == i, 3],
                              dat[dat$Target == i, 3] %*% expr[as.character(source), ])
                expr[as.character(i), ] <- as.numeric(val > 0.0) * 2 - 1 + rnorm(nSamp, 0, sigma)
                known <- append(known, i)
                unknown <- setdiff(unknown, i)
            }
        }
    }
    
    expr <- t(expr)
    colnames(expr) <- paste("X", colnames(expr), sep = "")
    write.csv(round(expr, 2), paste("data/", base.folder, "/", size, ".nSample", nSamp, ".sigma", noise, ".csv", sep = ""), row.names = FALSE)
}

###################  adjacency matrix  ####################
edge.to.adjacency <- function(edge.file) {
      file.name <- strsplit(edge.file, ".", fixed = TRUE)[[1]]
      adj.name <- paste(file.name[-length(file.name)], "matrix.csv")
      
      dat <- read.csv(edge.file)[,1:3]
      node.name <- sort(unique(c(dat$Source, dat$Target)))
      adj.mat <- matrix(0, length(node.name), length(node.name), dimnames = list(node.name, node.name))
      for (i in 1:nrow(dat)) {
            r.idx <- which(node.name == dat$Source[i])
            c.idx <- which(node.name == dat$Target[i])
            adj.mat[r.idx,c.idx] <- adj.mat[c.idx,r.idx] <- dat$regulation[i]
      }
      write.csv(adj.mat, adj.name)
}
