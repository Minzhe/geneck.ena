#################################################################
###    read the edge defination and simulate the regulation   ###
#################################################################

## size should be one of c("tiny","small","moderate","middle","large", "huge")
## sample can be any positive number.
## noise  can be any positive number (for our study, we should try 0.25, 0.5, 0.75, 1, 1.5 and 2)

network.simulation <- function(size, sample, noise) {
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
    dat$regulation <- round(sign(rbinom(N, 1, 0.5) - 0.5) * rnorm(N, 0.5, 0.2) , 2)
    dat$type <- sign(dat$regulation)
    
    write.csv(dat, paste("data/truth ", size, ".csv", sep = ""), row.names = F, quote = F)
    
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
                              expr[as.character(source), ] * dat[dat$Target == i, 3] ,
                              dat[dat$Target == i, 3] %*% expr[as.character(source), ])
                expr[as.character(i), ] <-
                    as.numeric(val > 0.0) * 2 - 1 + rnorm(nSamp, 0, sigma)
                known <- append(known, i)
                unknown <- setdiff(unknown, i)
            }
        }
    }
    
    expr <- t(expr)
    write.csv(round(expr, 2), paste("data/", size, ".nSamp", nSamp, ".Sigma", noise, ".csv", sep = ""), row.names = FALSE)
}

#### Call function ##

network.simulation(size = "huge", sample = 500, noise = 1.0)

#### Call function in loop to simulate all senarios ##

SIZE <- c("tiny", "small", "moderate", "middle", "large", "huge")
NP <- c(20, 50, 100, 200, 500, 1000)
NOISE <- c(0.25, 0.5, 0.75, 1.0, 1.5, 2.0)

set.seed(1234)

for (size in SIZE)
    for (sample in NP)
        for (noise in NOISE)
            network.simulation(size = size, sample = sample, noise = noise)
################