#################################################################
###        simulate data with gaussian graphical model        ###
#################################################################

suppressMessages(library(igraph))
suppressMessages(library(mvnfast))

### set
# huge(p=1344) large(p=612) middle(p=231) moderate(p=83) small(p=44) tiny(p=17)

##################  function  ###################
simulate.ggm <- function(base.folder, size, n.sample, sigma, seed) {
      M <- ggm_2(n = n.sample, sigma = sigma, size = size, seed = seed)
      expr.path <- paste("data/", base.folder, "/", size, ".nSample", n.sample, ".sigma", sigma, ".csv", sep = "")
      truth.path <- paste("data/", base.folder, "/", size, ".truth.csv", sep = "")
      ### write output
      write.csv(M$Y, file = expr.path, row.names = FALSE)
      write.csv(M$adj, file = truth.path, row.names = FALSE)
}



##################  core function  #########################

ggm_2 <- function(n, sigma, size, seed = 0) {
      
      # Read data
      data <- read.csv(paste0("data/", size, ".csv"))
      scv <- which(data$Source == data$Target)
      if (length(scv) > 0) {
            data <- data[-scv,] # Remove self-connected vertices
      }
      map <- sort(unique(c(data$Source, data$Target)))
      p <- length(map)
      
      # Build the network
      G <- make_empty_graph(p, directed = TRUE)
      for (i in 1:dim(data)[1]) {
            G <- add_edges(G, c(which(map == data[i, 1]), which(map == data[i, 2])))
      }
      
      # Plot the network
      # if (plot) {
      #   set.seed(seed);
      #   plot(G, vertex.size = 0.5, edge.arrow.size = 0.5);
      # }
      
      # Create the adjacent matrix
      adj <- as.matrix(as_adjacency_matrix(G))
      adj <- adj + t(adj)
      
      # Create the precision matrix
      if (seed > 0) set.seed(seed)
      Omega <- matrix(0, nrow = p, ncol = p)
      ## Step 1
      if (size == "huge") {
            diag(Omega) <- 2
      } else if (size == "tiny") {
            diag(Omega) <- 1
      } else {
            diag(Omega) <- 1.1
      }
      temp <- which(adj == 1, arr.ind = TRUE)
      temp <- temp[which(temp[, 1] < temp[, 2]),]
      temp_2 <- sample(c(-1, 1), dim(temp)[1], replace = TRUE)*runif(dim(temp)[1], 0.5, 1)
      Omega[temp] <- temp_2
      Omega[cbind(temp[, 2], temp[, 1])] <- temp_2
      ### Step 2
      A <- matrix(0, nrow = p, ncol = p)
      diag(A) <- diag(Omega)
      for (j in 1:p) {
            if (sum(adj[j,]) > 0) {
                  A[j, -j] <- Omega[j, -j]/1.5/sum(abs(Omega[j, -j]))
            }
      }
      A <- (A + t(A))/2
      ### Step 3
      temp <- which(abs(A) < 0.1 & abs(A) > 0)
      if (length(temp) > 0) {
            A[temp] <- 0.1
      }
      
      # Create the covariance matrix
      Ar <- solve(A)
      # print(which(diag(Ar) < 0))
      Sigma <- matrix(0, nrow = p, ncol = p)
      for (j in 1:p) {
            for (jj in 1:p) {
                  Sigma[j, jj] <- Ar[j, jj]/sqrt(Ar[j, j]*Ar[jj, jj])
            }
      }
      
      # Generate Y
      Y <- matrix(0, nrow = n, ncol = p)
      for(i in 1:n) {
            Y[i,] <- rmvn(1, rep(0, p), Sigma, ncores = 4) + rnorm(p, 0, sigma)
      }
      colnames(Y) <- 
      row.names(Sigma) <- colnames(Sigma) <- 
      row.names(adj) <- colnames(adj) <- paste("X", map, sep = "")
      
      return(list(Sigma = Sigma, adj = adj, Y = Y))
}




# ggm <- function(n, p, mu, sigma, network, plot) {
#       
#       
#       if (network %in% c("BA", "ER", "AR")) {
#             K <- 5;
#             group_index <- floor(seq(0, p, length.out = K + 1));
#             M <- list();
#             Y <- matrix(NA, nrow = n, ncol = p);
#             # Generate Sigma
#             if (network == "AR") {
#                   ar <- 2; # Order of the Autoregressive model
#                   for (k in 1:K) {
#                         set.seed(seed);
#                         if (network == "AR") {
#                               module <- make_empty_graph(group_index[k + 1] - group_index[k]);
#                               for (i in 1:(group_index[k + 1] - group_index[k])) {
#                                     for (ii in 1:ar) {
#                                           if (i - ii > 0) {
#                                                 module <- add_edges(module, c(i, i - ii));
#                                           } else {
#                                                 module <- add_edges(module, c(i, group_index[k + 1] - group_index[k] + i - ii));
#                                           }
#                                     }
#                               }
#                         } else if (network == "BA") {
#                               module <- barabasi.game(group_index[k + 1] - group_index[k], power = ba, m = 2, directed = FALSE);
#                         } else if (network == "ER") {
#                               module <- erdos.renyi.game(group_index[k + 1] - group_index[k], p = er, directed = FALSE);
#                         }
#                         V(module)$name <- as.character((group_index[k] + 1):group_index[k + 1]);
#                         M[[k]] <- module;
#                         if (k == 1) {
#                               g <- module;
#                         } else {
#                               g <- union(g, module, byname = "auto");
#                         }
#                   }
#                   if (k > 1) {
#                         for (k in 1:(k - 1)) {
#                               for (kk in (k + 1):K) {
#                                     for (i in 1:1) {
#                                           g <- g + edge(as.character(sample((group_index[k] + 1):group_index[k + 1], 1)), as.character(sample((group_index[kk] + 1):group_index[kk + 1], 1)))
#                                     }
#                               }
#                         }
#                   }
#             } else if (network == "BA") {
#                   ba <- 2; # Power of the Barabási-Albert model
#                   for (k in 1:K) {
#                         set.seed(seed);
#                         if (network == "AR") {
#                               module <- make_empty_graph(group_index[k + 1] - group_index[k]);
#                               for (i in 1:(group_index[k + 1] - group_index[k])) {
#                                     for (ii in 1:ar) {
#                                           if (i - ii > 0) {
#                                                 module <- add_edges(module, c(i, i - ii));
#                                           } else {
#                                                 module <- add_edges(module, c(i, group_index[k + 1] - group_index[k] + i - ii));
#                                           }
#                                     }
#                               }
#                         } else if (network == "BA") {
#                               module <- barabasi.game(group_index[k + 1] - group_index[k], power = ba, m = 5, directed = FALSE);
#                         } else if (network == "ER") {
#                               module <- erdos.renyi.game(group_index[k + 1] - group_index[k], p = er, directed = FALSE);
#                         }
#                         V(module)$name <- as.character((group_index[k] + 1):group_index[k + 1]);
#                         M[[k]] <- module;
#                         if (k == 1) {
#                               g <- module;
#                         } else {
#                               g <- union(g, module, byname = "auto");
#                         }
#                   }
#                   if (k > 1) {
#                         for (k in 1:(k - 1)) {
#                               for (kk in (k + 1):K) {
#                                     for (i in 1:1) {
#                                           g <- g + edge(as.character(sample((group_index[k] + 1):group_index[k + 1], 1)), as.character(sample((group_index[kk] + 1):group_index[kk + 1], 1)))
#                                     }
#                               }
#                         }
#                   }
#             } else if (network == "ER") {
#                   er <- 0.02; # Probability of the Erdos-Rényi model
#                   for (k in 1:K) {
#                         module <- erdos.renyi.game(group_index[k + 1] - group_index[k], p = er, directed = FALSE);
#                         V(module)$name <- as.character((group_index[k] + 1):group_index[k + 1]);
#                         M[[k]] <- module;
#                         if (k == 1) {
#                               g <- module;
#                         } else {
#                               g <- union(g, module, byname = "auto");
#                         }
#                   }
#                   if (k > 1) {
#                         for (k in 1:(k - 1)) {
#                               for (kk in (k + 1):K) {
#                                     for (i in 1:1) {
#                                           g <- g + edge(as.character(sample((group_index[k] + 1):group_index[k + 1], 1)), as.character(sample((group_index[kk] + 1):group_index[kk + 1], 1)))
#                                     }
#                               }
#                         }
#                   }
#             }
#             A <- matrix(0, nrow = p, ncol = p);
#             diag(A) <- 1;
#             for (i in 1:length(E(g))) {
#                   A[as.numeric(ends(g, i)[1]), as.numeric(ends(g, i)[2])] <- sample(c(-1, 1), 1)*runif(1, 0.5, 1);
#                   A[as.numeric(ends(g, i)[2]), as.numeric(ends(g, i)[1])] <- A[as.numeric(ends(g, i)[1]), as.numeric(ends(g, i)[2])];
#             }
#       } else {
#             if (network == "huge") {
#                   A <- huge;
#             } else if (network == "large") {
#                   A <- large;
#             } else if (network == "middle") {
#                   A  <- middle;
#             } else if (network == "moderate") {
#                   A <- moderate;
#             } else if (network == "small") {
#                   A <- small;
#             } else if (network == "tiny") {
#                   A <- tiny
#             }
#             g <- graph_from_adjacency_matrix(A != 0);
#             if (network == "huge") {
#                   diag(A) <- 2;
#             } else {
#                   diag(A) <- max(A)*1.5;
#             }
#             p <- dim(A)[2];
#             Y <- matrix(NA, nrow = n, ncol = p);
#       }
#       for (j in 1:p) {
#             if (sum(A[, j]) != 1) {
#                   A[j, -j] <- A[j, -j]/sum(abs(A[j, -j]))/1.5;
#             }
#       }
#       A <- (A + t(A))/2;
#       Ar <- solve(A);
#       # print(which(diag(Ar) < 0))
#       Sigma <- matrix(0, nrow = p, ncol = p);
#       for (j in 1:p) {
#             for (jj in 1:p) {
#                   Sigma[j, jj] <- Ar[j, jj]/sqrt(Ar[j, j]*Ar[jj, jj]);
#             }
#       }
#       
#       if (plot) {
#             par(mfrow = c(1, 1));
#             plot(g, vertex.size = 0, vertex.label = NA, edge.arrow.size = 0);
#             if (network %in% c("BA", "ER", "AR")) {
#                   par(mfrow = c(1, K));
#                   for (k in 1:K) {
#                         set.seed(seed);
#                         plot(M[[k]], vertex.size = 0, vertex.label = NA, edge.arrow.size = 0);
#                   }
#                   par(mfrow = c(1, 1));
#             }
#       }
#       
#       # Generate X
#       for(i in 1:n) {
#             Y[i,] <- rmvn(1, mu, Sigma, ncores = 4) + rnorm(p, 0, sigma);
#       }
#       if (network %in% c("large", "huge", "middle", "moderate", "small", "tiny")) {
#             colnames(Y) <- colnames(A);
#       }
#       
#       
#       # Save the results
#       return(list(g = g, Sigma = Sigma, A = A, Y = Y))
# }

