##################### grounp graph node ##################
library(igraph)
library(org.Hs.eg.db)

# pi3k.path <- read.table("http://rest.kegg.jp/link/hsa/hsa04151", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# pi3k.gene <- pi3k.path[,2]
# 
# KEGG.GP <- read.table("../Drug_combination/Report/cytoscape/KEGG.GP.net.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# idx <- c()
# for (i in 1:nrow(KEGG.GP)) {
#       if (KEGG.GP[i,1] %in% pi3k.gene & KEGG.GP[i,2] %in% pi3k.gene) {
#             idx <- c(idx,i)
#       }
# }
# pi3k.net <- KEGG.GP[idx,]
# 
# pi3k.id <- sapply(strsplit(pi3k.path[,2], ":", fixed = TRUE), "[[", 2)
# pi3k.gene <- sapply(mget(pi3k.id, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)
# 
# pi3k.net[,1] <- sapply(strsplit(pi3k.net[,1], ":", fixed = TRUE), "[[", 2)
# pi3k.net[,2] <- sapply(strsplit(pi3k.net[,2], ":", fixed = TRUE), "[[", 2)
# 
# pi3k.net[,1] <- sapply(1:nrow(pi3k.net), FUN = function(x) pi3k.gene[which(pi3k.net[x,1] == pi3k.id)])
# pi3k.net[,2] <- sapply(1:nrow(pi3k.net), FUN = function(x) pi3k.gene[which(pi3k.net[x,2] == pi3k.id)])
# 
# rm(pi3k.path, pi3k.id, i, idx, pi3k.gene)
# 
# ### network
# net <- graph_from_data_frame(d = pi3k.net, directed = TRUE)
# V(net)$community <- cluster_edge_betweenness(net)$membership

kegg.path <- read.table("http://rest.kegg.jp/link/hsa/hsa04010", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
kegg.gene <- kegg.path[,2]

KEGG.GP <- read.table("../Drug_combination/Report/cytoscape/KEGG.GP.net.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
idx <- c()
for (i in 1:nrow(KEGG.GP)) {
      if (KEGG.GP[i,1] %in% kegg.gene & KEGG.GP[i,2] %in% kegg.gene) {
            idx <- c(idx,i)
      }
}
kegg.net <- KEGG.GP[idx,]

kegg.id <- sapply(strsplit(kegg.path[,2], ":", fixed = TRUE), "[[", 2)
kegg.gene <- sapply(mget(kegg.id, org.Hs.egSYMBOL, ifnotfound = NA), "[[", 1)

kegg.net[,1] <- sapply(strsplit(kegg.net[,1], ":", fixed = TRUE), "[[", 2)
kegg.net[,2] <- sapply(strsplit(kegg.net[,2], ":", fixed = TRUE), "[[", 2)

kegg.net[,1] <- sapply(1:nrow(kegg.net), FUN = function(x) kegg.gene[which(kegg.net[x,1] == kegg.id)])
kegg.net[,2] <- sapply(1:nrow(kegg.net), FUN = function(x) kegg.gene[which(kegg.net[x,2] == kegg.id)])

rm(kegg.path, kegg.id, i, idx, kegg.gene)

### network
net <- graph_from_data_frame(d = kegg.net, directed = TRUE)
V(net)$community <- cluster_edge_betweenness(net)$membership
colnames(kegg.net) <- c("Node1", "Node2")

### output
write.csv(kegg.net, "report/MAPK.RAS.edge.csv", row.names = FALSE)
write.csv(as_data_frame(net, what = "vertices"), "report/MAPK.RAS.node.csv", row.names = FALSE)
