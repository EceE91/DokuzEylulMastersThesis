## Download and install the package
install.packages("igraph")

## Load package
library(igraph)

setwd("C:/Users/Ece/Desktop")
# http://www.ats.ucla.edu/stat/r/faq/snplot.htm
# http://cneurocvs.rmki.kfki.hu/igraphbook/igraphbook-foreign.html

xlist <- read.graph("graph.txt", format="ncol",directed = FALSE)
 
length(V(xlist))
[1] 10579
length(E(xlist))
[1] 200091

# plot.igraph(xlist)
# str(xlist)

xlist.PageRank <- page_rank(xlist, vids = V(xlist), directed=FALSE, damping = 0.85)$vector      # returns A numeric vector with the PageRank scores
xlist.PageRank.Top100 <- as.data.frame(head(sort(xlist.PageRank, decreasing=TRUE), 100))

# to get first,second vertex id - real Entrez id - from xlist.PageRank.Top100
v1 = as.numeric(rownames(xlist.PageRank.Top100)[1])
v2 = as.numeric(rownames(xlist.PageRank.Top100)[2])

write.table(xlist.PageRank.Top100, file = "Top100PageRank.txt", append = TRUE, sep = "\t")

# find degree of each vertex
dg = degree(xlist)

degree.Top100 <- as.data.frame(head(sort(dg, decreasing=TRUE), 100))

a =intersect(rownames(xlist.PageRank.Top100), rownames(degree.Top100))
length(a)