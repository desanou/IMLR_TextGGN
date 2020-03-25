rm(list=ls())
library(Matrix)
library(igraph)
library(rags2ridges)
source("~/FUNCs.R")

# load the text data corresponding to three data periods, and transform them into tf-idf matrices
load("allstagetext.RDATA") #newest text data of all 3 stages

# select intersecting terms in all three data periods
intersectSeveral <- function(...) { Reduce(intersect, list(...)) } 

voc <- intersectSeveral(out1$dm$word,out2$dm$word,out3$dm$word)
a1 <- sapply(voc, function(x) out1$dm$freq[which(out1$dm$word==x)], simplify = TRUE)
a2 <- sapply(voc, function(x) out2$dm$freq[which(out2$dm$word==x)], simplify = TRUE)
a3 <- sapply(voc, function(x) out3$dm$freq[which(out3$dm$word==x)], simplify = TRUE)
b <- sort(a1+a2+a3, decreasing = TRUE)
# select first 100 most frequent temrs intersecting in all periods
voc2  <- b[1:20]

# choose the number of documents per period
n <- 100
M <- t(out3$m)
M <- M[,voc2]
Y3 <- sweep(M,MARGIN=2,log(n/apply(M, 2,nnzero)),`*`)

M <- t(out2$m)
M <- M[,voc2]
Y2 <- sweep(M,MARGIN=2,log(n/apply(M, 2,nnzero)),`*`)

M <- t(out1$m)
M <- M[,voc2]
Y1 <- sweep(M,MARGIN=2,log(n/apply(M, 2,nnzero)),`*`)

Y.all <- list(Y1,Y2,Y3)
fit=BayesFusedGlassoMix(Y.all,iternum = 100,ex.data =NULL, maxclustern = 10, k0 = 1, a_lambda = 1, b_lambda = 0.2, burnin = 20, nmc=120, aa=n, bb=0.05)

m=pcor(symm(solve(fit$S.all[[1]])))
rownames(m)=colnames(m) <- voc2
m[abs(m)<0.3] = 0
diag(m) <- 0
m1 <- m[which(rowSums(m)!=0),which(colSums(m)!=0)]
diag(m1) <- 1
net <- igraph::graph.adjacency(m1,mode="undirected",weighted=TRUE,diag=FALSE)
V(net)$color <- "gray"

par(mar= c(.3,1.2,.3,1.2))
V(net)$shape="circle"
V(net)$label.cex=1.1
plot.igraph(net,vertex.label=V(net)$name,layout=layout_with_fr,
            edge.color="black",edge.width=10*E(net)$weight )

