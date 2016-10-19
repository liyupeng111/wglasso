# Warnings: you may get different results on each run, 
# but the conclusion (MCC.wglass > MCC.glasso) should be the same.

# load libraries

library(glasso)
library(igraph)
library(huge)


# read gene expression data
# gene expression file is downloaded from additional file 2
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC545783/
expr<-read.table(file="exprall.txt", header=T, sep=" ", stringsAsFactors=F)

# preprocessing
ids<-expr[,3:4]
expr<-expr[,-(1:4)]
expr<-aggregate(expr, list(Gene = ids[,2]), mean)
rownames(expr)<-expr[,1]
expr<-expr[,-1]

genes<-rownames(expr)
p<-nrow(expr)
n<-ncol(expr)


# correlation matrix (after log and nonparanormal transformation)
s<-cor(huge.npn(t(log(expr))))


# glasso

# set penalty list for glasso
rholist=seq(from=0.10,to=0.17,by=0.01)

# glasso, select rho based on BIC
# long time to run ################
BICs<-NULL
PLs<-NULL
ENs<-NULL
for(i in 1:(length(rholist))){
  g<-glasso(s,rholist[i])
  mat <- g$wi
  diag(mat) <- 0
  mat[which(mat != 0)] <- 1
  edge.num <- sum(mat)/2
  degree <- colSums(mat)
  ENs<-c(ENs, edge.num)
  ln <- n/2*(2/p*g$loglik+sum(abs(g$wi))*rholist[i])
  BIC<- -2*ln+edge.num*log(n)
  BICs<-c(BICs, BIC)
  power.law <- power.law.fit(degree, 2)$logLik
  PLs<-c(PLs, power.law)
  print(rholist[i])
}

# result of above code ####################
#  using multiple rholist to get an optimal rho

#> BICs
#[1] 84651.27 83035.72 82104.11 81511.36 81652.52 81766.23 81870.63 82259.18
#> rholist
#[1] 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17

#> BICs
# [1]     -Inf 79589.96 77022.26 79255.63 83178.84
#> rholist
#[1] 0.05 0.10 0.15 0.20 0.25 0.30

#> rholist
#[1] 0.10 0.12 0.14 0.16 0.18 0.20
#> BICs
#[1] 79589.96 77307.14 76830.33 77478.40 78032.62 79255.63
#> rholist
#[1] 0.13 0.15
#> BICs
#[1] 76850.70 77022.26

# plot
par(mfrow=c(2,2))
plot(rholist, BICs, type="l")
plot(rholist, PLs, type="l")
plot(rholist, ENs, type="l")

# which.min(BICs)
# rho<-rholist[which.min(BICs)]

# glasso, using rho=0.13 based on BIC

rho<-0.13
g<-glasso(s,rho)
mat <- g$wi
colnames(mat)<-genes
rownames(mat)<-genes
diag(mat) <- 0
mat[which(mat != 0)] <- 1
expr.graph <- graph.adjacency(mat, mode="undirected")

# plot network
plot(expr.graph, vertex.label=NA, vertex.size=1)

# graph info
# expr.graph
# IGRAPH UN-- 739 18162 -- 




# wglasso

# get prior information from AraNet
# http://www.functionalnet.org/aranet/download.html
arbnet<-read.table("AraNet.v1.join.txt", stringsAsFactors=F)

# prepare prior matrix
arbnet<-arbnet[,c(1,2,27)]
colnames(arbnet)<-c("V1","V2","weight")
arbg<-graph.data.frame(arbnet, directed=FALSE)
a<-E(arbg)$weight
b<-0.5*(max(a)-a)/(max(a)-min(a))+0.5
E(arbg)$weight<-b
gene<-genes
arbgene<-V(arbg)$name
addname<-setdiff(gene, arbgene)
arbg<-add.vertices(arbg, length(addname), name=addname)
arbgene<-V(arbg)$name
deletename<-setdiff(arbgene, gene)
arbg<-delete.vertices(arbg, deletename)
mat.arbg<-get.adjacency(arbg, attr="weight")
mat.arbg<-as.matrix(mat.arbg)

mat.arbg[mat.arbg==0]<-1

mat.prior<-mat.arbg
prior_part<-mat.prior[genes,genes]

# set penalty list for glasso
rholist=seq(from=0.16,to=0.18,by=0.01)

# wglasso, select rho based on BIC
# long time to run ################
BICs<-NULL
PLs<-NULL
ENs<-NULL
for(i in 1:(length(rholist))){
  g<-glasso(s,rholist[i]*prior_part)
  mat <- g$wi 
  diag(mat) <- 0
  mat[which(mat != 0)] <- 1
  edge.num <- sum(mat)/2
  degree <- colSums(mat)
  ENs<-c(ENs, edge.num)
  lamda<-1-log(n)/(4*log(p))
  ln <- n/2*(2/p*g$loglik+sum(abs(g$wi))*rholist[i])
  BIC<- -2*ln+edge.num*log(n)
  BICs<-c(BICs, BIC)

  power.law <- power.law.fit(degree, 2)$logLik
  PLs<-c(PLs, power.law)
  print(rholist[i])
}

# result of above code ####################
#  using multiple rholist to get an optimal rho

#> rholist
#[1] 0.10 0.11 0.12 0.13 0.14 0.15
#> BICs
#[1] 80096.19 77969.71 76928.04 76187.48 75671.33 75406.63

#> rholist
#[1] 0.16 0.17 0.18
#> BICs
#[1] 75190.91 75241.76 75477.56



# plot
par(mfrow=c(2,2))
plot(rholist, BICs, type="l")
plot(rholist, PLs, type="l")
plot(rholist, ENs, type="l")


# which.min(BICs)
# rho<-rholist[which.min(BICs)]

# wglasso, using rho=0.16 based on BIC

rho<-0.16
g<-glasso(s,rho*prior_part)
mat <- g$wi
colnames(mat)<-genes
rownames(mat)<-genes
diag(mat) <- 0
mat[which(mat != 0)] <- 1
expr.graph.prior <- graph.adjacency(mat, mode="undirected")

# plot network
plot(expr.graph.prior, vertex.label=NA, vertex.size=1)




# compare glasso vs wglasso

CalMCC <- function(mat, true.mat) {
  # given adjacency matrices from the inferred network and the true network
  # return a list with specificity(spec), sensitivity(sens), and MCC(mcc)

  tn <- length(intersect(which(mat == 0), which(true.mat == 0)))
  fp <- length(intersect(which(mat == 1), which(true.mat == 0)))
  tp <- length(intersect(which(mat == 1), which(true.mat == 1)))
  fn <- length(intersect(which(mat == 0), which(true.mat == 1)))
  tn <- as.double(tn)
  fp <- as.double(fp)
  tp <- as.double(tp)
  fn <- as.double(fn)

  spec <- tn / (tn + fp)
  sens <- tp / (tp + fn)
  mcc.d <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- (tp * tn - fp * fn) / mcc.d
  mcc.list <- list(spec=spec, sens=sens, mcc=mcc)

  return(mcc.list)
}



arbnet<-read.table("AraNet.v1.benchmark.txt", stringsAsFactors=F)

arbnet<-arbnet[,c(1,2)]
arbg<-graph.data.frame(arbnet, directed=FALSE)
gene<-genes
arbgene<-V(arbg)$name
deletename<-setdiff(arbgene, gene)
arbg<-delete.vertices(arbg, deletename)

mat.arbg<-get.adjacency(arbg)
mat.arbg<-as.matrix(mat.arbg)

mat<-get.adjacency(expr.graph)
mat<-as.matrix(mat)
mat.prior<-get.adjacency(expr.graph.prior)
mat.prior<-as.matrix(mat.prior)
mat<-mat[rownames(mat.arbg),colnames(mat.arbg)]
mat.prior<-mat.prior[rownames(mat.arbg),colnames(mat.arbg)]

# glasso MCC
mcc.list <- CalMCC(mat, mat.arbg)
# wglasso MCC
mcc.list.prior <- CalMCC(mat.prior, mat.arbg)

# MCC results
#> mcc.list
#$spec
#[1] 0.9373035
#
#$sens
#[1] 0.1146611
#
#$mcc
#[1] 0.0328489
#
#> mcc.list.prior
#$spec
#[1] 0.9510517

#$sens
#[1] 0.3214392
#
#$mcc
#[1] 0.1838942
