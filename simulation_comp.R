# simulations

# install required libraries
# install.packages("huge")
# install.packages("glasso")
# install.packages("glmnet")
# install.packages("doParallel")

# load required libraries and functions
library(huge)
library(glasso)
library(glmnet)
library(doParallel)
library(gplots)
library(igraph)
source("sim.func.R")

###single simulation

# generate a scale-free network and corresponding multi gaussian distn data
sample.size <- 50
node.num <- 100
sim.huge <- huge.generator(n=sample.size, d=node.num, graph="scale-free")

# multi gaussian distn data and the covariance matrix
a <- sim.huge$data
s <- cor(a)
#s<-cor( huge.npn(a))


# plot the network
#huge.plot(sim.huge$theta)


# the adjacency matrix
t.mat <- as.matrix(sim.huge$theta) # true network matrix
p.mat <- PriorMatrixGen(t.mat, support.ratio=0.7) # prior network matrix

# only used for generating figure 1
p.mat2<-p.mat
p.mat2[p.mat2<1]<-2
p.mat2[p.mat2<2]<-0
comb.mat<-t.mat+p.mat2

t.graph <- graph.adjacency(t.mat, mode="undirected")
E(t.graph)$width <- 2
E(t.graph)$color <- "black"
plot(t.graph,vertex.label=NA,vertex.size=5) # true network

# transfer the adjacency matrix to an igraph object
comb.graph <- graph.adjacency(comb.mat, mode="undirected",weighted=TRUE)
E(comb.graph)$color <- "gray"
E(comb.graph)$width <- 2
E(comb.graph)$color[E(comb.graph)$weight==2] <- "red"
E(comb.graph)$color[E(comb.graph)$weight==3] <- "black"
E(comb.graph)$width[E(comb.graph)$weight==3] <- 2
#E(comb.graph)$width[E(comb.graph)$weight==1] <- 2
plot(comb.graph,vertex.label=NA,vertex.size=5) # prior network

# heatmap of gene expression
heatmap.2(a, density.info="none", trace="none", labRow="",labCol="", dendrogram="none")

# penalty parameter
rholist<-seq(from=0.1,to=0.9,by=0.01)

# glasso, tune penalty parameter
noprior.mat <- t.mat
noprior.mat[] <- 1  # no prior information
sim.stat <- RunMultiGlasso(a, s, noprior.mat, t.mat, rholist=rholist)

# wglasso, tune penalty parameter
p.sim.stat <- RunMultiGlasso(a, s, p.mat, t.mat,rholist=rholist)

# glasso, select the best penalty parameter
g<-glasso(s,rholist[which.max(sim.stat$mcc)])

# only used for generating figure
mat<-g$wi
colnames(mat)<-rownames(s)
rownames(mat)<-rownames(s)
diag(mat)<-0
mat[which(mat != 0)] <-1

mat2<-mat
mat2[which(mat2 != 0)]<-2
comb.mat<-t.mat+mat2

# transfer the adjacency matrix to an igraph object
comb.graph <- graph.adjacency(comb.mat, mode="undirected",weighted=TRUE)
E(comb.graph)$color <- "gray"
E(comb.graph)$width <- 2
E(comb.graph)$color[E(comb.graph)$weight==2] <- "red"
E(comb.graph)$color[E(comb.graph)$weight==3] <- "black"
E(comb.graph)$width[E(comb.graph)$weight==3] <- 2
#E(comb.graph)$width[E(comb.graph)$weight==1] <- 2
plot(comb.graph,vertex.label=NA,vertex.size=5)

# wglasso, select the best penalty parameter
g<-glasso(s,rholist[which.max(p.sim.stat$mcc)]*p.mat)

# only used for generating figure
mat<-g$wi
colnames(mat)<-rownames(s)
rownames(mat)<-rownames(s)
diag(mat)<-0
mat[which(mat != 0)] <-1

mat2<-mat
mat2[which(mat2 != 0)]<-2
comb.mat<-t.mat+mat2

# transfer the adjacency matrix to an igraph object
comb.graph <- graph.adjacency(comb.mat, mode="undirected",weighted=TRUE)
E(comb.graph)$color <- "gray"
E(comb.graph)$width <- 2
E(comb.graph)$color[E(comb.graph)$weight==2] <- "red"
E(comb.graph)$color[E(comb.graph)$weight==3] <- "black"
E(comb.graph)$width[E(comb.graph)$weight==3] <- 2
#E(comb.graph)$width[E(comb.graph)$weight==1] <- 2
plot(comb.graph,vertex.label=NA,vertex.size=5)


### figure 2

par(mar=c(5, 4, 4, 6) + 0.1)
plot(sim.stat$rholist, sim.stat$mcc, type="l",col="blue", 
     xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",axes=F)
points(p.sim.stat$rholist, p.sim.stat$mcc,type="l", col="red")
points(sim.stat$rholist[which.max(sim.stat$mcc)],sim.stat$mcc[which.max(sim.stat$mcc)],col="blue",pch=17)
points(p.sim.stat$rholist[which.max(p.sim.stat$mcc)],p.sim.stat$mcc[which.max(p.sim.stat$mcc)],col="red",pch=17)

axis(2,col="black",las=1)  ## las=1 makes horizontal labels
mtext("MCC",side=2,line=2.5)
box()

par(new=TRUE)
plot(sim.stat$rholist, sim.stat$BIC, type="l",lty=2, col="blue", axes=F, xlab="", ylab="")
points(p.sim.stat$rholist,p.sim.stat$BIC,type="l",lty=2,col="red")
points(sim.stat$rholist[which.max(sim.stat$mcc)],sim.stat$BIC[which.min(sim.stat$BIC)],col="blue",pch=17)
points(p.sim.stat$rholist[which.max(p.sim.stat$mcc)],p.sim.stat$BIC[which.min(p.sim.stat$BIC)],col="red",pch=17)

mtext("BIC",side=4,line=4) 
axis(4,las=1)

axis(1,pretty(range(sim.stat$rholist),10))
mtext("Rho",side=1,line=2.5)  


## Add Legend
legend("topright",legend=c("MCC (glasso)","MCC (wglasso)","BIC (glasso)","BIC (wglasso)"),
       text.col=c("blue","red","blue","red"),
       lty=c(1,1,2,2),
       col=c("blue","red","blue","red"),
       bty = "n",y.intersp=0.6, cex=0.5)


#####


# systematic comparision of glasso and wglasso

## take time !!!
prec<-c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
sample<-c(20, 50, 100, 300)

sim.stat.all<-NULL
node.num <- 100
repeats.num<-100
rholist<-seq(from=0.1,to=0.9,by=0.01)

sim.stat.mat<-NULL

# graph simulation
for (p.ratio in prec){
  for (sample.size in sample){
    for(i in seq(1:repeats.num)){
      sim.huge <- huge.generator(n=sample.size, d=node.num, graph="scale-free")
      a <- sim.huge$data
      s <- cor(a)
      t.mat <- as.matrix(sim.huge$theta)
      
      if(p.ratio==0){
        noprior.mat <- t.mat
        noprior.mat[] <- 1
        sim.stat <- RunMultiGlasso(a, s, noprior.mat, t.mat, rholist=rholist)
      }else{
        p.mat <- PriorMatrixGen (t.mat, support.ratio=p.ratio)
        sim.stat <- RunMultiGlasso(a, s, p.mat, t.mat,rholist=rholist)
      }
      rep<-rep(i,length(rholist))
      samplenum<-rep(sample.size,length(rholist))
      precision<-rep(p.ratio,length(rholist))
      sim.stat.mat.temp<-cbind(do.call(cbind.data.frame, sim.stat),rep,samplenum,precision)
      sim.stat.mat<-rbind(sim.stat.mat, sim.stat.mat.temp)
    }
  }
}

# intead of running above code, load a result I ran before
load("20140818.RData")


mat.mcc.mean<-matrix(NA,nrow=length(prec),ncol=length(sample))
mat.mcc.ciw<-mat.mcc.mean
mat.bic.mean<-mat.mcc.mean
mat.bic.ciw<-mat.mcc.mean
mat.diff.mean<-mat.mcc.mean
mat.diff.ciw<-mat.mcc.mean

for (i in 1:length(prec)){
  for (j in 1:length(sample)){
    mcc.mean.temp<-NULL
    bic.mean.temp<-NULL
    diff.temp<-NULL
    for(k in 1:repeats.num){     
      sim.stat.temp<-sim.stat.mat[which(sim.stat.mat$samplenum==sample[j] & 
                                            sim.stat.mat$prec==prec[i] &
                                            sim.stat.mat$rep==k),]
      mcc.mean.temp<-c(mcc.mean.temp, max(sim.stat.temp$mcc, na.rm=T))
      bic.mean.temp<-c(bic.mean.temp, sim.stat.temp$mcc[which.min(sim.stat.temp$BIC)])
      diff.temp<-mcc.mean.temp-bic.mean.temp
    }
    mat.mcc.mean[i,j]<-mean(mcc.mean.temp, na.rm=T)
    mat.mcc.ciw[i,j]<-1.96*sd(mcc.mean.temp, na.rm=T)/sqrt(repeats.num)
    mat.bic.mean[i,j]<-mean(bic.mean.temp, na.rm=T)
    mat.bic.ciw[i,j]<-1.96*sd(bic.mean.temp, na.rm=T)/sqrt(repeats.num)
    mat.diff.mean[i,j]<-mean(diff.temp, na.rm=T)
    mat.diff.ciw[i,j]<-1.96*sd(diff.temp, na.rm=T)/sqrt(repeats.num)
  }
}

cols<-c("black","red","blue","purple")

# plot(prec, mat.mcc.mean[,1])


# figure 3, MCC
# main="Performance comparison of wglasso and glasso with p=100",
plotCI(prec, y= mat.mcc.mean[,1], uiw=mat.mcc.ciw[,1], gap = 0, type="l", pch=".", col=cols[1], ylim=c(0,1), 
        xlab="Precision ratio", ylab="maxMCC")

for (i in 2:(length(sample))){
    plotCI(prec, y= mat.mcc.mean[,i], uiw=mat.mcc.ciw[,i], gap = 0, type="l", col=cols[i], add=T)
}
legend("bottomright", horiz=F, cex=0.55, text.width=0.08,legend=c("n=20","n=50","n=100","n=300"),lty=rep(1,4),col=cols)

# figure 3_1, minBIC(MCC)
# main="Performance comparison of wglasso and glasso with p=100",
plotCI(prec, y= mat.bic.mean[,1], uiw=mat.bic.ciw[,1], gap = 0, type="l", pch=".", col=cols[1], ylim=c(0,1),
       xlab="Precision ratio", ylab="BICMCC")

for (i in 2:(length(sample))){
  plotCI(prec, y= mat.bic.mean[,i], uiw=mat.bic.ciw[,i], gap = 0, type="l", col=cols[i], add=T)
}
legend("bottomright", horiz=F, cex=0.55, text.width=0.08,legend=c("n=20","n=50","n=100","n=300"),lty=rep(1,4),col=cols)


# figure 4
# main="Comparison of two tuning parameter methods, BIC and MCC, with p=100", 
plotCI(prec, y= mat.diff.mean[,1], uiw=mat.diff.ciw[,1], gap = 0, type="l",pch=".", col=cols[1], ylim=c(0,0.15), 
       xlab="Precision ratio", ylab="maxMCC-BICMCC")

for (i in 2:(length(sample))){
    plotCI(prec, y= mat.diff.mean[,i], uiw=mat.diff.ciw[,i], gap = 0, type="l", col=cols[i], add=T)
}
legend("topleft", horiz=F, cex=0.55, text.width=0.08,legend=c("n=20","n=50","n=100","n=300"),lty=rep(1,4),col=cols)

