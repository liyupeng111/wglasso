#######
Calau <- function(mat, true.mat){
  library(caTools)
  prediction<-NULL
  outcome<-NULL
  diag(true.mat) <- 0
  T <-nrow(true.mat)*(nrow(true.mat)-1)/2
  P <-sum(true.mat)/2
  N <-T-P
  pred<-cbind(as.vector(mat[upper.tri(mat, diag=F)]),as.vector(true.mat[upper.tri(true.mat, diag=F)]))
  
  L<-sum(pred[,1]>0)
  TPL<-sum(pred[,1]>0 & pred[,2]==1)
  p<-0
  if(L<T){
    p<-(P - TPL) / (T - L)
  }
  pred.order<-pred[order(-pred[,1]),]
  
  TPk<-0;
  FPk<-0;
  REC<-NULL;
  PREC<-NULL;
  TPR<-NULL;
  FPR<-NULL;
  for(K in 1:T){
    if(pred.order[K,1]>0 & pred.order[K,2]==1){
      TPk<-TPk+1    
    }else if(pred.order[K,1]>0 & pred.order[K,2]==0){
      FPk<-FPk+1
    }else if(pred.order[K,1]==0){
      TPk<-TPk+p
      FPk<-FPk+1-p
    }
    REC <-c(REC, TPk/P);
    PREC<-c(PREC,TPk/K);
    TPR <-c(TPR, TPk/P);
    FPR <-c(FPR, FPk/N);
  }
  
 
  auc<-trapz(FPR,TPR)
  aupr<-trapz(REC, PREC)/(1-1/P)
  au<-list(auc=auc,aupr=aupr)
  return(au)
}
######

###############################################################################
# function to calculate specificity, sensitivity, and MCC
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
  #spec <- tp / (tp + fp)
  mcc.d <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- (tp * tn - fp * fn) / mcc.d
  #mcc <- 2*spec*sens/(spec+sens)
  mcc.list <- list(spec=spec, sens=sens, mcc=mcc)

  return(mcc.list)
}

###############################################################################
GetHub <- function(degree, top=0.05) {
  # get the hub: top 5% nodes in terms of degree
  len <- length(degree)
  hub.len <- as.integer(len*top)
  y <- cbind(degree, seq(1:len))
  z <- y[order(-y[, 1]), ]
  hub <- z[, 2][1:hub.len]
  nohub <- z[, 2][-(1:hub.len)]
  hub.list <- list(hub=hub, nohub=nohub)
  return(hub.list)
}

###############################################################################
# function to run multiple glasso
RunMultiGlasso <- function(a, s, prior.mat, true.mat, 
                           rholist=seq(from=0.01,to=1,by=0.01)) {
  # given multi gaussian data, prior info matrix, true matrix, and a rho list
  # rho is the penalty parameter of graphical lasso (see glasso package)
  # code explains itself

  powerlaw <- NULL
  spec     <- NULL
  sens     <- NULL
  mcc      <- NULL
  BIC      <- NULL
  mBIC     <- NULL
  aupr     <- NULL
  auc      <- NULL

  for (rho in rholist) {
    stat.list <- RunSingleGlasso(a, s, rho, prior.mat, true.mat)
    powerlaw  <- c(powerlaw , stat.list$powerlaw )
    spec      <- c(spec     , stat.list$spec     )
    sens      <- c(sens     , stat.list$sens     )
    mcc       <- c(mcc      , stat.list$mcc      )
    BIC       <- c(BIC      , stat.list$BIC      )
    mBIC      <- c(mBIC     , stat.list$mBIC     )
    aupr      <- c(aupr     , stat.list$aupr     )
    auc       <- c(auc      , stat.list$auc      )
  }

  stat.alllist <- list(powerlaw=powerlaw, 
                       spec=spec, 
                       sens=sens, 
                       mcc=mcc, 
                       BIC=BIC, 
                       mBIC=mBIC,
                       aupr=aupr,
                       auc=auc,
                       rholist=rholist)

  return(stat.alllist)
}

###############################################################################
# function to run single glasso
RunSingleGlasso <- function(a, s, rho, prior.mat, true.mat) {
  # given multi gaussian data, prior info matrix, true matrix, and a rho
  # code explains itself

  library(glasso)
  g <- glasso(s, prior.mat*rho)
  n <- nrow(a)
  p <- nrow(s)
  
  g.mat<-g$wi
  au<-Calau(g.mat, true.mat)
  
  mat <- g$wi
  diag(mat) <- 0
  mat[which(mat != 0)] <- 1
  edge.num <- sum(mat)/2
  
  degree <- colSums(mat)
  
  true.degree <- colSums(true.mat)

  power.law <- power.law.fit(degree, 2)$logLik # fit power law, log(mle)
  mcc.list <- CalMCC(mat, true.mat)
  
  aupr<-au$aupr
  auc <-au$auc
  ln <- n/2*(2/p*g$loglik+sum(abs(g$wi*prior.mat))*rho)

  BIC<- -2*ln+edge.num*log(n)
  mBIC<- -2*ln+edge.num*log(n)+4*0.5*edge.num*log(p)
  #lambda<-1-log(n)/(4*log(p))
  #BIC<- -2*ln+edge.num*log(n)+2*log(p/edge.num)
  #BIC<- -2*ln+2*edge.num

  #BIC <- -2*ln+edge.num*log(n)
  #BIC<- -2*g$loglik+log(n)*edge.num+4*edge.num*rho*log(ncol(a))
  #BIC<- -log(det(g$wi))+sum(diag(s%*%g$wi))+log(n)/n*edge.num
  #mBIC <- BIC+max(5-mean(degree),0)*log(n)
  
  stats <- list(powerlaw=power.law, 
                spec=mcc.list$spec, 
                sens=mcc.list$sens, 
                mcc=mcc.list$mcc, 
                BIC=BIC, 
                mBIC=mBIC,
                aupr=aupr,
                auc=auc)

  return(stats)
}

###############################################################################
# function to plot the glasso results
WGlassoPlot <- function(sim.stat, t.graph) {
  # given the simulation stats from RunMultiGlasso, and true graph
  # code explains itself

  par(mfrow=c(2,2))

  plot(sim.stat$rholist, sim.stat$powerlaw, type="l")
  abline(h=power.law.fit(igraph::degree(t.graph)+1, 2)$logLik, col="red")

  plot(sim.stat$rholist, sim.stat$BIC,type="l")
  plot(sim.stat$rholist, sim.stat$aupr,type="l")
  
  plot(sim.stat$rholist, sim.stat$spec, type="l",col="blue", 
                                        xlim=c(0,1), ylim=c(0,1))
  points(sim.stat$rholist, sim.stat$sens,type="l", col="red")
  points(sim.stat$rholist, sim.stat$mcc,type="l", col="black")
}

###############################################################################
# function to generate prior infor matrix
PriorMatrixGen <- function(true.mat,
                           support.ratio=0.8) {
  # given:
  #   true matrix, 
  #   support ratio, e.g. 0.8=80% true edges are supported by prior info
  #   wrong support ratio, e.g. 0.1=10% true un-edges are mis-supported by prior
  # code explains itself

  p.mat <- true.mat
  p.mat[] <- 1
  t.edge   <- which(true.mat == 1)
  t.unedge <- which(true.mat == 0)
  p.edge <- sample(t.edge, size=length(t.edge) * support.ratio)
  p.edge <- sort(c(p.edge, sample(t.unedge, 
                                  size=length(t.edge) * (1-support.ratio))))
  p.mat[p.edge] <- runif(length(p.edge), min=0, max=1)
  p.mat<-as.matrix(forceSymmetric(p.mat))
  diag(p.mat) <- 1

  return(p.mat)
}

###############################################################################
# function to run multiple glasso_mb
RunMultiGlassoMB <- function(a, prior.mat, true.mat, core.num=2, 
                           rholist=seq(from=0.1,to=1,by=0.1)) {
  # given multi gaussian data, prior info matrix, true matrix, and a rho list
  # rho is the penalty parameter of graphical lasso (see glasso package)
  # code explains itself

  powerlaw <- NULL
  spec     <- NULL
  sens     <- NULL
  mcc      <- NULL
  BIC      <- NULL
  mBIC     <- NULL
  
  for (rho in rholist) {
    stat.list <- RunSingleGlassoMB(a, rho, prior.mat, true.mat, core.num=core.num)
    powerlaw  <- c(powerlaw , stat.list$powerlaw )
    spec      <- c(spec     , stat.list$spec     )
    sens      <- c(sens     , stat.list$sens     )
    mcc       <- c(mcc      , stat.list$mcc      )
    BIC       <- c(BIC      , stat.list$BIC      )
    mBIC      <- c(mBIC     , stat.list$mBIC     )
  }

  stat.alllist <- list(powerlaw=powerlaw, 
                       spec=spec, 
                       sens=sens, 
                       mcc=mcc, 
                       BIC=BIC, 
                       mBIC=mBIC, 
                       rholist=rholist)
                       
  return(stat.alllist)
}

###############################################################################
# function to run single glasso
RunSingleGlassoMB <- function(a, rho, prior.mat, true.mat, core.num=2) {
  # given multi gaussian data, prior info matrix, true matrix, and a rho
  # code explains itself

  library(igraph)

  betas <- GlassoParMB(a, rho, prior.mat, core.num=core.num)
  graph <- GlassoMB2igraph(a, betas)
  mat <- get.adjacency(graph)

  edge.num <- ecount(graph)  # count the edge num
  
  power.law <- power.law.fit(degree(graph), 2)$logLik # fit power law, log(mle)
  
  mcc.list <- CalMCC(mat, true.mat)
  spec <- mcc.list$spec
  sens <- mcc.list$sens
  mcc  <- mcc.list$mcc
  BIC  <- 0
  mBIC <- 0
  
  stats <- list(powerlaw=power.law, 
                spec=mcc.list$spec, 
                sens=mcc.list$sens, 
                mcc=mcc.list$mcc, 
                BIC=BIC, 
                mBIC=mBIC)

  return(stats)
}

###############################################################################
# function to run single glasso using mb method
GlassoParMB <- function(a, rho, prior.mat, core.num=2) {
  library(glmnet)
  library(doParallel)

  registerDoParallel(cores=core.num)
  genes <- colnames(a)

  b <- foreach(i=1:ncol(a), .combine='cbind', .packages='glmnet') %dopar% {
    penalty <- as.vector(prior.mat[, i])[-i]
    b.lasso <- glmnet(a[, -i], a[, i], lambda=rho, penalty.factor=penalty)
    as.matrix(b.lasso$beta)
  }

  rownames(b) <- seq(1:(ncol(a)-1))
  
  return(b)
}

###############################################################################
CorrectEdge <- function (x) {
  if (x[1] >= x[2]) {
    x[1] <- x[1]+1
  }
  
  return(x)
}

###############################################################################
# function to parse the glasso(MB) results into an igraph object
GlassoMB2igraph <- function(a, b) {
  library(igraph)
  
  node.num <- ncol(a)
  
  if (is.null(colnames(a))) {
    node.name <- seq(1:node.num)
  } else {
    node.name <- colnames(a)
  }
  
  edges.F.ind <- which(b != 0, arr.ind = TRUE)
  
  if (length(edges.F.ind) == 0) {
    mat <- matrix(0, nrow=node.num,ncol=node.num)
    g <- graph.adjacency(mat, mode="undirected")
  } else {
    edges.T.ind <- apply(edges.F.ind, 1, CorrectEdge)
    g <- graph(edges.T.ind)
    g <- as.undirected(g, mode="collapse")
  }
  
  V(g)$name<-node.name[V(g)]
  
  return(g)
}

