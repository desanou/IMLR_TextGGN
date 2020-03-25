
#library(statmod)
library(Matrix)
library(mvtnorm)
library(MASS)
library(rags2ridges)
library(statmod)
library("clusterGeneration")
#################################################################

Folds <- function (n, nfold = 5) {
  split(sample(1:n), rep(1:nfold, length = n))
}

#################################################################
LL.fused <- function(Ylist, s.num,Ks,pis, cMeans, cPs){
  LL <- N <- numeric()
  c=0
  p <- dim(Ylist[[1]])[2] 
  for (s in 1:s.num){
    N[s] <- dim(Ylist[[s]])[1]
    cWeights <- numeric()
    for (k in 1:Ks[s]){
      # calculate individual contributions on log-scale
      c=c+1
      cWeights <- cbind(cWeights,log(pis[[s]][k]) - 0.5*(p * log(2 * pi) - log(det(cPs[[c]])) + 
                                                           apply(sweep(Ylist[[s]], 2, cMeans[[c]], "-"), 1, function(Yi, cPk){ t(Yi) %*% cPk %*% Yi }, cPk=cPs[[c]]))); 
      
    }
    # subtract maximum value (to avoid infinity when taking exponential)
    maxW <- max(cWeights);
    cWeights <- cWeights - maxW;
    LL[s] <- sum(log(rowSums(exp(cWeights)))[log(rowSums(exp(cWeights)))!=-Inf])  + N[s] * maxW;
    #cWeights.list[[s]] <- t(matrix(apply(cWeights, 1, function(cWi){ exp(cWi) / sum(exp(cWi)) }),nrow=Ks[s])); 
  }
  return(LL)
}
###############################################################
LL <- function(Y, pis, cPs, cMeans){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  K <- length(pis)  
  cWeights <- matrix(NaN, n,K)
  
  for (k in 1:K){
    # calculate individual contributions on log-scale
    cWeights[,k] <- log(pis[k]) - p * log(2 * pi) / 2 + determinant(cPs[[k]])$modulus / 2 - 
      apply(sweep(Y, 2, cMeans[k,], "-"), 1, function(Yi, cPk){ t(Yi) %*% cPk %*% Yi }, cPk=cPs[[k]]) / 2; 
  }
  # subtract maximum value (to avoid infinity when taking exponential)
  maxW <- max(cWeights);
  cWeights <- cWeights - maxW;
  ll <- sum(log(rowSums(exp(cWeights)))[log(rowSums(exp(cWeights)))!=-Inf]) + n * maxW;
  
  return(ll)
  
}

#############################################

lili <- function(Y,mean,cP){
  sum( - 0.5*(p * log(2 * pi) - log(det(cP)) + 
                apply(sweep(Y, 2, mean, "-"), 1, function(Yi){ t(Yi) %*% cP %*% Yi })));
}

assess.errors <- function(cov.hat,mu.hat, mutruth, covtruth) {
  if (is.matrix(mu.hat)) K=dim(mutruth)[1] else {K=1 ; mu.hat=matrix(mu.hat, nrow=1); mutruth=matrix(mutruth, nrow=1);cov.hat=list(cov.hat)}
  #if (is.list(mu.hat)) {K=length(mu.hat);mu.hat=t(simplify2array(mu.hat))}
  sparseP <- list()
  mu.err = cov.err <-  TP <- FP <-TN <- FN <- FPR <- TPR <- f1_score <- PPV <- rep(0,K)
 
  a <- sort(rowSums(mu.hat), index.return=T)$ix # to account for the switching problem
  mu.hat <- matrix(mu.hat[a,] ,nrow = K)
 
  
  M <- cov.hat
  
  for (k in 1:K) {
    mu.err[k] <- f.norm(mutruth[k,],mu.hat[k,])
    cov.hat[[k]] <- M[[a[k]]] # to account for the switching problem
    cov.err[k] <- f.norm(covtruth[[k]],cov.hat[[k]]) 
    sparseP[[k]] <- sparsify(symm(cov.hat[[k]]),threshold="localFDR", FDRcut = 0.5,verbose = FALSE)$sparsePrecision
    TP[k] <- length(intersect(which(sparseP[[k]]!=0), which(covtruth[[k]]!=0)))
    FP[k] <- length(intersect(which(sparseP[[k]]!=0), which(covtruth[[k]]==0)))
    FN[k] <- length(intersect(which(sparseP[[k]]==0), which(covtruth[[k]]!=0)))
    TN[k] <- length(intersect(which(sparseP[[k]]==0), which(covtruth[[k]]==0)))
    f1_score[k] <- (2*TP[k])/(2*TP[k]+FP[k]+FN[k]) # F1 score
    TPR[k] <- TP[k]/(TP[k]+FN[k]) # sensitivity
    FPR[k] <- FP[k]/(FP[k]+TN[k]) # 1- SPC
    PPV[k] <- TP[k]/(TP[k]+FP[k]) # precision
  }
  
  return(list(thet_err=cov.err,cov_hat=cov.hat, mu_err=mu.err,TP=TP,FP=FP,FN=FN,TN=TN,
              fscore=f1_score, TPR=TPR, FPR=FPR,PPV=PPV))
}


assess.errors2 <- function(sparseP, cov.hat,mu.hat, mutruth, covtruth) {
  if (is.matrix(mu.hat)) K=dim(mutruth)[1] else {K=1 ; mu.hat=matrix(mu.hat, nrow=1); mutruth=matrix(mutruth, nrow=1);cov.hat=list(cov.hat);sparseP=list(sparseP)}
  #if (is.list(mu.hat)) {K=length(mu.hat);mu.hat=t(simplify2array(mu.hat))}
  #sparseP <- list()
  mu.err = cov.err <-  TP <- FP <-TN <- FN <- FPR <- TPR <- f1_score <- PPV <- rep(0,K)
  
  a <- sort(rowSums(mu.hat), index.return=T)$ix # to account for the switching problem
  mu.hat <- matrix(mu.hat[a,] ,nrow = K)
  
  
  M <- cov.hat
  
  for (k in 1:K) {
    mu.err[k] <- f.norm(mutruth[k,],mu.hat[k,])
    cov.hat[[k]] <- M[[a[k]]] # to account for the switching problem
    cov.err[k] <- f.norm(covtruth[[k]],cov.hat[[k]]) 
    sparseP[[k]] <- sparseP[[a[k]]] #sparsify(symm(cov.hat[[k]]),threshold="localFDR", FDRcut = 0.5,verbose = FALSE)$sparsePrecision
    TP[k] <- length(intersect(which(sparseP[[k]]!=0), which(covtruth[[k]]!=0)))
    FP[k] <- length(intersect(which(sparseP[[k]]!=0), which(covtruth[[k]]==0)))
    FN[k] <- length(intersect(which(sparseP[[k]]==0), which(covtruth[[k]]!=0)))
    TN[k] <- length(intersect(which(sparseP[[k]]==0), which(covtruth[[k]]==0)))
    f1_score[k] <- (2*TP[k])/(2*TP[k]+FP[k]+FN[k]) # F1 score
    TPR[k] <- TP[k]/(TP[k]+FN[k]) # sensitivity
    FPR[k] <- FP[k]/(FP[k]+TN[k]) # 1- SPC
    PPV[k] <- TP[k]/(TP[k]+FP[k]) # precision
  }
  
  return(list(thet_err=cov.err,cov_hat=cov.hat, mu_err=mu.err,TP=TP,FP=FP,FN=FN,TN=TN,
              fscore=f1_score, TPR=TPR, FPR=FPR,PPV=PPV))
}

###################################
#### SIMULATE DATA
###################################
sim.Gmix <- function(n, cov, mu, pi1=c(0.5,0.5),pi2=c(1/3,1/3,1/3) ){
  # simulate with p=10
  Y <- Z <- list()
  Y[[1]] <- mvrnorm(n,mu[[1]],cov[[1]])
  dat <- sim.mvmnorm(n,p,pi1,mu[[2]],list(cov[[2]],cov[[3]]))
  Y[[2]] <- dat[[1]]
  Z[[2]] <- dat[[2]]
  dat <- sim.mvmnorm(n,p,pi2,mu[[3]],list(cov[[4]],cov[[5]],cov[[6]]))
  Y[[3]] <- dat[[1]]
  Z[[3]] <- dat[[2]]
  return(list(Y=Y,Z=Z))
}


###############################################################
#generates sample from alpha based on Escobar and West (1995).#
###############################################################
randalpha<-function(aa,bb,n,k,alpha0,itr=5){
  ### ======================================================================
  ### Purpose: generates sample from alpha based on Escobar and West (1995).
  ### ======================================================================
  ### Arguments: aa and bb (hypermarameters); n (sample size of data), k (nu
  ### -mber of occupied clusters); alpha0 (initial alpha); itr (iterations).
  ### ======================================================================
  alpha<-rep(1,itr)
  for (l in 1:itr){
    xx<-rbeta(1,alpha0+1,n)
    weights<-c((aa+k-1)/(bb-log(xx)),n)/sum(c((aa+k-1)/(bb-log(xx)),n))
    pisample<-sample(1:2,1,p=weights)
    if (pisample==1) alpha[l]<-rgamma(1,aa+k,scale=1/(bb-log(xx)))
    else alpha[l]<-rgamma(1,aa+k-1,scale=1/(bb-log(xx)))
  }
  return(mean(alpha))
}


########################
## generate a GMM data##
########################
sim.mvmnorm <- function(n,p,prob,mu,sd){
  
  
  ## Purpose: generate n samples from a GMM
  ## ----------------------------------------------------------------------
  ## Arguments: 
  ## n: nr. of samples
  ## p: dimension
  ## mu: K*p matrix of mean values for a GMM with K components.
  ## sd: K*p*p array consisting K matrix of p*p each presenting covariance matrix of a component.
  ## prob: vector of length K, each element representing mixture probability corresponding to one component.
  ## ----------------------------------------------------------------------
  
  k <- length(prob)
  Y <- matrix(0,n,p)
  Z <- sample(1:k,n,replace=TRUE,p=prob) # generate component memberships
  for (i in 1:n){
    Y[i,] <- mvrnorm(1,mu[Z[i],],sd[[Z[i]]]) # given Z, generate Y from multivariate normal dist.
  }
  return(list(Y,Z))
}

############################
## Calculate Forinobis loss#
############################
f.norm <- function (P,P_hat) {
  d <- P-P_hat
  r <- sqrt(sum(sapply(d, function(X){ sum(X^2) }, simplify=TRUE)))
  return(r)
}
###########################################
## Bayesian covariance estimation function#
###########################################
BayesGLasso_Columnwise <- function(S,n,Sig,C,a_lambda,b_lambda,burnin,nmc){
  
  p=dim(S)[1];
  indmx <- matrix(1:p^2 , ncol=p)
  upperind <- indmx[upper.tri(indmx,diag=FALSE)]
  lowerind <- indmx[lower.tri(indmx,diag=FALSE)]
  
  
  C_save <- Sig_save <- array(0,c(p,p,nmc)); 
  lambda_save <- rep(0,nmc);
  tau <- matrix(0,p,p);
  
  ind_noi_all <- matrix(0,p-1,p);
  
  for (i in 1:p){
    if (i==1){ind_noi = 2:p} else if (i==p) {ind_noi = 1:(p-1)} else  ind_noi = c(1:(i-1),(i+1):p);i
    ind_noi_all[,i] = ind_noi;
  }
  
  apost <- a_lambda + p*(p+1)/2; 
  
  for (iter in 1: (burnin+nmc)){ 
    
    # Sample lambda 
    bpost <- b_lambda + sum(abs(C))/2;    
    lambda <- rgamma(1,apost,rate=bpost)
    
    # sample tau off-diagonal        
    Cadjust = pmax(abs(C[upperind]),10^-6);        
    lambda_prime = lambda^2;  mu_prime = pmin(lambda/Cadjust,10^12);
    
    
    tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
    tau[upperind] <- tau_temp;
    tau[lowerind] <- tau_temp;
    
    
    # sample Sig and C
    for (i in 1:p){
      ind_noi <- ind_noi_all[,i];
      tau_temp <- tau[ind_noi,i];
      
      Sig11 <- Sig[ind_noi,ind_noi]; Sig12 <- Sig[ind_noi,i];
      invC11 <- Sig11 - Sig12%*%t(Sig12)/Sig[i,i];
      
      Ci <- (S[i,i]+lambda)*invC11+diag(1/tau_temp);
      Ci[mapply(is.infinite, Ci)] <- 10000000000
      Ci_chol = chol(Ci);
      mu_i = -solve(Ci,S[ind_noi,i]);
      beta = mu_i+ solve(Ci_chol,rnorm(p-1));
      C[ind_noi,i] = beta;
      C[i,ind_noi] = beta;
      gam = rgamma(1,n/2+1,rate=1/solve((S[i,i]+lambda),2)) 
      
      C[i,i] = gam+t(beta)%*%invC11%*%beta;
      
      
      # Below updating Covariance matrix according to one-column change of precision matrix
      invC11beta = invC11%*%beta;
      
      Sig[ind_noi,ind_noi] = invC11+invC11beta%*%t(invC11beta)/gam;
      Sig12 = -invC11beta/gam;
      Sig[ind_noi,i] = Sig12;
      Sig[i,ind_noi] = t(Sig12);
      Sig[i,i] = 1/gam;
      
    }
    
    
    if (iter >burnin){           
      Sig_save[,,iter-burnin] <- Sig 
      C_save[,,iter-burnin] <- C
      lambda_save[iter-burnin] <- lambda
    };
  }# end iter
  
  return(list(Sig_save=Sig_save,C_save=C_save,lambda_save=lambda_save))
  
}


BayesGlassoMix_baseFunc <- function(Y,iternum = 100,ex.data = NULL, maxclustern = 20, k0 = 1, 
                                    a_lambda = 1, b_lambda = 0.2, burnin = 20, nmc=120, aa=2*dim(Y)[1],bb=0.05) {
  # Bayesian inference for the Mixture of Graphical models with Lasso penalty as a counterpart for frequentist Mixture of 
  # Graphical models with the Ridge penalties
  # Bayesian Stagewise method: look at each stage separately. Given z of each stage, optimize penalized log-likelihood
  # with lasso penalties
  # Bayesian Fused method: Not possible for now as there are no appropriate priors, that can capture fused type errors, around.
  
  H<-function(x,clusters) {
    ### ============================================================================================================
    ### Purpose: calculates similarity of a data point with other data points a component: integrated IMLR "I-IMLR".
    ###=============================================================================================================
    ### Arguments: x (data point indicator); clusters (vector of cluster indicators)
    ###=============================================================================================================
    c=0
    output<-rep(0,length(clusters))
    for (f in clusters){
      c=c+1
      ind <- which(z==f)
      
      #if (is.element(x,ind)&&length(ind)>1) ind<-ind[-which(ind==x)]
      #if (length(ind)==1) {output[c]<-1+SIM[x,ind]}
     # else{
        #bb<-SIM[ind,ind]
        #T<-median(bb[lower.tri(bb)])
        # when SIM is a binary similarity matrix:
        nStar=sum(SIM[x,ind])
        output[c] <- (nStar+1)*nStar
        
        # when SIM is a distance matrix use function below (as in text data analysis where SIM is distance of enteries e.g. SIM[2,9]=abs(2-9))
        #T<-sum(SIM[ind,ind][lower.tri(SIM[ind,ind], diag = FALSE)])/(1/2*length(ind)*(length(ind)-1))#Threshold
        #nStar<-length(which(SIM[x,ind]<=T))
        #output[c]<-(1/(min(sum(SIM[x,ind]),10^(-10))))*(nStar)+1      # used when distances are given
      #}
    }
    return(output)
  }
  H1<-function(x,clusters) {
    ### ==========================================================================================================
    ### Purpose: is used instead of H when there are no additional similarity/distance data, and calculates number
    ### of existing members in a component (i.e. component size): "IMLR".
    ###===========================================================================================================
    ### Arguments: x (data point indicator); clusters (vector of cluster indicators).
    ###===========================================================================================================
    
    return(N[clusters])
  }
  
  if(is.null(ex.data)) {  H<-H1  } else {SIM<-ex.data ; rm(ex.data)   }
  
  
  
  n<-dim(Y)[1]
  p<-dim(Y)[2]
  alpha<-rep(2*n,iternum+1)
  LMBDs <- rep(1,maxclustern)
  #########
  # Inits #
  #########
  # Latent variables
  z<-sample(1:maxclustern, size=n, replace = TRUE)  
  # component sizes
  A=as.data.frame(table(z))
  N<-rep(0,maxclustern)
  N[1:maxclustern %in% A$z]=A$Freq
  
  # mixture covariance matrices and means
  SIGs <-  replicate(maxclustern,genPositiveDefMat(p)$Sigma, simplify = FALSE)
  MUs <- sapply(SIGs, function(x) {mvrnorm(1,rep(0,p),x)}, simplify = FALSE)
  # likelihoods
  DNORMs <- sapply(1:maxclustern, function(nr) dmvnorm(Y, mean = MUs[[nr]], sigma =SIGs[[nr]]))
  ac1=ac2=ac3=ac4=0
  
  for (m in 1:iternum){
    for (i in 1:n){
      
      if (N[z[i]]>1 && length(N[N!=0])<maxclustern) {   #### Find non-singelton obserations and try to open a new cluster
        
        newsig <- genPositiveDefMat(p)$Sigma ;  newmu <- mvrnorm(1,rep(0,p),newsig);
        a1=dmvnorm(Y[i,], mean = newmu, sigma =newsig);
        a0<-DNORMs[i,z[i]];
        
        if (a0==0 || runif(1)<(alpha[m]/(n-1))*(a1/a0)){
          N[z[i]]<-N[z[i]]-1; z[i]<-which(N==0)[1]; N[z[i]]<-1;  SIGs[[z[i]]]<-newsig; MUs[[z[i]]]<-newmu;  ac1=ac1+1
        } #end if 
        
      } #end if (N[z[i]]>1)
      
      else if (N[z[i]]==1) { ### if i-th observarion is a singelton send it to an alreay occupied cluster with probability prop to N[i]/n-1.
        BB<-N; BB[z[i]]<-0; index<-which(BB!=0);
        b0<- DNORMs[i,z[i]];
        prob<-c(H(i,index)/sum(H(i,index)))
        o<-0
        repeat{
          o=o+1
          v<-sample(1:length(index),size=1,prob=prob) 
          sss<-prob[v] 
          z1<-index[v]
          if(DNORMs[i,z1]!=0 || sss==1 || o>100) {
            break
          }
        }#end repeat
        b1<- DNORMs[i,z1]
        if (b0==0 || runif(1)<((n-1)/alpha[m])*b1/b0) {
          SIGs[[z[i]]]=MUs[[z[i]]]=0
          N[z[i]]<-0
          ac2=ac2+1
          z[i]<-z1
          N[z1]<-N[z1]+1 
        }# end if
        
        
      }# end else if (N[z[i]]==1)      
    }# end for i
    
    
    ######################################################################
    # Ignore singelton obs, exchange obs among those clusters with size>1#
    ######################################################################
    
    for (i in 1:n){ 
      ind2<-which(N>1)
      if (length(ind2)>1){
        
        d=H(i,ind2)*DNORMs[i,ind2] 
        if (sum(d)!=0){
          o<-0
          repeat{
            o=o+1
            w<-d/sum(d)
            newz<-sample(ind2,size=1, prob=w)
            if (DNORMs[i,newz]!=0 || w[which(ind2==newz)]==1 || o>100) 
            {break}
          }#end repeat
          #c1<-H(i,newz)*dnorm(y[i],beta[newz,]%*%X[i,],sqrt(sigma[newz]))
          #c1<-dnorm(y[i],beta[newz,]%*%X[i,],sqrt(sigma[newz]))
          
          #if (c0==0 || runif(1)<(c1)/(c0)){
          N[z[i]]<-N[z[i]]-1
          z[i]<-newz
          N[newz]<-N[newz]+1
          ac3=ac3+1
          #} #end if
        }#end if
      }#end if
    }# end for i
    
    
    
    
    
    
    ###############################################################
    # Parameter sampling from posteriors
    ###############################################################
    
    for (k in which(N>0)) {
      Yk=matrix(Y[z==k,], ncol=p);
      nk=dim(Yk)[1]
      A=BayesGLasso_Columnwise(t(Yk-MUs[[k]])%*%(Yk-MUs[[k]]),nk,SIGs[[k]],solve(SIGs[[k]]),a_lambda,b_lambda,burnin,nmc)
      
      SIGs[[k]]= apply(A$Sig_save,1:2, mean)
      MUs[[k]] <- mvrnorm(1,colMeans(Yk)*nk/(k0+nk),SIGs[[k]]/(k0+nk))
      DNORMs[,k] <-  dmvnorm(Y, mean = MUs[[k]], sigma =SIGs[[k]])
    } #end for k
    alpha[m+1]<-randalpha(aa,bb,n,length(N[N>0]),alpha0=alpha[m],itr = 5)  
    
    
  } # end for m
  
  # remove & merge clusters that are less populated
  for (i in 1:n){
    if (N[z[i]]< 0.08*n) {
      N[z[i]]=N[z[i]]-1; SIGs[[z[i]]]=0; MUs[[z[i]]]=0; DNORMs[i,z[i]]=0
      z[i] <- which(N> 0.08*n)[which.max(DNORMs[i,which(N> 0.08*n)])];  N[z[i]]= N[z[i]]+1
    }
  }
  # last update
  namana=which(N>0.08*n)
  sparseCs = Cs = replicate(length(namana), list())
  c=0
  for (k in namana) {
    c=c+1
    Yk=matrix(Y[z==k,], ncol=p);
    nk=dim(Yk)[1]
    A=BayesGLasso_Columnwise(t(Yk-MUs[[k]])%*%(Yk-MUs[[k]]),nk,SIGs[[k]],solve(SIGs[[k]]),a_lambda,b_lambda,100,1000)
    SIGs[[k]]= apply(A$Sig_save,1:2, mean)
    Cs[[c]]= apply(A$C_save,1:2, mean)
    quant_low=apply(A$C_save, 1:2, function(x) quantile(x, 0.05))
    quant_up=apply(A$C_save, 1:2, function(x) quantile(x, 0.95))
    #x=sum(abs(C.mean[[i]]))/max(apply(M$C_save,1:2, function(x) sum(abs(x))))
    #sparsify using credible intervals
    sparseCs[[c]] <- Cs[[c]]
    sparseCs[[c]][mapply(function(x,y) {findInterval(0,c(x,y))==1 }, x=quant_low, y=quant_up, SIMPLIFY = TRUE)]=0
    MUs[[k]] <- mvrnorm(1,colMeans(Yk)*nk/(k0+nk),SIGs[[k]]/(k0+nk))
    DNORMs[,k] <-  dmvnorm(Y, mean = MUs[[k]], sigma =SIGs[[k]])
  } 
  
  
  Out<-list(ac1=ac1,ac2=ac2,ac3=ac3, ac4=ac4,groups=N,Y=Y,z=z, SIGs=SIGs, Cs=Cs, sparseCs=sparseCs, MUs=MUs, ALPHAs=alpha)
  return(Out)
  
  print("##################################################################################################")
}



BayesFusedGlassoMix <- function(Y.all,iternum = 100,ex.data = NULL, maxclustern = 20, k0 = 1, 
                                a_lambda = 1, b_lambda = 0.2, burnin = 20, nmc=120, aa=dim(Y.all[[1]])[1],bb=0.05) {
  # Bayesian inference for the Mixture of Graphical models with Lasso penalty as a counterpart for frequentist Mixture of 
  # Graphical models with the Ridge penalties
  # Bayesian Stagewise method: look at each stage separately. Given z of each stage, optimize penalized log-likelihood
  # with lasso penalties
  # Bayesian Fused method: Not possible for now as there are no appropriate priors, that can capture fused type errors, around.
  
  H<-function(x,clusters) {
    ### ============================================================================================================
    ### Purpose: calculates similarity of a data point with other data points a component: integrated IMLR "I-IMLR".
    ###=============================================================================================================
    ### Arguments: x (data point indicator); clusters (vector of cluster indicators)
    ###=============================================================================================================
    c=0
    output<-rep(0,length(clusters))
    for (f in clusters){
      c=c+1
      ind <- which(z==f)
      
      #if (is.element(x,ind)&&length(ind)>1) ind<-ind[-which(ind==x)]
      #if (length(ind)==1) {output[c]<-1+SIM[x,ind]}
      # else{
      #bb<-SIM[ind,ind]
      #T<-median(bb[lower.tri(bb)])
      # when SIM is a binary similarity matrix:
      nStar=sum(SIM[x,ind])
      output[c] <- (nStar+1)*nStar
      
      # when SIM is a distance matrix use function below (as in text data analysis where SIM is distance of enteries e.g. SIM[2,9]=abs(2-9))
      #T<-sum(SIM[ind,ind][lower.tri(SIM[ind,ind], diag = FALSE)])/(1/2*length(ind)*(length(ind)-1))#Threshold
      #nStar<-length(which(SIM[x,ind]<=T))
      #output[c]<-(1/(min(sum(SIM[x,ind]),10^(-10))))*(nStar)+1      # used when distances are given
      #}
    }
    return(output)
  }
  H1<-function(x,clusters) {
    ### ==========================================================================================================
    ### Purpose: is used instead of H when there are no additional similarity/distance data, and calculates number
    ### of existing members in a component (i.e. component size): "IMLR".
    ###===========================================================================================================
    ### Arguments: x (data point indicator); clusters (vector of cluster indicators).
    ###===========================================================================================================
    
    return(N[clusters])
  }
  
  
  
  
  G <- length(Y.all)
  ns <- sapply(Y.all, dim)[1,]
  ps <- sapply(Y.all, dim)[2,]
  N.all = z.all = alpha.all = SIGs.all = MUs.all = DNORMs.all  <- list()
  # mixture covariance matrices and means
  
  #########
  # Inits #
  #########
  # Latent variables
  for (g in 1:G){
    Y <- Y.all[[g]]
    n<-ns[g]
    p<-ps[g]
    z.all[[g]] <-sample(1:maxclustern, size=n, replace = TRUE)  
    # component sizes
    R=as.data.frame(table(z.all[[g]]))
    N.all[[g]] <- rep(0,maxclustern)
    N.all[[g]][1:maxclustern %in% R$Var1]=R$Freq
    SIGs.all[[g]] <-  replicate(maxclustern,genPositiveDefMat(p)$Sigma, simplify = FALSE)
    MUs.all[[g]] <- sapply(SIGs.all[[g]], function(x) {mvrnorm(1,rep(0,p),x)}, simplify = FALSE)
    # likelihoods
    DNORMs.all[[g]] <- sapply(1:maxclustern, function(nr) dmvnorm(Y, mean = MUs.all[[g]][[nr]], sigma =SIGs.all[[g]][[nr]]))
    alpha.all[[g]] <- rep(2*n,iternum+1)
    
    
  }
  
  for (m in 1:iternum){  
    for (g in 1:G){
      
      Y <- Y.all[[g]]
      n <- ns[g]
      p <- ps[g]
      z <- z.all[[g]]
      N <- N.all[[g]]
      SIGs <- SIGs.all[[g]]
      MUs <- MUs.all[[g]]
      DNORMs <- DNORMs.all[[g]]
      alpha <- alpha.all[[g]]
      if(is.null(ex.data[[g]])) {  H<-H1  } else {SIM <-ex.data[[g]]} 
    
      for (i in 1:n){
        
        if (N[z[i]]>1 && length(N[N!=0])<maxclustern) {   #### Find non-singelton obserations and try to open a new cluster
          
          newsig <- genPositiveDefMat(p)$Sigma ;  newmu <- mvrnorm(1,rep(0,p),newsig);
          a1=dmvnorm(Y[i,], mean = newmu, sigma =newsig);
          a0<-DNORMs[i,z[i]];
          
          if (a0==0 || runif(1)<(alpha[m]/(n-1))*(a1/a0)){
            N[z[i]]<-N[z[i]]-1; z[i]<-which(N==0)[1]; N[z[i]]<-1;  SIGs[[z[i]]]<-newsig; MUs[[z[i]]]<-newmu;
          } #end if 
          
        } #end if (N[z[i]]>1)
        
        else if (N[z[i]]==1) { ### if i-th observarion is a singelton send it to an alreay occupied cluster with probability prop to N[i]/n-1.
          BB<-N; BB[z[i]]<-0; index<-which(BB!=0);
          b0<- DNORMs[i,z[i]];
          prob<-c(H(i,index)/sum(H(i,index)))
          o<-0
          repeat{
            o=o+1
            v<-sample(1:length(index),size=1,prob=prob) 
            sss<-prob[v] 
            z1<-index[v]
            if(DNORMs[i,z1]!=0 || sss==1 || o>100) {
              break
            }
          }#end repeat
          b1<- DNORMs[i,z1]
          if (b0==0 || runif(1)<((n-1)/alpha[m])*b1/b0) {
            SIGs[[z[i]]]=MUs[[z[i]]]=0
            N[z[i]]<-0
            z[i]<-z1
            N[z1]<-N[z1]+1 
          }# end if
          
          
        }# end else if (N[z[i]]==1)      
      }# end for i
      
      
      ######################################################################
      # Ignore singelton obs, exchange obs among those clusters with size>1#
      ######################################################################
      
      for (i in 1:n){ 
        ind2<-which(N>1)
        if (length(ind2)>1){
          
          d=H(i,ind2)*DNORMs[i,ind2] 
          if (sum(d)!=0){
            o<-0
            repeat{
              o=o+1
              w<-d/sum(d)
              newz<-sample(ind2,size=1, prob=w)
              if (DNORMs[i,newz]!=0 || w[which(ind2==newz)]==1 || o>100) 
              {break}
            }#end repeat
            #c1<-H(i,newz)*dnorm(y[i],beta[newz,]%*%X[i,],sqrt(sigma[newz]))
            #c1<-dnorm(y[i],beta[newz,]%*%X[i,],sqrt(sigma[newz]))
            
            #if (c0==0 || runif(1)<(c1)/(c0)){
            N[z[i]]<-N[z[i]]-1
            z[i]<-newz
            N[newz]<-N[newz]+1
            #} #end if
          }#end if
        }#end if
      }# end for i
      z.all[[g]] <- z
      N.all[[g]] <- N
      SIGs.all[[g]] <-  SIGs
      MUs.all[[g]] <-   MUs 
      DNORMs.all[[g]] <-  DNORMs
      alpha.all[[g]][m+1]<-randalpha(aa,bb,n,length(N[N>0]),alpha0=alpha[m],itr = 5)  
      
    }#end g
    
    # remove & merge clusters that are less populated
    for (g in 1:G){
      for (i in 1:ns[g]){
        if (N.all[[g]][z.all[[g]][i]]< 0.08*ns[g]) {
          N.all[[g]][z.all[[g]][i]]=N.all[[g]][z.all[[g]][i]]-1; SIGs.all[[g]][[z.all[[g]][i]]]=0; MUs.all[[g]][[z.all[[g]][i]]]=0; DNORMs.all[[g]][i,z.all[[g]][i]]=0
          z.all[[g]][i] <- which(N.all[[g]]> 0.08*ns[g])[which.max(DNORMs.all[[g]][i,which(N.all[[g]]> 0.08*ns[g])])];  N.all[[g]][z.all[[g]][i]]= N.all[[g]][z.all[[g]][i]]+1
        }
      }
    }
    ###############################################################
    # Parameter sampling from posteriors
    ###############################################################
    Y.list= S.list = Sig.list = C.list= Mu.list <- list()
    Ns <- indicator <- numeric()
    c=0
    for (g in 1:G){
      for (k in which(N.all[[g]]>0)) {
        indicator=append(indicator,g) # to indicate the stages
        c=c+1
        Y.list[[c]]=matrix(Y.all[[g]][z.all[[g]]==k,], ncol=p);
        Ns[c]=N.all[[g]][k]
        Sig.list[[c]]= SIGs.all[[g]][[k]]
        Mu.list[[c]] <- MUs.all[[g]][[k]] <- mvrnorm(1,colMeans(Y.list[[c]])*Ns[c]/(k0+Ns[c]),SIGs.all[[g]][[k]]/(k0+Ns[c]))
        S.list[[c]] <- t(Y.list[[c]]-Mu.list[[c]])%*%(Y.list[[c]]-Mu.list[[c]])
        C.list[[c]] <-solve(Sig.list[[c]])
      } #end for k
    } # end g
    
    #S.hat=BayesFGLasso_Columnwise(Y.list, S.list, Ns, replicate(length(Ns),list(diag(p))), replicate(length(Ns),list(diag(p)))
    #,Mu.list, indicator=indicator, a_lambda=a_lambda, b_lambda=b_lambda, burnin=burnin, nmc=nmc)$S
    # following calculates sparseC at the last iteration
    indicator=matrix(1, length(C.list),length(C.list))
    A=BayesFGLasso_ColumnwiseNew(Y.list, S.list, Ns, replicate(length(Ns),list(diag(p))), replicate(length(Ns),list(diag(p)))
                              ,Mu.list, ind.matrix = indicator, a_lambda=a_lambda, b_lambda=b_lambda, burnin=burnin, nmc=nmc)
    S.hat=A$S
    if (m==iternum){
      sparseCs=Cs=replicate(length(S.hat), list())
      for (k in 1:length(S.hat)){
        Cs[[k]]= apply(A$post.C[[k]],1:2, mean); # precision matrix
        quant_low=apply(A$post.C[[k]], 1:2, function(x) quantile(x, 0.05))
        quant_up=apply(A$post.C[[k]], 1:2, function(x) quantile(x, 0.95))
        #sparsify precision matrix using credible intervals
        sparseCs[[k]] <- Cs[[k]]
        sparseCs[[k]][mapply(function(x,y) {findInterval(0,c(x,y))==1 }, x=quant_low, y=quant_up, SIMPLIFY = TRUE)]=0

       }
    }
     c=0
    for (g in 1:G){
      for (k in which(N.all[[g]]>0)) {
        c=c+1
        SIGs.all[[g]][[k]] <- S.hat[[c]]
        DNORMs.all[[g]][,k] <-  dmvnorm(Y.all[[g]], mean = MUs.all[[g]][[k]], sigma=SIGs.all[[g]][[k]])
      } #end for k
    } # end g
    #print(paste("iteration nr= ",m))
    
  } # end for m
  
  
  
  
  Out<-list(groups=N.all,z.all=z.all, S.all=S.hat, Mu.list=Mu.list, sparseCs=sparseCs, Cs=Cs)
  return(Out)
  
  print("##################################################################################################")
}


######################################################
## Bayesian fused graphical lasso estimation function#
######################################################


BayesFGLasso_Columnwise <- function(Y.list, S.list, Ns, Sig.list, C.list, Mu.list, indicator=c(1,2,2,3,3,3),a_lambda=1, b_lambda=0.2, burnin=20, nmc=100){
  
  K=length(Ns);
  p=dim(Y.list[[1]])[2];
  indmx <- matrix(1:p^2 , ncol=p)
  upperind <- indmx[upper.tri(indmx,diag=FALSE)]
  lowerind <- indmx[lower.tri(indmx,diag=FALSE)]
  #apost <- a_lambda + p*(p+1)/2;
  # initiate lambda and tau for all groups 
  lambda1 <- rep(1,K)
  tau.list <- replicate(K,matrix(1,p,p), simplify = FALSE);
  add <- function(x) Reduce("+", x)
  ind_noi_all <- matrix(0,p-1,p); 
  for (i in 1:p){
    if (i==1){ind_noi = 2:p} else if (i==p) {ind_noi = 1:(p-1)} else  ind_noi = c(1:(i-1),(i+1):p);i
    ind_noi_all[,i] = ind_noi;
  }
  post.Sig <- post.C <- replicate(K, array(0,c(p,p,nmc)), simplify = FALSE)
  
  post.lmd <- matrix(0,K,nmc);
  
  for (iter in 1: (burnin+nmc)){ 
    for (k in 1:K){ 
      # Sample lambda 
      bpost1 <- b_lambda + sum(abs(C.list[[k]]))/2;    
      lambda1[k] <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost1)
      # sample tau off-diagonal        
      Cadjust = pmax(abs(C.list[[k]][upperind]),10^-6);        
      lambda_prime = lambda1[k]^2;  mu_prime = pmin(lambda1[k]/Cadjust,10^12);
      tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
      tau.list[[k]][upperind] <- tau_temp;
      tau.list[[k]][lowerind] <- tau_temp;
      if(indicator[k]==1) {
        k_sibling  <- sample(size=1, which(indicator==2));
        bpost <- b_lambda + sum(abs(C.list[[k]]-C.list[[k_sibling]]))/2;
        lambdaf <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost)
        Cadjust = pmax(abs(C.list[[k_sibling]][upperind]),10^-6);        
        lambda_prime = lambdaf^2;  mu_prime = pmin(lambdaf/Cadjust,10^12);
        tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
        tau.list[[k_sibling]][upperind] <- tau_temp;
        tau.list[[k_sibling]][lowerind] <- tau_temp;
        A <- 1/tau.list[[k]] +1/tau.list[[k_sibling]]
        B <- C.list[[k_sibling]]/tau.list[[k_sibling]]
      } else if(indicator[k]==max(indicator)) {
        k_parent  <- sample(size=1,which(indicator==max(indicator)-1));
        bpost <- b_lambda + sum(abs(C.list[[k]]-C.list[[k_parent]]))/2;
        lambdaf <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost);
        Cadjust = pmax(abs(C.list[[k_parent]][upperind]),10^-6);        
        lambda_prime = lambdaf^2;  mu_prime = pmin(lambdaf/Cadjust,10^12);
        tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
        tau.list[[k_parent]][upperind] <- tau_temp;
        tau.list[[k_parent]][lowerind] <- tau_temp;
        A <- 1/tau.list[[k]] +1/tau.list[[k_parent]]
        B <- C.list[[k_parent]]/tau.list[[k_parent]]
      }else {
        k_sibling  <- sample(size=1,which(indicator==indicator[k]+1));
        k_parent  <- sample(size=1,which(indicator==indicator[k]-1));
        bpost <- b_lambda + sum(abs(C.list[[k]]-C.list[[k_parent]]))/2;
        lambdaf1 <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost);
        bpost <- b_lambda + sum(abs(C.list[[k]]-C.list[[k_sibling]]))/2;
        lambdaf2 <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost);
        Cadjust = pmax(abs(C.list[[k_sibling]][upperind]),10^-6);        
        lambda_prime = lambdaf2;  mu_prime = pmin(lambdaf2/Cadjust,10^12);
        tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
        tau.list[[k_sibling]][upperind] <- tau_temp;
        tau.list[[k_sibling]][lowerind] <- tau_temp;
        
        Cadjust = pmax(abs(C.list[[k_parent]][upperind]),10^-6);        
        lambda_prime = lambdaf1^2;  mu_prime = pmin(lambdaf1/Cadjust,10^12);
        tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
        tau.list[[k_parent]][upperind] <- tau_temp;
        tau.list[[k_parent]][lowerind] <- tau_temp;
        A <- 1/tau.list[[k]] + 1/tau.list[[k_sibling]] + 1/tau.list[[k_parent]]
        B <- add(mapply('/',  list(C.list[[k_parent]],C.list[[k_sibling]]),
                        list(tau.list[[k_parent]],tau.list[[k_sibling]]), SIMPLIFY = FALSE))
        
      }
      
      
      # sample Sig and C
      for (i in 1:p){
        ind_noi <- ind_noi_all[,i];
        A_temp <- A[ind_noi,i];
        B_temp <- B[ind_noi,i];
        #tau_temp <- tau[ind_noi,i];
        # added :lambda1, A_temp, B_temp, nk
        Sig11 <- Sig.list[[k]][ind_noi,ind_noi]; Sig12 <- Sig.list[[k]][ind_noi,i];
        invC11 <- Sig11 - Sig12%*%t(Sig12)/Sig.list[[k]][i,i];
        
        Ci <- (S.list[[k]][i,i]+lambda1[k])*invC11+diag(A_temp);
        Ci[mapply(is.infinite, Ci)] <- 100000000
        Ci_chol = chol(Ci);
        mu_i = -solve(Ci,S.list[[k]][ind_noi,i]+B_temp);
        beta = mu_i+ solve(Ci_chol,rnorm(p-1));
        C.list[[k]][ind_noi,i] = beta;
        C.list[[k]][i,ind_noi] = beta;
        gam = rgamma(1,Ns[k]/2+1,rate=1/solve((S.list[[k]][i,i]+lambda1[k]),2)) 
        
        C.list[[k]][i,i] = gam+t(beta)%*%invC11%*%beta;
        
        
        # Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11%*%beta;
        
        Sig.list[[k]][ind_noi,ind_noi] = invC11+invC11beta%*%t(invC11beta)/gam;
        Sig12 = -invC11beta/gam;
        Sig.list[[k]][ind_noi,i] = Sig12;
        Sig.list[[k]][i,ind_noi] = t(Sig12);
        Sig.list[[k]][i,i] = 1/gam;
        
      }
      
      
      if (iter >burnin){           
        post.Sig[[k]][,,iter-burnin] <- Sig.list[[k]] 
        post.C[[k]][,,iter-burnin] <- C.list[[k]]
        post.lmd[k,iter-burnin] <- lambda1[k]
      };
    }# end k
  }#end iter
  S=list()
  for (k in 1:K){ S[[k]]= apply(post.Sig[[k]],1:2, mean);}
  return(list(post.Sig=post.Sig,post.C=post.C,post.lmd=post.lmd, S=S))
}


BayesFGLasso_ColumnwiseNew <- function(Y.list, S.list, Ns, Sig.list, C.list, Mu.list, ind.matrix=matrix(1,length(Ns),length(Ns)),a_lambda=1, b_lambda=0.2, burnin=20, nmc=100){
  # not completed!
  # not applied: another idea is to choose parents and siblings for a given network based on probability
  #proportional to the l_1 distances between every possible parent/sibling networks.
  K=length(Ns);
  p=dim(Y.list[[1]])[2];
  indmx <- matrix(1:p^2 , ncol=p)
  upperind <- indmx[upper.tri(indmx,diag=FALSE)]
  lowerind <- indmx[lower.tri(indmx,diag=FALSE)]
  #apost <- a_lambda + p*(p+1)/2;
  # initiate lambda and tau for all groups 
  lambda1 <- rep(1,K)
  tau.list <- replicate(K,matrix(1,p,p), simplify = FALSE);
  add <- function(x) Reduce("+", x)
  ind_noi_all <- matrix(0,p-1,p); 
  for (i in 1:p){
    if (i==1){ind_noi = 2:p} else if (i==p) {ind_noi = 1:(p-1)} else  ind_noi = c(1:(i-1),(i+1):p);i
    ind_noi_all[,i] = ind_noi;
  }
  post.Sig <- post.C <- replicate(K, array(0,c(p,p,nmc)), simplify = FALSE)
  
  post.lmd <- matrix(0,K,nmc);
  
  for (iter in 1: (burnin+nmc)){ 
    for (k in 1:K){ 
      # Sample lambda 
      bpost1 <- b_lambda + sum(abs(C.list[[k]]))/2;    
      lambda1[k] <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost1)
      # sample tau off-diagonal        
      Cadjust = pmax(abs(C.list[[k]][upperind]),10^-6);        
      lambda_prime = lambda1[k]^2;  mu_prime = pmin(lambda1[k]/Cadjust,10^12);
      tau_temp <-  1/rinvgauss(length(mu_prime),mu_prime,lambda_prime);
      tau.list[[k]][upperind] <- tau_temp;
      tau.list[[k]][lowerind] <- tau_temp;
      
      k_prime  <- which(ind.matrix[k,]==1);
      if(length(k_prime)>1){
        
        b_post <- b_lambda + sapply(C.list[k_prime], function (x,y) sum(abs(x-y)/2), x=C.list[[k]])
        lambdaf=sapply(b_post, function(x) rgamma(1,a_lambda + p*(p+1)/2,rate=x))
        Cadjusts = sapply(C.list[k_prime], function(x) pmax(abs(x[upperind]),10^-6), simplify = FALSE);
        lambda_prime = lambdaf^2;  mu_prime = mapply(function(x,y) pmin(x/y,10^12) ,x=lambdaf,y=Cadjusts, SIMPLIFY = FALSE)  
        tau_temp <- mapply(function(x,y) rinvgauss(length(x),x,y), x=mu_prime,y=lambda_prime, SIMPLIFY = FALSE);
        c=0; for (kk in k_prime) { c=c+1; tau.list[[kk]][upperind]= tau.list[[kk]][lowerind]=tau_temp[[c]]}
        A <- add(mapply("/",1,tau.list[c(k,k_prime)],SIMPLIFY = FALSE))
        B <- add(mapply('/',  C.list[k_prime],tau.list[k_prime], SIMPLIFY = FALSE))
      } 
      if(length(k_prime)==1){
        bpost <- b_lambda + sum(abs(C.list[[k]]-C.list[[k_prime]]))/2;
        lambdaf <- rgamma(1,a_lambda + p*(p+1)/2,rate=bpost);
        Cadjust = pmax(abs(C.list[[k_prime]][upperind]),10^-6);        
        lambda_prime = lambdaf^2;  mu_prime = pmin(lambdaf/Cadjust,10^12);
        tau_temp <-  rinvgauss(length(mu_prime),mu_prime,lambda_prime);
        tau.list[[k_prime]][upperind] <- tau.list[[k_prime]][lowerind] <- tau_temp;
        A <- 1/tau.list[[k]] +1/tau.list[[k_prime]]
        B <- C.list[[k_prime]]/tau.list[[k_prime]]}
      if(length(k_prime)==0) {
        A <- B <- diag(p)
        
      }
      
      
      # sample Sig and C
      for (i in 1:p){
        ind_noi <- ind_noi_all[,i];
        A_temp <- A[ind_noi,i];
        B_temp <- B[ind_noi,i];
        #tau_temp <- tau[ind_noi,i];
        # added :lambda1, A_temp, B_temp, nk
        Sig11 <- Sig.list[[k]][ind_noi,ind_noi]; Sig12 <- Sig.list[[k]][ind_noi,i];
        invC11 <- Sig11 - Sig12%*%t(Sig12)/Sig.list[[k]][i,i];
        
        Ci <- (S.list[[k]][i,i]+lambda1[k])*invC11+diag(A_temp);
        Ci[mapply(is.infinite, Ci)] <- 100000000
        Ci_chol = chol(Ci);
        mu_i = -solve(Ci,S.list[[k]][ind_noi,i]+B_temp);
        beta = mu_i+ solve(Ci_chol,rnorm(p-1));
        C.list[[k]][ind_noi,i] = beta;
        C.list[[k]][i,ind_noi] = beta;
        gam = rgamma(1,Ns[k]/2+1,rate=1/solve((S.list[[k]][i,i]+lambda1[k]),2)) 
        
        C.list[[k]][i,i] = gam+t(beta)%*%invC11%*%beta;
        
        
        # Below updating Covariance matrix according to one-column change of precision matrix
        invC11beta = invC11%*%beta;
        
        Sig.list[[k]][ind_noi,ind_noi] = invC11+invC11beta%*%t(invC11beta)/gam;
        Sig12 = -invC11beta/gam;
        Sig.list[[k]][ind_noi,i] = Sig12;
        Sig.list[[k]][i,ind_noi] = t(Sig12);
        Sig.list[[k]][i,i] = 1/gam;
        
      }
      
      
      if (iter >burnin){           
        post.Sig[[k]][,,iter-burnin] <- Sig.list[[k]] 
        post.C[[k]][,,iter-burnin] <- C.list[[k]]
        post.lmd[k,iter-burnin] <- lambda1[k]
      };
    }# end k
  }#end iter
  S=list()
  for (k in 1:K){ S[[k]]= apply(post.Sig[[k]],1:2, mean);}
  return(list(post.Sig=post.Sig,post.C=post.C,post.lmd=post.lmd, S=S))
}

