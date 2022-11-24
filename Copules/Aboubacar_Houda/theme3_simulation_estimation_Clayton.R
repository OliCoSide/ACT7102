
## densité de la copule Clayton
Clayc<-function(u1, u2, alph){
  ( (1+alph)/((u1*u2)^(alph+1)) )*(u1^(-alph)+u2^(-alph)-1)^(-2 -1/alph)
}

## densité de la copule Frank
Frankc <- function(u1, u2, alph){
  (-1/alph)*log(1+( (exp(-alph*u1)-1)*(exp(-alph*u2)-1)/(exp(-alph)-1) ))
}

##################################  le cas bien spécifié (marginales suivent des normales)
############## Méthode de MLE 
## fonction qui calcule le log vraisemblance afin d'estimer apres les parametres alpha et theta
log_like <- function(param, X1, X2,n){
  theta <- param[1]
  mean1 <- param[2]
  sd1 <- param[3]
  mean2 <- param[4]
  sd2 <- param[5]
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  loglik  <- sum(sapply(1:n, function(i)
    log(Clayc(pnorm(X1[i],mean=mean1,sd=sd1),pnorm(X2[i],mean=mean2,sd=sd2),theta))
    +log(dnorm(X1[i],mean=mean1,sd=sd1)) +log(dnorm(X2[i],mean=mean2,sd=sd2)) ))
  return(-loglik)
}

###fonction qui simule les distributions X1 et X2 et estimer les alpha et le theta de la copule de clayton
simulMLE1 <- function(alph, nsim, N){
  
  theta_MLcc <- rep(0,N)
  set.seed(2022)
  for (i in 1:N){
  vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
  vTheta<-qgamma(vV[,1], 1/alph, 1)
  vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
  vU<-(1+vY)**(-1/alph)
  #### X1 et X2 suivent une normale
  vX<-sapply(1:2, function(t) qnorm(vU[,t], mean=0,sd=1))
  colnames(vX) <- c("X1","X2")
  #### estimation des parametres alpha et theta
  MLE_estimates <- optim(fn=log_like,
                         par=c(0.001,0,1,0,1),
                         X1 = vX[,1],
                         X2 = vX[,2],
                         n=100)
  theta_MLcc[i] <- MLE_estimates$par[1]
  }
  MSE <- mean((theta_MLcc - alph)^2)
  return(MSE)
}

simulMLE1(0.5,100,250)#0.03266325
simulMLE1(2,100,250)#0.1470203
simulMLE1(8,100,250)#1.533016

########## Méthode de IFM
########### déclaration de la fonction pour ML afin d'estimer les paramètres alpha1 et alpha2
X1.lik <- function(alpha1, X){
  mean1 <- alpha1[1]
  sd1 <- alpha1[2]
  X1 <- as.matrix(X[,1])
  logl <-  sum(sapply(1:length(X1), function(i) log(dnorm(X1[i],mean= mean1, sd= sd1)) ))
  return(-logl)
}


X2.lik <- function(alpha2, X){
  mean2 <- alpha2[1]
  sd2 <- alpha2[2]
  X2 <- as.matrix(X[,2])
  logl <-  sum(sapply(1:length(X2), function(i) log(dnorm(X2[i],mean= mean2, sd= sd2)) ))
  return(-logl)
}

####fonction de log vraisemblance pour estimer de theta de la copule
theta.lik<- function(theta, X){
  mean1 <- X1.estim$par[1]
  sd1 <- X1.estim$par[2]
  mean2 <- X2.estim$par[1]
  sd2 <- X2.estim$par[2]
  X1 <- as.matrix(X[,1])
  X2 <- as.matrix(X[,2])
  
  #F1 <- pnorm(X1,mean=mean1,sd=sd1)
  #F2 <- pnorm(X2,mean=mean2,sd=sd2)
  #logl <-  sum(log(Clayc(F1,F2,theta) ) )
  
  logl <-  sum(sapply(1:100, function(i)
    log(Clayc(pnorm(X1[i],mean=mean1,sd=sd1),pnorm(X2[i],mean=mean2,sd=sd2),theta) ) ))
  return(-logl)
}

########## fonction qui simule les X1 et X2 et estime le paramètre theta de la copule
simul_IFM1 <- function(alph, nsim, N, ic){
  theta_IFMcc <- rep(0,N)
  set.seed(2022)
  
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent une normale
    vX<-sapply(1:2, function(t) qnorm(vU[,t], mean=0,sd=1))
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    ### estimation des paramètres alpha1 et alpha2
    X1.estim <- optim(par= c(0.01,0.5), X1.lik, X=vX, method="BFGS")
    X2.estim <- optim(par= c(0.01,0.5), X2.lik, X=vX, method="BFGS")
    
    theta.estim <- optimize(theta.lik, ic, X=vX)
    theta_IFMcc[i] <- theta.estim$minimum
    
  }
  return(mean((theta_IFMcc - alph)^2))
}

simul_IFM1(0.5,100,250,c(0,4))#0.03165538
simul_IFM1(2,100,250,c(0,6))#0.141516
simul_IFM1(8,100,250,c(0,12))#1.314598

##### Méthode PML

########## definition de la fonction ML pour estimer le paramètre theta en utilisant la méthode SP
theta.liksp<- function(theta, X){
  #X1 <- as.matrix(X[,1])
  #X2 <- as.matrix(X[,2])
  
  F1 <- ecdf(X[,1])
  F2 <- ecdf(X[,2])
  logl <-  sum(log(Clayc(F1(X[,1]),F2(X[,2]),theta) ) )
  
  #logl <-  sum(sapply(1:100, function(i) log(Clayc(knots(ecdf(X1))[i],knots(ecdf(X2))[i],theta) ) ))
  return(-logl)
}

##### fonction qui estime le theta de la copule Clayton par la méthode SP
simul_PML1 <- function(alph, nsim, N){
  theta_SPcc <- rep(0,N)
  set.seed(2022)
  
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent une normale
    vX<-sapply(1:2, function(t) qnorm(vU[,t], mean=0,sd=1))
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    theta.estim <- optim(par= 0.01, theta.liksp, X=vX)
    theta_SPcc[i] <-theta.estim$par
  }
  return(mean((theta_SPcc - alph)^2))
}

simul_PML1(0.5,100,250)#0.0450031
simul_PML1(2,100,250)#0.1717364
simul_PML1(8,100,250)#1.511804

##### Méthode Bench-Mark

######## fonction de ML 
theta.lik<- function(theta, X){
  X1 <- as.matrix(X[,1])
  X2 <- as.matrix(X[,2])
  
  #F1 <- pnorm(X1,mean=mean1,sd=sd1)
  #F2 <- pnorm(X2,mean=mean2,sd=sd2)
  #logl <-  sum(log(Clayc(F1,F2,theta) ) )
  
  logl <-  sum(sapply(1:100, function(i)
    log(Clayc(pnorm(X1[i],mean=0,sd=1),pnorm(X2[i],mean=0,sd=1),theta) ) ))
  return(-logl)
}

###### fonction qui simulent les X1 et X2 afin d'estimer le theta de la copule en utilisant le BM
simul_BM1 <- function(alph, nsim, N, ic){
  theta_BMcc <- rep(0,N)
  set.seed(2022)
  
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent une normale
    vX<-sapply(1:2, function(t) qnorm(vU[,t], mean=0,sd=1))
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    theta.estim <- optimize(theta.lik, X=vX, ic)
    theta_BMcc[i] <- theta.estim$minimum
  }
  return(mean((theta_BMcc - alph)^2))
}

simul_BM1(0.5,100,250, c(0,2))#0.02541772
simul_BM1(2,100,250, c(0,5))#0.07962646
simul_BM1(8,100,250, c(0,9.2))#0.4538143

###################################### le cas mal spécifié ou X1 et X2 suivent des student de df=3
######## Méthode MLE

#### fonction qui estime les alpha et le theta de la copule de Clayton par la méthode MLE
simul_MLE2 <- function(alph, nsim, N){
  
  theta_MLcc <- rep(0,N)
  set.seed(2022)
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    vX<-sapply(1:2, function(t) qt(vU[,t], df=3))
    colnames(vX) <- c("X1","X2")
    #### estimation des parametres alpha et theta
    MLE_estimates <- optim(fn=log_like,
                           par=c(0.00001,0,1,0,1),
                           X1 = vX[,1],
                           X2 = vX[,2],
                           n=100)
    theta_MLcc[i] <- MLE_estimates$par[1]
  }
  MSE <- mean((theta_MLcc - alph)^2)
  EB <- mean((theta_MLcc - alph))
  SD <- sd(theta_MLcc)
  #fig <- hist(theta_MLcc,
  #           main="Histogramme des estimations de theta",
  #           xlab="theta estimé",
  #           ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD,fig)
  return(resultat)
}

simul_MLE2(0.5,100,250)#MSE=0.1085084 et EB=-0.1113198 et sd=0.3106482
simul_MLE2(2,100,250)#MSE=1.047656 et EB=0.4045985 et sd=0.9420753
simul_MLE2(8,100,250)#MSE= 25.64878 et EB=0.81771 et sd=1.827799

########### Méthode IFM
simul_IFM2 <- function(alph, nsim, N, ic){
  theta_IFMcc <- rep(0,N)
  set.seed(2022)
  
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent une normale
    vX<-sapply(1:2, function(t) qt(vU[,t], df=3))
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    ### estimation des paramètres alpha1 et alpha2
    X1.estim <- optim(par= c(-10,10), X1.lik, X=vX)
    X2.estim <- optim(par= c(-10,10), X2.lik, X=vX)
    
    theta.estim <- optimize(theta.lik, ic, X=vX)
    theta_IFMcc[i] <- theta.estim$minimum
    
  }
  MSE <- mean((theta_IFMcc - alph)^2)
  EB <- mean((theta_IFMcc - alph))
  SD <- sd(theta_IFMcc)
  fig <- hist(theta_IFMcc,
              main="Histogramme des estimations de theta",
              xlab="theta estimé",
              ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD,fig)
  return(resultat)
}

simul_IFM2(0.5,100,250,c(0,2))#MSE=0.06704396 et EB=-0.1487564 et sd=0.2123579
simul_IFM2(2,100,250,c(0,5))#MSE=0.3047202 et EB=-0.06101557 et sd=0.5497328
simul_IFM2(8,100,250,c(5,10))#MSE= 4.304784et EB=-0.9296558 et sd=1.858586

########## Méthode PML
########## definition de la fonction ML pour estimer le paramètre theta en utilisant la méthode SP
theta.liksp<- function(theta, X){
  #X1 <- as.matrix(X[,1])
  #X2 <- as.matrix(X[,2])
  
  F1 <- ecdf(X[,1])
  F2 <- ecdf(X[,2])
  logl <-  sum(log(Clayc(F1(X[,1]),F2(X[,2]),theta) ) )
  
  #logl <-  sum(sapply(1:100, function(i) log(Clayc(knots(ecdf(X1))[i],knots(ecdf(X2))[i],theta) ) ))
  return(-logl)
}

simul_PML2 <- function(alph, nsim, N){
  theta_SPcc_stu <- rep(0,N)
  set.seed(2022)
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent des students de df=3
    vX<-sapply(1:2, function(t) qt(vU[,t], df=3))
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    #theta.estim <- optim(par= 0.01,lower = 0,upper=5, theta.liksp, X=vX, method= "L-BFGS-B")
    theta.estim <- optim(par= 0.0001, theta.liksp, X=vX)
    theta_SPcc_stu[i] <-theta.estim$par
  }
  MSE <- mean((theta_SPcc_stu - alph)^2)
  EB <- mean((theta_SPcc_stu - alph))
  SD <- sd(theta_SPcc_stu)
  fig <- hist(theta_SPcc_stu,
              main="Histogramme des estimations de theta",
              xlab="theta estimé",
              ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD,fig)
  return(resultat)
}

simul_PML2(0.5,100,250)#MSE=0.04500363 et EB=0.05322555 et sd=0.2057669
simul_PML2(2,100,250)#MSE=0.1717355 et EB=0.05424089 et sd=0.4116689
simul_PML2(8,100,250)#MSE=1.511793 et EB=-0.3108083 et sd=1.192004

################################## le cas mal spécifié ou X1 suit une student de df=3 et X2 suit une khi-carré de df=2
######## Méthode MLE
#### fonction qui estime les alpha et le theta de la copule de Clayton par la méthode MLE
simul_MLE3 <- function(alph, nsim, N){
  
  theta_MLcc <- rep(0,N)
  set.seed(2022)
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #### X1 qui suit une student de df=3 et X2 qui suit une khi-carré de df=2
    vX <- matrix(NA,nsim,2)
    vX[,1] <- qt(vU[,1], df=3)
    vX[,2] <- qchisq(vU[,2], df=2)
    colnames(vX) <- c("X1","X2")
    #### estimation des parametres alpha et theta
    MLE_estimates <- optim(fn=log_like,
                           par=c(0.001,0,1,0,1),
                           X1 = vX[,1],
                           X2 = vX[,2],
                           n=100)
    theta_MLcc[i] <- MLE_estimates$par[1]
  }
  MSE <- mean((theta_MLcc - alph)^2)
  EB <- mean((theta_MLcc - alph))
  SD <- sd(theta_MLcc)
  #fig <- hist(theta_MLcc,
  #           main="Histogramme des estimations de theta",
  #           xlab="theta estimé",
  #           ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD)
  return(resultat)
}

simul_MLE3(0.5,100,250)#MSE=0.1873836 et EB=-0.005964708 et sd=0.4337055
simul_MLE3(2,100,250)#MSE=0.6271307 et EB=-0.2572067 et sd=0.7504854
simul_MLE3(8,100,250)#MSE=18.18497 et EB=-4.062702 et sd=1.298527

########### Méthode IFM
simul_IFM3 <- function(alph, nsim, N, ic){
  theta_IFMcc <- rep(0,N)
  set.seed(2022)
  
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 et X2 supposé quils suivent une normale
    vX <- matrix(NA,nsim,2)
    vX[,1] <- qt(vU[,1], df=3)
    vX[,2] <- qchisq(vU[,2], df=2)
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    ### estimation des paramètres alpha1 et alpha2
    X1.estim <- optim(par= c(-10,10), X1.lik, X=vX)
    X2.estim <- optim(par= c(-10,10), X2.lik, X=vX)
    
    theta.estim <- optimize(theta.lik, ic, X=vX)
    theta_IFMcc[i] <- theta.estim$minimum
    
  }
  MSE <- mean((theta_IFMcc - alph)^2)
  EB <- mean((theta_IFMcc - alph))
  SD <- sd(theta_IFMcc)
  fig <- hist(theta_IFMcc,
              main="Histogramme des estimations de theta",
              xlab="theta estimé",
              ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD,fig)
  return(resultat)
}

simul_IFM3(0.5,100,250,c(0,4))#MSE=0.11025 et EB=-0.09229241 et sd=0.3195946
simul_IFM3(2,100,250,c(0,7))#MSE=0.9748311 et EB=-0.8093383 et sd=0.5666454
simul_IFM3(8,100,250,c(0,11))#MSE=33.42069 et EB=-5.701224 et sd=0.9593836

simul_PML3(0.5,100,250)#MSE=0.04500363 et EB=0.05322555 et sd=0.2057669
simul_PML3(2,100,250)#MSE=0.1717355 et EB=0.05424089 et sd=0.4116689
simul_PML3(8,100,250)#MSE=1.511793 et EB=-0.3108083 et sd=1.192004
####### Méthode de PML
simul_PML3 <- function(alph, nsim,N){
  theta_SPcc_stukhi <- rep(0,N)
  set.seed(2022)
  for (i in 1:N){
    vV<-matrix(runif(nsim*3), nsim, 3, byrow=T)
    vTheta<-qgamma(vV[,1], 1/alph, 1)
    vY<-sapply(1:2, function(t) qexp(vV[,t+1], vTheta))
    vU<-(1+vY)**(-1/alph)
    #plot(vU)
    #### X1 qui suit une student de df=3 et X2 qui suit une khi-carré de df=2
    vX <- matrix(NA,nsim,2)
    vX[,1] <- qt(vU[,1], df=3)
    vX[,2] <- qchisq(vU[,2], df=2)
    colnames(vX) <- c("X1","X2")
    X1 <- vX[,1]
    X2 <- vX[,2]
    
    #theta.estim <- optim(par= 0.01,lower = 0,upper=5, theta.liksp, X=vX, method= "L-BFGS-B")
    theta.estim <- optim(par= 0.0001, theta.liksp, X=vX)
    theta_SPcc_stukhi[i] <-theta.estim$par
  }
  MSE <- mean((theta_SPcc_stukhi - alph)^2)
  EB <- mean((theta_SPcc_stukhi - alph))
  SD <- sd(theta_SPcc_stukhi)
  fig <- hist(theta_SPcc_stukhi,
              main="Histogramme des estimations de theta",
              xlab="theta estimé",
              ylab="Fréquence par rapport au nombre d'échantillons")
  resultat <- list(MSE,EB,SD, fig)
  return(resultat)
}

simul_PML3(0.5,100,250)#MSE=0.04500363 et EB=0.05322555 et sd=0.2057669
simul_PML3(2,100,250)#MSE=0.1717355 et EB=0.05424089 et sd=0.4116689
simul_PML3(8,100,250)#MSE=1.511793 et EB=-0.3108083 et sd=1.192004
