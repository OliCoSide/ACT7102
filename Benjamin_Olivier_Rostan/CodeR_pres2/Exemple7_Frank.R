###Exemple 7 de Cossette et al. (2011)

densM1 <- function(xx) dpois(xx, lambda = 4)
densM2 <- function(xx) dnbinom(xx, size = 4, prob = 0.5)

repM1 <- function(xx) ppois(xx, lambda = 4)
repM2 <- function(xx) pnbinom(xx, size = 4, prob = 0.5)

Frank <- function(u1,u2, alpha){
  if (alpha ==0)
   u1*u2
  else
  -1/alpha * log(1+(exp(-alpha*u1)-1)*(exp(-alpha*u2)-1)/(exp(-alpha)-1))
}

repM1M2 <- function(xx1, xx2, alpha) Frank(repM1(xx1),repM2(xx2),alpha)

densM1M2 <- function(xx1,xx2,alpha) repM1M2(xx1,xx2,alpha)-repM1M2(xx1-1,xx2,alpha)-repM1M2(xx1,xx2-1,alpha)+repM1M2(xx1-1,xx2-1,alpha)

aB1 <- 0.5
bB1 <- 0.1
aB2 <- 0.25
bB2 <- 0.1

i <- 0:200

EM1 <- sum(i*densM1(i)) #Espérance de M1
EM2 <- sum(i*densM2(i)) #Espérance de M2

DMM1 <- sum(i^2*densM1(i)) #Deuxième moment de M1
DMM2 <- sum(i^2*densM2(i)) #Deuxième moment de M2

VarM1 <- DMM1-EM1^2 #Variance de M1
VarM2 <- DMM2-EM2^2 #Variance de M2

EX1 <- EM1*aB1/bB1   #Espérance de X1
EX2 <- EM2*aB2/bB2   #Espérance de X2
EXS <- EX1+EX2

VarX1 <- EM1*aB1/bB1^2+aB1^2/bB1^2*VarM1 #Variance de X1
VarX2 <- EM2*aB2/bB2^2+aB2^2/bB2^2*VarM2 #Variance de X2

table1 <- data.frame(alpha = NULL, CovM1M2 = NULL, corrM1M2 = NULL, CovX1X2 = NULL, VarS = NULL)

for (alph in c(-20,0,20)){
  i <- 0:200
  EM1M2 <- i%*%outer(i,i,densM1M2, alpha = alph)%*%i
  CovM1M2 <- round(EM1M2 - EM1*EM2,4)
  corrM1M2 <- round(CovM1M2/(sqrt(VarM1*VarM2)),4)
  CovX1X2 <- round(CovM1M2*aB1/bB1*aB2/bB2,4)
  VarS <- round(VarX1+VarX2+2*CovX1X2,4)
  ttt <- data.frame(data.frame(alpha = alph, CovM1M2 = CovM1M2, corrM1M2 = corrM1M2, CovX1X2 = CovX1X2, VarS = VarS))
  table1 <- rbind(table1, ttt)
}
table1

repX1 <- function(x) densM1(0)+sum(densM1(1:500)*pgamma(x,shape = (1:500)*0.5,rate = 0.1))
repX2 <- function(x) densM2(0)+sum(densM2(1:500)*pgamma(x,shape = (1:500)*0.25,rate = 0.1))

VaRkX1 <- function(k) optimize(function(x,k) abs(repX1(x)-k), c(0,500), k=k)$minimum
VaRkX2 <- function(k) optimize(function(x,k) abs(repX2(x)-k), c(0,500), k=k)$minimum

TVaRkX1 <- function(k) 1/(1-k)*sum(densM1(0:500)*(0:500)*0.5/0.1*(1-pgamma(VaRkX1(k),(0:500)*0.5+1,0.1)))
TVaRkX2 <- function(k) 1/(1-k)*sum(densM2(0:500)*(0:500)*0.25/0.1*(1-pgamma(VaRkX2(k),(0:500)*0.25+1,0.1)))

                                   
table2 <- data.frame(kappa = NULL, VaRX1 = NULL, VaRX2 = NULL, TVaRX1 = NULL, TVaRX2 = NULL, sumVaR = NULL, sumTVaR = NULL)

for (k in c(0.25,0.50,0.95,0.99,0.995)){
  VaRX1 <- round(VaRkX1(k),4)
  VaRX2 <- round(VaRkX2(k),4)
  TVaRX1 <- round(TVaRkX1(k),4)
  TVaRX2 <- round(TVaRkX2(k),4)
  sumVaR <- round(VaRX1 + VaRX2,4)
  sumTVaR <- round(TVaRX1 + TVaRX2,4)
  ttt <- data.frame(kappa = k, VaRX1 = VaRX1, VaRX2 = VaRX2, TVaRX1 = TVaRX1, TVaRX2 = TVaRX2, sumVaR = sumVaR, sumTVaR = sumTVaR)
  table2 <- rbind(table2, ttt)
}
table2

dens <- function(alpha) outer(0:100, 0:100, densM1M2, alpha = alpha)
sum.des.ma <- outer((0:100)*0.5,(0:100)*0.25,"+")
des.ma1 <- matrix((0:100)*0.5, nrow=101,ncol = 101)
des.ma2 <- matrix((0:100)*0.25, byrow = TRUE, nrow = 101, ncol = 101)

repS <- function(x, alpha) sum(dens(alpha)*pgamma(x, sum.des.ma, rate = 0.1))

VaRkS <- function(k, alpha) optimize(function(x,k,alpha) abs(repS(x, alpha)-k), c(5,200), k=k, alpha=alpha)$minimum

TVaRkS <- function(k,alpha) 1/(1-k)*sum(dens(alpha)*sum.des.ma/0.1*(1-pgamma(VaRS,sum.des.ma+1,rate = 0.1)))
TVaRkX1.S <- function(k,alpha) 1/(1-k)*sum(dens(alpha)*des.ma1/0.1*(1-pgamma(VaRS,sum.des.ma+1,rate = 0.1)))
TVaRkX2.S <- function(k,alpha) 1/(1-k)*sum(dens(alpha)*des.ma2/0.1*(1-pgamma(VaRS,sum.des.ma+1,rate = 0.1)))

table3 <- data.frame(alpha = NULL, kappa = NULL, VaRS = NULL, TVaRS = NULL, TVaRX1.S = NULL, TVaRX2.S = NULL)

for (alph in c(-20,0,20)){
for (k in c(0.25,0.50,0.95,0.99,0.995)){
  VaRS <- VaRkS(k, alph)
  TVaRS <- TVaRkS(k, alph)
  TVaRX1.S <- TVaRkX1.S(k, alph)
  TVaRX2.S <- TVaRkX2.S(k, alph)
  ttt <- data.frame(alpha = alph, kappa = k, VaRS = VaRS, TVaRS = TVaRS, TVaRX1.S = TVaRX1.S, TVaRX2.S = TVaRX2.S)
  table3 <- rbind(table3, ttt)
}
}
##Temps de calcul est ~2min. 

table3
