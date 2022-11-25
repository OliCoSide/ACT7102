### ARticle Ben Cossette Al: 2002 -------

###EXEMPLE 3 ###


##COOK-JONHSON

##Définissons d'abord la copule de Cook-Johnson (a = param.dép., u = vecteur)

cook.johnson <- function(a, u) (sum(u^(-a))-19)^(-1/a)

##La fonction de répartition conjointe est donc
rep.conj <- function(a, x) cook.johnson(a, pbinom(x, size=1, prob=0.05))

##La fonction de densité conjointe est donc (on emploie équation 16+1 puisque les I sont i.d.) (a = param.dép., x = vecteur de 0 et de 1)
dens.conj <- function(a, x) {
  j <- sum(x)
  
  ToSum <- numeric(j+1)
  for (k in (0:j)){
    ToSum[k+1] <- (-1)^k*factorial(j)/(factorial(k)*factorial(j-k))*rep.conj(a,c(rep(1,j-k),rep(0,20+k-j)))
  }
  sum(ToSum)
}

##La fonction génératrice de M ou M = somme(I) détermine le nombre de fois qu'il faut sommer (puisque le B_i sont iid)

pgf.M <- function(t, dens){
  sum(dens[1:21]*t^(0:20))
}

##Évaluation de la somme composée à l'aide de FFT (on a fgp de M et fdm de B)

ffst <- function(densM){
  
  #Discrétisation de B
  nfft <- 2^15
  
  fB <- c(densiteB, rep(0, nfft-length(densiteB)))
  ftB <- fft(fB)
  ftS <- sapply(ftB, pgf.M, dens=densM)
  fS <- Re(fft(ftS,inverse = TRUE))/nfft
  fS
}

i <- (0:5000)/100 
densite.upperB <- sapply(i, function(x) pgamma(x+0.01, shape = 1, rate = 0.5)-pgamma(x, shape = 1, rate = 0.5))
densite.lowerB <- c(0,densite.upperB)

i <- (0:(2^15-1))/100

##alpha = 1
##upper
densiteB <- densite.upperB
fMa1u <- sapply(0:20, dens.conj, a=1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa1u <- ffst(fMa1u)
sum(fSa1u)          #bel et bien une densité
Ea1u <- sum(fSa1u*i)      #espérance                     
E2a1u <- sum(fSa1u*i^2)    #deuxième moment
Vara1u <- E2a1u-Ea1u^2    #variance
repa1u <- cumsum(fSa1u)   #Fonction de répartition
plot(repa1u, xlim = c(0,5000), ylim = c(0,1)) 
sla1u <- function(d) sum(fSa1u*pmax(i-d,0))
Esla1u <- sapply((0:500)/10,sla1u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa1l <- sapply(0:20, dens.conj, a=1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa1l <- ffst(fMa1l)
sum(fSa1l)          #bel et bien une densité
Ea1l <- sum(fSa1l*i)      #espérance                     
E2a1l <- sum(fSa1l*i^2)    #deuxième moment
Vara1l <- E2a1l-Ea1l^2    #variance
repa1l <- cumsum(fSa1l)   #Fonction de répartition
plot(repa1l, xlim = c(0,5000), ylim = c(0,1)) 
sla1l <- function(d) sum(fSa1l*pmax(i-d,0))
Esla1l <- sapply((0:500)/10,sla1l) #Stop-Loss
##moyenne des deux
Ea1 <- mean(c(Ea1u,Ea1l))      #espérance
E2a1 <- mean(c(E2a1u,E2a1l))    #deuxième moment
Vara1 <- mean(c(Vara1u,Vara1l)) #variance
Esla1 <- (Esla1u+Esla1l)/2      #stop-loss
repa1 <- ((repa1u+repa1l)/2)[1:5001]     #repartition
rm(repa1u,repa1l,fSa1u,fSa1l)

##alpha = 10
##upper
densiteB <- densite.upperB
fMa10u <- sapply(0:20, dens.conj, a=10)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa10u <- ffst(fMa10u)
sum(fSa10u)          #bel et bien une densité
Ea10u <- sum(fSa10u*i)      #espérance                     
E2a10u <- sum(fSa10u*i^2)    #deuxième moment
Vara10u <- E2a10u-Ea10u^2    #variance
repa10u <- cumsum(fSa10u)   #Fonction de répartition
plot(repa10u, xlim = c(0,5000), ylim = c(0,1)) 
sla10u <- function(d) sum(fSa10u*pmax(i-d,0))
Esla10u <- sapply((0:500)/10,sla10u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa10l <- sapply(0:20, dens.conj, a=10)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa10l <- ffst(fMa10l)
sum(fSa10l)          #bel et bien une densité
Ea10l <- sum(fSa10l*i)      #espérance                     
E2a10l <- sum(fSa10l*i^2)    #deuxième moment
Vara10l <- E2a10l-Ea10l^2    #variance
repa10l <- cumsum(fSa10l)   #Fonction de répartition
plot(repa10l, xlim = c(0,5000), ylim = c(0,1)) 
sla10l <- function(d) sum(fSa10l*pmax(i-d,0))
Esla10l <- sapply((0:500)/10,sla10l) #Stop-Loss
##moyenne des deux
Ea10 <- mean(c(Ea10u,Ea10l))      #espérance
E2a10 <- mean(c(E2a10u,E2a10l))    #deuxième moment
Vara10 <- mean(c(Vara10u,Vara10l)) #variance
Esla10 <- (Esla10u+Esla10l)/2      #stop-loss
repa10 <- ((repa10u+repa10l)/2)[1:5001]     #repartition
rm(repa10u,repa10l,fSa10u,fSa10l)

##alpha = 30
##upper
densiteB <- densite.upperB
fMa30u <- sapply(0:20, dens.conj, a=30)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa30u <- ffst(fMa30u)
sum(fSa30u)          #bel et bien une densité
Ea30u <- sum(fSa30u*i)      #espérance                     
E2a30u <- sum(fSa30u*i^2)    #deuxième moment
Vara30u <- E2a30u-Ea30u^2    #variance
repa30u <- cumsum(fSa30u)   #Fonction de répartition
plot(repa30u, xlim = c(0,5000), ylim = c(0,1)) 
sla30u <- function(d) sum(fSa30u*pmax(i-d,0))
Esla30u <- sapply((0:500)/10,sla30u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa30l <- sapply(0:20, dens.conj, a=30)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa30l <- ffst(fMa30l)
sum(fSa30l)          #bel et bien une densité
Ea30l <- sum(fSa30l*i)      #espérance                     
E2a30l <- sum(fSa30l*i^2)    #deuxième moment
Vara30l <- E2a30l-Ea30l^2    #variance
repa30l <- cumsum(fSa30l)   #Fonction de répartition
plot(repa30l, xlim = c(0,5000), ylim = c(0,1)) 
sla30l <- function(d) sum(fSa30l*pmax(i-d,0))
Esla30l <- sapply((0:500)/10,sla30l) #Stop-Loss
##moyenne des deux
Ea30 <- mean(c(Ea30u,Ea30l))      #espérance
E2a30 <- mean(c(E2a30u,E2a30l))    #deuxième moment
Vara30 <- mean(c(Vara30u,Vara30l)) #variance
Esla30 <- (Esla30u+Esla30l)/2      #stop-loss
repa30 <- ((repa30u+repa30l)/2)[1:5001]     #repartition
rm(repa30u,repa30l,fSa30u,fSa30l)


##Indépendance
#La fonction de répartition conjointe devient
rep.conj <- function(a, x) prod(pbinom(x, size=1, prob=0.05)) #paramètre a ne change rien, je le mets pour pouvoir réutiliser les fonctions précédentes
##upper
densiteB <- densite.upperB
fMIndu <- sapply(0:20, dens.conj, a=1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSIndu <- ffst(fMIndu)
sum(fSIndu)          #bel et bien une densité
EIndu <- sum(fSIndu*i)      #espérance                     
E2Indu <- sum(fSIndu*i^2)    #deuxième moment
VarIndu <- E2Indu-EIndu^2    #variance
repIndu <- cumsum(fSIndu)   #Fonction de répartition
plot(repIndu, xlim = c(0,5000), ylim = c(0,1)) 
slIndu <- function(d) sum(fSIndu*pmax(i-d,0))
EslIndu <- sapply((0:500)/10,slIndu) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMIndl <- sapply(0:20, dens.conj, a=1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSIndl <- ffst(fMIndl)
sum(fSIndl)          #bel et bien une densité
EIndl <- sum(fSIndl*i)      #espérance                     
E2Indl <- sum(fSIndl*i^2)    #deuxième moment
VarIndl <- E2Indl-EIndl^2    #variance
repIndl <- cumsum(fSIndl)   #Fonction de répartition
plot(repIndl, xlim = c(0,5000), ylim = c(0,1)) 
slIndl <- function(d) sum(fSIndl*pmax(i-d,0))
EslIndl <- sapply((0:500)/10,slIndl) #Stop-Loss
##moyenne des deux
EInd <- mean(c(EIndu,EIndl))      #espérance
E2Ind <- mean(c(E2Indu,E2Indl))    #deuxième moment
VarInd <- mean(c(VarIndu,VarIndl)) #variance
EslInd <- (EslIndu+EslIndl)/2      #stop-loss
repInd <- ((repIndu+repIndl)/2)[1:5001]     #repartition
rm(repIndu,repIndl,fSIndu,fSIndl)


##Tableau des résultats
TableEx3 <- data.frame(Alpha = c('Indépendante', 1,10,30), Espérance = c(EInd,Ea1,Ea10,Ea30), DeuxiemeMoment = c(E2Ind,E2a1,E2a10,E2a30), Variance = c(VarInd,Vara1,Vara10,Vara30))
TableEx3


###Exemple 4

##Cook-Johnson
##La fonction de répartition conjointe est donc
rep.conj <- function(a, x) cook.johnson(a, pbinom(x, size=1, prob=0.05))

##alpha = 13.38
##upper
densiteB <- densite.upperB
fMa13u <- sapply(0:20, dens.conj, a=13.38)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa13u <- ffst(fMa13u)
sum(fSa13u)          #bel et bien une densité
Ea13u <- sum(fSa13u*i)      #espérance                     
E2a13u <- sum(fSa13u*i^2)    #deuxième moment
Vara13u <- E2a13u-Ea13u^2    #variance
repa13u <- cumsum(fSa13u)   #Fonction de répartition
plot(repa13u, xlim = c(0,5000), ylim = c(0,1)) 
sla13u <- function(d) sum(fSa13u*pmax(i-d,0))
Esla13u <- sapply((0:500)/10,sla13u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa13l <- sapply(0:20, dens.conj, a=13.38)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa13l <- ffst(fMa13l)
sum(fSa13l)          #bel et bien une densité
Ea13l <- sum(fSa13l*i)      #espérance                     
E2a13l <- sum(fSa13l*i^2)    #deuxième moment
Vara13l <- E2a13l-Ea13l^2    #variance
repa13l <- cumsum(fSa13l)   #Fonction de répartition
plot(repa13l, xlim = c(0,5000), ylim = c(0,1)) 
sla13l <- function(d) sum(fSa13l*pmax(i-d,0))
Esla13l <- sapply((0:500)/10,sla13l) #Stop-Loss
##moyenne des deux
Ea13 <- mean(c(Ea13u,Ea13l))      #espérance
E2a13 <- mean(c(E2a13u,E2a13l))    #deuxième moment
Vara13 <- mean(c(Vara13u,Vara13l)) #variance
Esla13 <- (Esla13u+Esla13l)/2      #stop-loss
repa13 <- ((repa13u+repa13l)/2)[1:5001]     #repartition
rm(repa13u,repa13l,fSa13u,fSa13l)



##alpha = 2.71
##upper
densiteB <- densite.upperB
fMa2u <- sapply(0:20, dens.conj, a=2.71)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa2u <- ffst(fMa2u)
sum(fSa2u)          #bel et bien une densité
Ea2u <- sum(fSa2u*i)      #espérance                     
E2a2u <- sum(fSa2u*i^2)    #deuxième moment
Vara2u <- E2a2u-Ea2u^2    #variance
repa2u <- cumsum(fSa2u)   #Fonction de répartition
plot(repa2u, xlim = c(0,5000), ylim = c(0,1)) 
sla2u <- function(d) sum(fSa2u*pmax(i-d,0))
Esla2u <- sapply((0:500)/10,sla2u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa2l <- sapply(0:20, dens.conj, a=2.71)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa2l <- ffst(fMa2l)
sum(fSa2l)          #bel et bien une densité
Ea2l <- sum(fSa2l*i)      #espérance                     
E2a2l <- sum(fSa2l*i^2)    #deuxième moment
Vara2l <- E2a2l-Ea2l^2    #variance
repa2l <- cumsum(fSa2l)   #Fonction de répartition
plot(repa2l, xlim = c(0,5000), ylim = c(0,1)) 
sla2l <- function(d) sum(fSa2l*pmax(i-d,0))
Esla2l <- sapply((0:500)/10,sla2l) #Stop-Loss
##moyenne des deux
Ea2 <- mean(c(Ea2u,Ea2l))      #espérance
E2a2 <- mean(c(E2a2u,E2a2l))    #deuxième moment
Vara2 <- mean(c(Vara2u,Vara2l)) #variance
Esla2 <- (Esla2u+Esla2l)/2      #stop-loss
repa2 <- ((repa2u+repa2l)/2)[1:5001]     #repartition
rm(repa2u,repa2l,fSa2u,fSa2l)



##alpha = 37.72
##upper
densiteB <- densite.upperB
fMa37u <- sapply(0:20, dens.conj, a=37.72)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa37u <- ffst(fMa37u)
sum(fSa37u)          #bel et bien une densité
Ea37u <- sum(fSa37u*i)      #espérance                     
E2a37u <- sum(fSa37u*i^2)    #deuxième moment
Vara37u <- E2a37u-Ea37u^2    #variance
repa37u <- cumsum(fSa37u)   #Fonction de répartition
plot(repa37u, xlim = c(0,5000), ylim = c(0,1)) 
sla37u <- function(d) sum(fSa37u*pmax(i-d,0))
Esla37u <- sapply((0:500)/10,sla37u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMa37l <- sapply(0:20, dens.conj, a=37.72)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSa37l <- ffst(fMa37l)
sum(fSa37l)          #bel et bien une densité
Ea37l <- sum(fSa37l*i)      #espérance                     
E2a37l <- sum(fSa37l*i^2)    #deuxième moment
Vara37l <- E2a37l-Ea37l^2    #variance
repa37l <- cumsum(fSa37l)   #Fonction de répartition
plot(repa37l, xlim = c(0,5000), ylim = c(0,1)) 
sla37l <- function(d) sum(fSa37l*pmax(i-d,0))
Esla37l <- sapply((0:500)/10,sla37l) #Stop-Loss
##moyenne des deux
Ea37 <- mean(c(Ea37u,Ea37l))      #espérance
E2a37 <- mean(c(E2a37u,E2a37l))    #deuxième moment
Vara37 <- mean(c(Vara37u,Vara37l)) #variance
Esla37 <- (Esla37u+Esla37l)/2      #stop-loss
repa37 <- ((repa37u+repa37l)/2)[1:5001]     #repartition
rm(repa37u,repa37l,fSa37u,fSa37l)



###Gumbel

##Définissons d'abord la copule de Gumbel (a = param.dép., u = vecteur)
gumbel <- function(a, u) exp(-1*((sum((-1*log(u))^(a)))^(1/a)))
##La fonction de répartition conjointe est donc
rep.conj <- function(a, x) gumbel(a, pbinom(x, size=1, prob=0.05))

##alpha = 1.1
##upper
densiteB <- densite.upperB
fMg11u <- sapply(0:20, dens.conj, a=1.1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg11u <- ffst(fMg11u)
sum(fSg11u)          #bel et bien une densité
Eg11u <- sum(fSg11u*i)      #espérance                     
E2g11u <- sum(fSg11u*i^2)    #deuxième moment
Varg11u <- E2g11u-Eg11u^2    #variance
repg11u <- cumsum(fSg11u)   #Fonction de répartition
plot(repg11u, xlim = c(0,5000), ylim = c(0,1)) 
slg11u <- function(d) sum(fSg11u*pmax(i-d,0))
Eslg11u <- sapply((0:500)/10,slg11u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMg11l <- sapply(0:20, dens.conj, a=1.1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg11l <- ffst(fMg11l)
sum(fSg11l)          #bel et bien une densité
Eg11l <- sum(fSg11l*i)      #espérance                     
E2g11l <- sum(fSg11l*i^2)    #deuxième moment
Varg11l <- E2g11l-Eg11l^2    #variance
repg11l <- cumsum(fSg11l)   #Fonction de répartition
plot(repg11l, xlim = c(0,5000), ylim = c(0,1)) 
slg11l <- function(d) sum(fSg11l*pmax(i-d,0))
Eslg11l <- sapply((0:500)/10,slg11l) #Stop-Loss
##moyenne des deux
Eg11 <- mean(c(Eg11u,Eg11l))      #espérance
E2g11 <- mean(c(E2g11u,E2g11l))    #deuxième moment
Varg11 <- mean(c(Varg11u,Varg11l)) #variance
Eslg11 <- (Eslg11u+Eslg11l)/2      #stop-loss
repg11 <- ((repg11u+repg11l)/2)[1:5001]     #repartition
rm(repg11u,repg11l,fSg11u,fSg11l)



##alpha = 1.5
##upper
densiteB <- densite.upperB
fMg15u <- sapply(0:20, dens.conj, a=1.5)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg15u <- ffst(fMg15u)
sum(fSg15u)          #bel et bien une densité
Eg15u <- sum(fSg15u*i)      #espérance                     
E2g15u <- sum(fSg15u*i^2)    #deuxième moment
Varg15u <- E2g15u-Eg15u^2    #variance
repg15u <- cumsum(fSg15u)   #Fonction de répartition
plot(repg15u, xlim = c(0,5000), ylim = c(0,1)) 
slg15u <- function(d) sum(fSg15u*pmax(i-d,0))
Eslg15u <- sapply((0:500)/10,slg15u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMg15l <- sapply(0:20, dens.conj, a=1.5)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg15l <- ffst(fMg15l)
sum(fSg15l)          #bel et bien une densité
Eg15l <- sum(fSg15l*i)      #espérance                     
E2g15l <- sum(fSg15l*i^2)    #deuxième moment
Varg15l <- E2g15l-Eg15l^2    #variance
repg15l <- cumsum(fSg15l)   #Fonction de répartition
plot(repg15l, xlim = c(0,5000), ylim = c(0,1)) 
slg15l <- function(d) sum(fSg15l*pmax(i-d,0))
Eslg15l <- sapply((0:500)/10,slg15l) #Stop-Loss
##moyenne des deux
Eg15 <- mean(c(Eg15u,Eg15l))      #espérance
E2g15 <- mean(c(E2g15u,E2g15l))    #deuxième moment
Varg15 <- mean(c(Varg15u,Varg15l)) #variance
Eslg15 <- (Eslg15u+Eslg15l)/2      #stop-loss
repg15 <- ((repg15u+repg15l)/2)[1:5001]     #repartition
rm(repg15u,repg15l,fSg15u,fSg15l)


##alpha = 2,5
##upper
densiteB <- densite.upperB
fMg25u <- sapply(0:20, dens.conj, a=2.5)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg25u <- ffst(fMg25u)
sum(fSg25u)          #bel et bien une densité
Eg25u <- sum(fSg25u*i)      #espérance                     
E2g25u <- sum(fSg25u*i^2)    #deuxième moment
Varg25u <- E2g25u-Eg25u^2    #variance
repg25u <- cumsum(fSg25u)   #Fonction de répartition
plot(repg25u, xlim = c(0,5000), ylim = c(0,1)) 
slg25u <- function(d) sum(fSg25u*pmax(i-d,0))
Eslg25u <- sapply((0:500)/10,slg25u) #Stop-Loss
##lower
densiteB <- densite.lowerB
fMg25l <- sapply(0:20, dens.conj, a=2.5)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
fSg25l <- ffst(fMg25l)
sum(fSg25l)          #bel et bien une densité
Eg25l <- sum(fSg25l*i)      #espérance                     
E2g25l <- sum(fSg25l*i^2)    #deuxième moment
Varg25l <- E2g25l-Eg25l^2    #variance
repg25l <- cumsum(fSg25l)   #Fonction de répartition
plot(repg25l, xlim = c(0,5000), ylim = c(0,1)) 
slg25l <- function(d) sum(fSg25l*pmax(i-d,0))
Eslg25l <- sapply((0:500)/10,slg25l) #Stop-Loss
##moyenne des deux
Eg25 <- mean(c(Eg25u,Eg25l))      #espérance
E2g25 <- mean(c(E2g25u,E2g25l))    #deuxième moment
Varg25 <- mean(c(Varg25u,Varg25l)) #variance
Eslg25 <- (Eslg25u+Eslg25l)/2      #stop-loss
repg25 <- ((repg25u+repg25l)/2)[1:5001]     #repartition
rm(repg25u,repg25l,fSg25u,fSg25l)


TableEx4 <- data.frame(Alpha = c("Cook-Jonhson 13.38","Cook-Jonson 2.71", "Cook-Jonhson 37.72", "Gumbel 1.1", "Gumbel 1.5", "Gumbel 2.5"), Espérance = c(Ea13,Ea2,Ea37,Eg11,Eg15,Eg25), DeuxièmeMoment = c(E2a13,E2a2,E2a37,E2g11,E2g15,E2g25), Variance = c(Vara13,Vara2,Vara37,Varg11,Varg15,Varg25))
TableEx4

xsl <- (0:500)/10
xrep <- (0:5000)/100

donnees.stoplosscook <- rbind(data.frame(alpha = 1, stop.loss = Esla1, absc = xsl),data.frame(alpha = 2, stop.loss = Esla2, absc = xsl),data.frame(alpha = 10, stop.loss = Esla10, absc = xsl),data.frame(alpha = 13, stop.loss = Esla13, absc = xsl),data.frame(alpha = 30, stop.loss = Esla30, absc = xsl),data.frame(alpha = 37, stop.loss = Esla37, absc = xsl))

donnees.stoplossgumbel <- rbind(data.frame(alpha = 1.1, stop.loss = Eslg11, absc = xsl),data.frame(alpha = 1.5, stop.loss = Eslg15, absc = xsl),data.frame(alpha = 2.5, stop.loss = Eslg25, absc = xsl))

donnees.repcook <- rbind(data.frame(alpha = 1, repartition = repa1, absc = xrep),data.frame(alpha = 2, repartition = repa2, absc = xrep),data.frame(alpha = 10, repartition = repa10, absc = xrep),data.frame(alpha = 13, repartition = repa13, absc = xrep),data.frame(alpha = 30, repartition = repa30, absc = xrep),data.frame(alpha = 37, repartition = repa37, absc = xrep))

donnees.repgumbel <- rbind(data.frame(alpha = 1.1, repartition = repg11, absc = xrep),data.frame(alpha = 1.5, repartition = repg15, absc = xrep),data.frame(alpha = 2.5, repartition = repg25, absc = xrep))

library(readr)
library(tidyverse)
library(latex2exp)
library(ggplot2)

col <- hcl.colors(5, "SunsetDark")

to_graph <- donnees.repcook
to_graph %>%  ggplot(aes(x = absc)) + 
  geom_line(aes(y = repartition, color = alpha, group = alpha), alpha = 0.6, lwd = 1.25) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),low = col[1], high = tail(col, 1)) + 
  theme_bw() + 
  labs(x = "x", y = TeX("$\\F_S(x)$"),
       title = TeX("Fonction de répartition $F_S(x)$ selon le paramètre de dépendance $\\alpha$ )"),
       subtitle = "Copule de Cook-Johnson") 

to_graph <- donnees.stoplosscook
to_graph %>%  ggplot(aes(x = absc)) + 
  geom_line(aes(y = stop.loss, color = alpha, group = alpha), alpha = 0.6, lwd = 1.25) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),low = col[1], high = tail(col, 1)) + 
  theme_bw() + 
  labs(x = "d", y = TeX("$\\pi_S(d)$"),
       title = TeX("Stop-Loss $\\pi_S(d)$ selon le paramètre de dépendance $\\alpha$ )"),
       subtitle = "Copule de Cook-Johnson") 

to_graph <- donnees.repgumbel
to_graph %>%  ggplot(aes(x = absc)) + 
  geom_line(aes(y = repartition, color = alpha, group = alpha), alpha = 0.6, lwd = 1.25) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),low = col[1], high = tail(col, 1)) + 
  theme_bw() + 
  labs(x = "x", y = TeX("$\\F_S(x)$"),
       title = TeX("Fonction de répartition $F_S(x)$ selon le paramètre de dépendance $\\alpha$ )"),
       subtitle = "Copule de Gumbel") 

to_graph <- donnees.stoplossgumbel
to_graph %>%  ggplot(aes(x = absc)) + 
  geom_line(aes(y = stop.loss, color = alpha, group = alpha), alpha = 0.6, lwd = 1.25) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),low = col[1], high = tail(col, 1)) + 
  theme_bw() + 
  labs(x = "d", y = TeX("$\\pi_S(d)$"),
       title = TeX("Stop-Loss $\\pi_S(d)$ selon le paramètre de dépendance $\\alpha$ )"),
       subtitle = "Copule de Gumbel") 



## Article Oli et Rostan Cossette et al 2018 ----
## Var
VaR <- function(dens, x, kappa){
  repart <-  cumsum(dens)
  min(x[repart >= kappa])
}

## tVaR
TVaR <- function(dens,  x, kappa){
  Varrrr <- VaR(dens, x, kappa)
  x_zero <- x * (x > Varrrr) 
  esp <- sum(x_zero * dens)
  
  x_var_id <- which(x == min(x[x >= Varrrr]))
  Fx <- sum(dens[1:x_var_id])
  
  return((esp + Varrrr * (Fx - kappa))/(1 - kappa))
}

setwd("C:/Users/olico/OneDrive - Université Laval/ULaval/10 - H22/ACT-7102 Théorie du risque avancée/Pres1")

## expectation
esp <- sum(to_graph$f_s * to_graph$x)
var <- sum(to_graph$f_s * (to_graph$x - esp)^2)
VaR_9 <- VaR(dens = to_graph$f_s, to_graph$x, kappa = 0.9)
VaR_999 <- VaR(dens = to_graph$f_s, to_graph$x, kappa = 0.999)
TVaR_9 <- TVaR(dens = to_graph$f_s, x = to_graph$x, kappa = 0.9)
TVaR_999 <- TVaR(dens = to_graph$f_s,
                 x = to_graph$x, kappa = 0.999)

to_show <- data.frame("E_S" = esp,
                      "var_S" = var,
                      "VaR_9_S" = VaR_9,
                      "VaR_999_S" = VaR_999,
                      "TVaR_9_S" = TVaR_9,
                      "TVaR_999_S" = TVaR_999)       
knitr::kable(to_show)
### EXEMPLE 7

i <- 1:4
qi <- 0.1*i

## copule de Franck
franck.copule <- function(alpha, u) {
  ## u est le vecteur (u1, u2, u3, u4)
  (-1 / alpha) * log(1 - prod(sapply(u, function(ui) 1 - exp(-ui * alpha))) / (1 - exp(- alpha))^3)
}

## Fonction de répartition conjointe :
conj.fun <- function(alpha, x) {
  ## x = (x1, x2, x3, x4)
  franck.copule(alpha = alpha, u = pbinom(x, size = 10, prob = qi))
}

## Fonction de densité conjointe
dens.conj <- function(alpha, K, h) {
  ## K = (k1, k2, k3, k4)
  i <- c(0, 1)
  sum(sapply(i, FUN = function(i1) {
    sapply(i, FUN = function(i2) {  
      sapply(i, FUN = function(i3) {
        sapply(i, FUN = function(i4) {
          I <- c(i1, i2, i3, i4)
          x <- sapply(1:4, FUN = function(p) (K[p] - I[p]) * h)
          (- 1)^sum(i1, i2, i3, i4) * conj.fun(alpha = alpha, x = x)
        })
      })
    })  
  }))
  
}

## Fonction de densité de S
dens.S <- function(alpha, k, h) {
  ## k != (k1, k2, k3, k4)
  sum(unlist(sapply(0:k, FUN = function(k1) {
    u1 <- k - k1
    sapply(0:u1, FUN = function(k2) {
      u2 <- k - k1 - k2
      sapply(0:u2, FUN = function(k3) {
        K <- c(k1, k2, k3, k - k1 - k2 - k3)
        dens.conj(alpha = alpha, K = K, h = h)
      })
    })
  })))
}


## YESSSSS CA MARCHE :))))
dens.S(1, 10, 1)
dens.S(3, 10, 1)
dens.S(6, 10, 1)

## Espérance d'ordre r de S : E[S^r] 
esp.S <- function(alpha, r, h = 1) {
  i <- 0:10
  sum(unlist(sapply(i, FUN = function(k1) {
    sapply(i, FUN = function(k2) {  
      sapply(i, FUN = function(k3) {
        sapply(i, FUN = function(k4) {
          K <- c(k1, k2, k3, k4)
          {(sum(K) * h)^r}  * dens.conj(alpha = alpha, K = K, h = h)
        })
      })
    })  
  })))
}

## YESSSSS CA MARCHE :))))
esp.S(1, 1)
esp.S(3)
esp.S(6)

## Variance de S
esp.S(alpha = 1, r = 2) - (esp.S(alpha = 1, r = 1))^2
esp.S(alpha = 3, r = 2) - (esp.S(alpha = 3, r = 1))^2

vect_Sf <- function(alpha, kmax = 50, h = 1) cumsum(sapply(0:kmax, FUN = function(k) dens.S(alpha = alpha, k = k, h = h)))
vect_S <- vect_Sf(alpha = 1)
VaR(vect_S, kappa = 0.9)
VaR(vect_S, kappa = 0.999)

##============================================ TVaR de S
## Espérance modifiée de S   
esp.S.mod <- function(alpha, kappa, h = 1) {
  i <- 0:10
  vrs <- VaR(vect_S, kappa)
  sum(unlist(sapply(i, FUN = function(k1) {
    sapply(i, FUN = function(k2) {  
      sapply(i, FUN = function(k3) {
        sapply(i, FUN = function(k4) {
          K <- c(k1, k2, k3, k4)
          s <- sum(K) * h
          max(s - vrs, 0) * dens.conj(alpha = alpha, K = K, h = h)
        })
      })
    })  
  })))
}

# ## VRAI :)
# TVaR <- function(vect_S, alpha, kappa) {
#   es <- esp.S.mod(alpha, kappa)
#   vr <- VaR(vect_S, kappa = kappa)
#   
#   vr + es / (1 - kappa)
#   #(es + VaR(vect_S, kappa = kappa) * (fvr - kappa)) / (1 - kappa)
# }

TVaR(vect_S, alpha = 1, kappa = 0.9)

## Ex 7 ---------------------------------------------------------------------------------


ff1 <- dnbinom(0:200, 2, 0.8)
ff2 <- dnbinom(0:200, 2, 0.9)
ff3 <- dnbinom(0:200, 3, 0.85)
ff4 <- dnbinom(0:200, 3, 0.95)
matff <- rbind(ff1, ff2, ff3, ff4)
v.n <- c(60, 40, 70, 30) 

fft.nrisks<-function(matff,v.n,m=14)
{
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx)))
  fts<-(ftx)^v.n[1]
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx)))
    fts<-fts*(ftx^v.n[i])
  }
  ffs <- Re(fft(fts, TRUE))/aa
  return(ffs)
}


ffs <- fft.nrisks(matff = matff, v.n = v.n)
ffs[c(71, 81, 91, 101)]

#====== step 1  
## Trouvons theta.star
epsi <- 1e-10

## Fonction de densité de Theta
f.theta <- function(t, alpha) (1 - exp(- alpha))^t / (t * alpha)

## Fonction de répartition de Theta
F.theta <- function(thmax = 200, alpha) {
  cumsum(sapply(1:thmax, FUN = function(t) f.theta(t, alpha)))
}

theta.star <- function(alph) VaR(dens = f.theta(1:1e6, alpha = alph),
                                 x = 1:1e6,
                                 kappa = 1 - epsi)

#====== step 2
theta <- 1
alpha <- 6

#====== step 3
## Fonction de laplace de Theta
Laplace <- function(t, alpha) {
  (-1/alpha) * log(1 - (1 - exp(- alpha)) * exp(- t))
}

## Fonction inverse de Laplace de Theta
Laplace.inv <- function(u, alpha) {
  -log((1 - exp(- alpha * u)) / (1 - exp(- alpha)))
}
## Fonction de repartition de X sachant theta
F_Xi_theta <- function(theta, alpha, q) {
  cum <- cumsum(dbinom(0:10, size = 10, prob = q))
  exp(- theta * Laplace.inv(cum, alpha = alpha))
}

F_Xi_theta(theta, alpha, q = 0.1)

F_X_theta <- function(theta, alpha, qi) {
  lapply(qi, FUN = function(q) F_Xi_theta(theta, alpha, q))
}

q <- as.list(1:4/10)

F_X_theta(theta, alpha, q)

#====== step 4
## Fonction de densité de X sachant theta
f_xi_theta <-  function(i, theta, alpha) {
  Ft <- F_X_theta(theta, alpha, q)[[i]]
  c(Ft[1], diff(Ft))
}

f_x_theta <- function(theta, alpha) {
  lapply(1:4, FUN = function(i) f_xi_theta(i, theta, alpha))
}

## Vérification
f_x_theta(theta, alpha) %>% lapply(sum)

#====== step 5
## Fonction de densité de S sachant Theta

## Algorithme FFT : Somme de n v.a. discrètes indépendantes
fft.nrisks<-function(matff, v.n, m = 14) {
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx)))
  fts<-(ftx)^v.n[1]
  
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx)))
    fts <- fts*(ftx^v.n[i])
  }
  
  ffs <- Re(fft(fts, TRUE))/aa
  return(ffs)
}

matff <- function(theta, alpha) matrix(data = unlist(f_x_theta(theta, alpha)),
                                       nrow = 4,
                                       ncol = 11,
                                       byrow = T)
v.n <- rep(1, 4)

f_s_theta <- function(theta, alpha) fft.nrisks(matff(theta, alpha), v.n)
sum(f_s_theta(theta, alpha))

vect_theta <- 1:theta.star(alpha)

##========== step 7 
vect_f_theta <- sapply(vect_theta, FUN = function(theta) f.theta(theta, alpha))


fs_thet <- sapply(vect_theta, function(t) {f_s_theta(t, alpha)})
fs <- apply(fs_thet, 1, function(f) sum(f * vect_f_theta))
sum(fs)  

vect_s <- seq(fs) - 1
es <- sum(vect_s * fs)
es2 <- sum(vect_s^2 * fs)

alpha <- as.list(1:7)
names(alpha) = sapply(alpha, function(a) paste0("alph_", a))

fs7 <- lapply(alpha, function(a){
tic()
  vect_theta <- 1:theta.star(a)
  vect_f_theta <- sapply(vect_theta, FUN = function(theta) f.theta(theta, a))
  fs_thet <- sapply(vect_theta, function(t) {f_s_theta(t, a)})
  fs <- apply(fs_thet, 1, function(f) sum(f * vect_f_theta))
time <- toc()
return(list("f" = fs,
            "alpha" = a,
            "time" = time$toc - time$tic))
})
fs_7 <- lapply(fs7, function(f) f$f)
times <- sapply(1:length(fs7), function(i) fs7[[i]]$time)
alphas <- sapply(1:length(fs7), function(i) fs7[[i]]$a)

library(scales)
library(latex2exp)

# "#FF3300" red
# "#0099FF" blue
  ggsave("plot/ex7_time.png",
         data.frame(times, alphas) %>%
           ggplot(aes(x = alphas, y = times)) +
           geom_line(lwd = 2, alpha = 0.4, col = "#0099FF") + 
           geom_point(size = 4, col = "#0099FF") +
           scale_y_time(labels = time_format("%M:%S")) + 
           theme_classic() + 
           theme_bw() + 
           labs(x = TeX("$\\alpha$"), y = "Temps de calculs (minutes)",
                title = TeX("Temps requis pour obtenir $f_S(s)$ en fonction de $\\alpha$"),
                subtitle = "La copule Frank est utilisée") + 
           scale_x_continuous(breaks = 1:7)
  )

stats_7 <- lapply(fs7, function(l){
  fs <- l$f
  
  to_keep <- c(1, 3, 6)
  
  if(!(l$alpha %in% to_keep)) return(NULL)
  vect_s <- seq(fs) - 1
  
  esp <- sum(vect_s * fs)
  esp2 <- sum(vect_s^2 * fs)
  var <- esp2 - esp^2
  VaR_09 <- VaR(fs, vect_s, kappa = 0.9)
  VaR_099 <- VaR(fs, vect_s, kappa = 0.999)
  TVaR_09 <- TVaR(fs, vect_s, kappa = 0.9)
  TVaR_099 <- TVaR(fs, vect_s, kappa = 0.999)
  
  to_ret <- c(esp, var, VaR_09, VaR_099, TVaR_09, TVaR_099)
  names(to_ret) <- c("esp", "var", "VaR09", "VaR0999", "TVaR09", "TVaR0999")
  return(to_ret)
})

stats_7

stats_7[sapply(stats_7, is.null)] <- NULL
stats_7

to_output <- as.data.frame.list(stats_7)

library(xtable)

xtable(to_output, digits = 5)

# ### Code Chris ----------------
# kmax <- 2^12
# eps <- 1e-10
# alph <- 0.9
# f_theta <- function(k) (1 - alph) * alph^(k - 1)
# F_theta <- function(k) 1 - alph^k
# laplace_inv <- function(u) log((1 - alph)/u + alph)
# theta_star <- ceiling(log(eps) / log(alph))
# fs <- rep(0, kmax)
# for(theta in 1:theta_star) {
#   Fx_theta <- exp(-theta * laplace_inv(pbinom(0:10, 10, 0.1)))
#   fx_theta <- diff(c(0, Fx_theta))
#   phix <- fft(c(fx_theta, rep(0, kmax - 11)))
#   phis <- phix^100
#   fs_theta <- Re(fft(phis, inverse = TRUE))/kmax
#   fs <- fs + f_theta(theta) * fs_theta
# }
# sum(fs)
# sum((1:kmax - 1) * fs)
# sum((1:kmax - 1)^2 * fs) - (sum((1:kmax - 1) * fs))^2
# (VaR9 <- min(which(cumsum(fs) > 0.9)) - 1)
# (VaR999 <- min(which(cumsum(fs) > 0.999)) - 1)
# round((sum((VaR9 + 1):(kmax - 1) * fs[(VaR9 + 2):kmax]) + VaR9 * (sum(fs[1:(VaR9 + 1)]) - 0.9))/(1 - 0.9), 3)
# round((sum((VaR999 + 1):(kmax - 1) * fs[(VaR999 + 2):kmax]) + VaR999 * (sum(fs[1:(VaR999 + 1)]) - 0.999))/(1 - 0.999), 3)

# > sum((1:kmax - 1) * fs)
# [1] 100
# > sum((1:kmax - 1)^2 * fs) - (sum((1:kmax - 1) * fs))^2
# [1] 2792.839
# >
#   > (VaR9 <- min(which(cumsum(fs) > 0.9)) - 1)
# [1] 172
# > (VaR999 <- min(which(cumsum(fs) > 0.999)) - 1)
# [1] 242
# >
#   > round((sum((VaR9 + 1):(kmax - 1) * fs[(VaR9 + 2):kmax]) + VaR9 * (sum(fs[1:(VaR9 + 1)]) - 0.9))/(1 - 0.9), 3)
# [1] 192.113
# > round((sum((VaR999 + 1):(kmax - 1) * fs[(VaR999 + 2):kmax]) + VaR999 * (sum(fs[1:(VaR999 + 1)]) - 0.999))/(1 - 0.999), 3)
# [1] 250.154


## Ex 8 ---------------------------------------------------------------------------------

# theta is a geom
# X are bin(10, 0.1) 

library(tictoc)
library(tidyverse)

## Laplace de theta pour l'exemple 9 (inutilisé)
psi <- function(t, alpha) (1 - alpha)/(exp(t) - alpha)
## Laplace de theta pour l'exemple 9 (inutilisé)
psi_m1 <- function(u, alpha) log((1 - alpha)/u + alpha)

## param 
ni <- 10
qi <- 0.1


psi_m1(3, 0.5)
psi(psi_m1(3, 0.5), 0.5)

## next fonction de répartition de X
p_X_8 <- function(x) pbinom(x, ni, qi)

## next fonction de répartition de X sachant theta
p_X_theta_8 <- function(x,
                      theta,
                      alpha) {
  
  exp(-theta * psi_m1(p_X_8(x), alpha))
}

## test
p_X_theta_8(4, 3, 0.5)

## fonction de densité de x sachant theta discrètisé
f_xtheta_8 <- function(xmax, h, theta, alpha, method = "lower"){
    return(
      sapply(0:(xmax/h), function(t){
      if(t == 0) return(p_X_theta_8(0, theta, alpha))
      return(p_X_theta_8(t * h, theta, alpha) - p_X_theta_8((t - 1)*h, theta, alpha))
    })
    )
}

alpha <- as.list(seq(0, 0.99, length.out = 1e2))
h <- 1
n_risk <- 100

fs8 <- lapply(alpha, function(a){

tic()
# on cherche theta*
espilon <- 1e-10

theta_star <- qnbinom(1 - espilon, 1, 1 - a) + 1

## theta étoile
values_theta <- 1:theta_star
names(values_theta) <- sapply(values_theta, function(t) paste0("S_theta_",
                                                               ifelse(nchar(t) <=1, "0", ""),
                                                               t))
## On optimise pour avoir xmax
xmax <- 10

## step6
f_s_theta <- sapply(values_theta, function(thet){
  # step4
  f_xi_theta <- f_xtheta_8(xmax, h, thet, a, method = "lower")
  
  # step 5
  aa <- 2^10 # length vector
  nx <- length(f_xi_theta) # length right now
  ftx <- fft(c(f_xi_theta, rep(0, aa - nx))) # fonction caractéristique
  f_x_cum <- Re(fft(ftx^n_risk, TRUE))/aa # On inverse et conserve la partie réelle
  return(f_x_cum)
})

f_s_theta %>% apply(2, sum)
values_x <- seq(from = 0, by = h, length.out = nrow(f_s_theta))

f_theta <- dnbinom(values_theta - 1, 1, 1 - a)
## Étape 7 
f_s <- apply(f_s_theta, 1, function(ft) sum(ft * f_theta))
time <- toc()
return(list("f" = f_s,
            "alpha" = a,
            "time" = time$toc - time$tic))
})


 
times <- sapply(1:length(fs8), function(i) fs8[[i]]$time)
alphas <- sapply(1:length(fs8), function(i) fs8[[i]]$a)

data.frame(times, alphas) %>%
  ggplot(aes(x = alphas, y = times)) +
  geom_line()

library(scales)
library(latex2exp)
library(tidyverse)
ggsave("plot/ex8_time.png",
       data.frame(times, alphas) %>%
         ggplot(aes(x = alphas, y = times)) +
         geom_line(lwd = 2, alpha = 0.4, col = "#0099FF") + 
         geom_point(size = 2, col = "#0099FF") +
         scale_y_time(labels = time_format("%M:%S")) + 
         scale_x_continuous(labels = scales::percent) + 
         theme_classic() + 
         theme_bw() + 
         labs(x = TeX("$\\alpha$"), y = "Temps de calculs (minutes)",
              title = TeX("Temps requis pour obtenir $f_S(s)$ en fonction de $\\alpha$"),
              subtitle = "La copule AMH est utilisée")
)

to_keep <- c(0, 0.5, 0.9)
names(fs8) <- sapply(alpha, function(a) paste0("a_", a))

stats_8 <- lapply(fs8, function(dens){
 f <- dens$f
if(!(dens$alpha %in% to_keep)) return(NULL)
values <- seq(f) - 1
 
 esp <- sum(values * f)
 esp2 <- sum(values^2 * f)
 var <- esp2 - esp^2
 VaR_09 <- VaR(f, values, kappa = 0.9)
 VaR_099 <- VaR(f, values, kappa = 0.999)
 TVaR_09 <- TVaR(f, values, kappa = 0.9)
 TVaR_099 <- TVaR(f, values, kappa = 0.999)
 
 to_ret <- c(esp, var, VaR_09, VaR_099, TVaR_09, TVaR_099)
 names(to_ret) <- c("esp", "var", "VaR09", "VaR0999", "TVaR09", "TVaR0999")
 return(to_ret)
})

stats_8[sapply(stats_8, is.null)] <- NULL
stats_8
str(stats_8)
to_output <- as.data.frame.list(stats_8)
to_output
library(xtable)

xtable(to_output, digits = 3)


## Ex 9 ---------------------------------------------------------------------------------
# theta is a geom
# X are exp (0.1) 

## Laplace de theta pour l'exemple 9 (inutilisé)
psi <- function(t, alpha) (1 - alpha)/(exp(t) - alpha)

## next fonction de répartition de X sachant theta
p_X_theta <- function(x,
                      theta,
                      alpha) {
  
  if(x <= 0) return(0)
  
   ((1 - exp(-0.1 * x))/(1 - alpha * exp(-0.1 * x)))^theta
}

## test
p_X_theta(4, 3, 0.5)

## fonction de densité de x sachant theta discrètisé
f_xtheta <- function(xmax, h, theta, alpha, method = "lower"){
  if(method == "lower"){
    return(sapply(0:(xmax/h), function(t){
        p_X_theta(t * h, theta, alpha) - p_X_theta((t - 1)*h, theta, alpha)
      }))
  }
  if(method == "upper"){
    return(sapply(0:(xmax/h), function(t){
        p_X_theta((t + 1) *h, theta, alpha) - p_X_theta(t * h, theta, alpha)
      }))
  }
}

alpha <- 0.5
# on cherche theta*
espilon <- 1e-8
(theta_star <- qnbinom(1 - espilon, 1, 1 - alpha) + 1)

h <- 0.1

## theta étoile
values_theta <- 1:theta_star
names(values_theta) <- sapply(values_theta, function(t) paste0("S_theta_",
                                                               ifelse(nchar(t) <=1, "0", ""),
                                                               t))

## fonction à optimiser pour obtenir une bonne valeur xmax qui n'est pas trop grosse
to_optim <- function(x){
  min(
    sapply(1:theta_star,
           function(t) abs(ifelse(p_X_theta(x, theta = t, alpha) == 0.999,
                                  0,
                                  log(p_X_theta(x, theta = t, alpha) - 0.999))))
      )
}
## On optimise pour avoir xmax
xmax <- optimize(to_optim, interval = c(20, 5e2))$minimum

## step6
f_s_theta <- sapply(values_theta, function(thet){
  # step4
  f_xi_theta <- f_xtheta(xmax, h, thet, alpha, method = "lower")
  
  # step 5
  aa <- 2^14 # length vector
  nx <- length(f_xi_theta) # length right now
  ftx <- fft(c(f_xi_theta, rep(0, aa - nx))) # fonction caractéristique
  f_x_cum <- Re(fft(ftx^40, TRUE))/aa # On inverse et conserve la partie réelle
  
  return(f_x_cum)
})

library(tidyverse)


f_s_theta %>% apply(2, sum)
values_x <- seq(from = 0, by = h, length.out = nrow(f_s_theta))

f_theta <- dnbinom(values_theta - 1, 1, 1 - alpha)
## Étape 7 
f_s <- apply(f_s_theta, 1, function(a) sum(a * f_theta))


## tests
sum(f_s)
F_s <- cumsum(f_s)
c(length(f_s), length(values_x))



## graphs
to_graph <- data.frame(f_s_theta,
                       "f_s" = f_s,
                       "x" = values_x)



## expectation
(esp <- sum(to_graph$f_s * to_graph$x))
(var <- sum(to_graph$f_s * (to_graph$x - esp)^2))
(VaR_9 <- VaR(dens = to_graph$f_s, to_graph$x, kappa = 0.9))
(VaR_999 <- VaR(dens = to_graph$f_s, to_graph$x, kappa = 0.999))
(TVaR_9 <- TVaR(dens = to_graph$f_s, x = to_graph$x, kappa = 0.9))
(TVaR_999 <- TVaR(dens = to_graph$f_s,
                  x = to_graph$x, kappa = 0.999))


to_show <- data.frame("E_S" = esp,
                      "var_S" = var,
                      "VaR_9_S" = VaR_9,
                      "VaR_999_S" = VaR_999,
                      "TVaR_9_S" = TVaR_9,
                      "TVaR_999_S" = TVaR_999)     

library(xtable)

xtable(to_show)
to_graph %>% ggplot(aes(x = x, y = f_s)) + 
  geom_line()

to_graph %>% ggplot(aes(x = x)) + 
  geom_line(aes(y = S_theta_01), col = "red") + 
  geom_line(aes(y = S_theta_14), col = "blue") + 
  geom_line(aes(y = S_theta_27), col = "green") +
  geom_line(aes(y = f_s), col = "black")

to_graph2 <- to_graph %>% reshape2::melt(id.vars = c("x", "f_s")) %>% 
  mutate(thet = variable %>% substr(9, 10) %>% as.numeric) 
  
#install.packages("latex2exp")
library(latex2exp)

# couleurs 


col <- hcl.colors(length(values_theta), "SunsetDark")

ggsave("plot/ex9_agg.png",
to_graph2 %>%  ggplot(aes(x = x)) + 
  geom_line(aes(y = value, color = thet, group = thet), alpha = 0.6, lwd = 1.25) +
  scale_colour_gradient(name = TeX("Valeur de $\\theta$"),low = col[1], high = tail(col, 1)) + 
  geom_line(aes(y = f_s), col = "black", lwd = 2) +
  theme_bw() + 
  labs(x = "s", y = "Densité",
       title = TeX("fonctions de densité $f_S(s)$ (ou $f_{S|\\Theta = \\theta}(s)$ )")) + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(labels = scales::dollar)
)

write_csv(to_graph2, "graph1.csv")
## pour un graph plus cool 
values_alpha <- seq(0.05, 0.95, length.out = 50)
names(values_alpha) <- sapply(values_alpha, function(n){
  tag <- substr(as.character(n), 3, 4)
  paste0("fs_a_", tag)
})

espilon <- 1e-4
h <- 0.1

f_s_a <- sapply(values_alpha, function(a){
  # on cherche theta*

  (theta_star <- qnbinom(1 - espilon, 1, 1 - a) + 1)
  

  
  values_theta <- 1:theta_star
  
  to_optim <- function(x){
    min(
      sapply(1:theta_star,
             function(t) abs(ifelse(p_X_theta(x, theta = t, a) == 0.99,
                                    0,
                                    log(p_X_theta(x, theta = t, a) - 0.99))))
    )
  }
  xmax <- optimize(to_optim, interval = c(20, 5e2))$minimum
  
## step6
f_s_theta <- sapply(values_theta, function(thet){
    # step4
    f_xi_theta <- f_xtheta(xmax, h, thet, a, method = "lower")
    
    # step 5
    aa <- 2^14 # length vector
    nx <- length(f_xi_theta) # length right now
    ftx <- fft(c(f_xi_theta, rep(0, aa - nx))) # fonction caractéristique
    f_x_cum <- Re(fft(ftx^40, TRUE))/aa # On inverse et conserve la partie réelle
    
    return(f_x_cum)
  })
  
  
  
  f_theta <- dnbinom(values_theta - 1, 1, 1 - a)
  f_s <- apply(f_s_theta, 1, function(fa) sum(fa * f_theta))
  return(f_s)
})

apply(f_s_a, 2, sum)

values_x2 <- seq(from = 0, by = h, length.out = nrow(f_s_a))

dat <- data.frame("x" = values_x2,
                  f_s_a)

dat2 <- dat %>% reshape2::melt(id.vars = "x") %>% 
  mutate(alpha = (variable %>% substr(6, 8) %>% as.numeric)/100) 

write_csv(dat2, "table2.csv")

col_vi <- hcl.colors(5, "Viridis")

ggsave("plot/ex9_dens.png",
dat2 %>%  ggplot(aes(x = x,
                     color = alpha,
                     group = factor(variable))) + 
  geom_line(aes(y = value), lwd = 2, alpha = 0.6) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                        low = col_vi[1],
                        high = tail(col_vi, 1)) + 
  theme_bw() + 
  labs(x = "s",
       y = "Densité",
       title = TeX("$f_S(s)$ selon la valeur $\\alpha$"),
       subtitle = "La copule AMH est utilisée") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(labels = scales::dollar)
)
F_s_a <- apply(f_s_a, 2, cumsum)

dat3 <- data.frame("x" = values_x2,
                   F_s_a) %>% reshape2::melt(id.vars = "x") %>% 
  mutate(alpha = (variable %>% substr(6, 8) %>% as.numeric)/100) 

write_csv(dat3, "table3.csv")

col_pl <- hcl.colors(5, "Plasma")

ggsave("plot/ex9_repart.png",
dat3 %>% 
  ggplot(aes(x = x,
             color = alpha,
             group = factor(variable))) + 
  geom_line(aes(y = value), lwd = 2, alpha = 0.8) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                        low = col_pl[1],
                        high = tail(col_pl, 1)) + 
  theme_bw() + 
  labs(x = "s",
       y = TeX("$F_S(s)$"),
       title = TeX("$F_S(s)$ selon la valeur $\\alpha$"),
       subtitle = "La copule AMH est utilisée") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(labels = scales::dollar)
)

values_kappa <- seq(0.001, 0.999, length.out = 1e3)

VaRs <- apply(f_s_a, 2, function(f){
  sapply(values_kappa, function(k) VaR(f, values_x2, k))
})

TVaRs <- apply(f_s_a, 2, function(f){
  sapply(values_kappa, function(k) TVaR(f, values_x2, k))
})

values_d <- seq(0, 1e3, by = 10)
stop_losses <- apply(f_s_a, 2, function(f){
  sapply(values_d, function(d){
    stop_loss_values <- pmax(values_x2 - d, 0)
    return(sum(f * stop_loss_values))
  })
})


VaRs2 <- VaRs %>% as.data.frame %>% mutate(kap = values_kappa) %>%
  reshape::melt(id.vars = "kap") %>% 
  mutate(VaR = value,
         alpha = (variable %>% substr(6, 8) %>% as.numeric)/100) %>% 
  select(-variable, -value)

TVaRs2 <- TVaRs %>% as.data.frame %>% mutate(kap = values_kappa) %>%
  reshape::melt(id.vars = "kap") %>% 
  mutate(TVaR = value,
         alpha = (variable %>% substr(6, 8) %>% as.numeric)/100)%>% 
  select(-variable, -value)

stop_losses2 <- stop_losses %>% as.data.frame %>% mutate(d = values_d) %>%
  reshape::melt(id.vars = "d") %>% 
  mutate(stoploss = value,
         alpha = (variable %>% substr(6, 8) %>% as.numeric)/100)%>% 
  select(-variable, -value)

tables <- merge(VaRs2, TVaRs2, by = c("kap", "alpha"))

write_csv(tables, "table_var1.csv")

col2 <- hcl.colors(length(values_theta), "Purple-Orange")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark"
ggsave("plot/ex9_VaR.png",
tables %>%  ggplot(aes(x = kap,
                     color = alpha,
                     group = factor(alpha))) + 
  geom_line(aes(y = VaR), lwd = 2, alpha = 0.8) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                        low = col2[1],
                        high = tail(col2, 1)) + 
  theme_bw() + 
  labs(x = TeX("$\\kappa$"),
       y = "Quantile",
       title = TeX("$VaR_{\\kappa}(S)$ la valeur $\\alpha$ dans la copule AMH")) + 
  scale_y_continuous(labels = scales::dollar) + 
  scale_x_continuous(labels = scales::percent)
)


col3 <- hcl.colors(length(values_theta), "Zissou1")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark"
ggsave("plot/ex9_TVaR.png",
tables %>%  ggplot(aes(x = kap,
                       color = alpha,
                       group = factor(alpha))) + 
  geom_line(aes(y = TVaR), lwd = 2, alpha = 0.8) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                        low = col3[1],
                        high = tail(col3, 1)) + 
  theme_bw() + 
  labs(x = TeX("$\\kappa$"),
       y = "TVaR",
       title = TeX("$TVaR_{\\kappa}(S)$ selon la valeur $\\alpha$ dans la copule AMH")) + 
  scale_y_continuous(labels = scales::dollar) + 
  scale_x_continuous(labels = scales::percent)
)

col4 <- hcl.colors(length(values_theta), "Batlow")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark",  "Spectral"
ggsave("plot/ex9_stoploss.png",
stop_losses2 %>% 
  ggplot(aes(x = d,
             color = alpha,
             group = factor(alpha))) + 
  geom_line(aes(y = stoploss), lwd = 2, alpha = 0.8) +
  scale_colour_gradient(name = TeX("Valeur de $\\alpha$"),
                        low = col4[1],
                        high = tail(col4, 1)) + 
  theme_bw() + 
  labs(x = TeX("d"),
       y = TeX("$\\pi_S(d)$"),
       title = TeX("$\\pi_{S}(d)$ selon la valeur $\\alpha$ dans la copule AMH")) + 
  scale_y_continuous(labels = scales::dollar) + 
  scale_x_continuous(labels = scales::dollar)
)
write_csv(stop_losses2, "table_stoploss9.csv")

## EX BEN

###Code pour copier
###EXEMPLE 3 ###


##COOK-JONHSON

##Définissons d'abord la copule de Cook-Johnson (a = param.dép., u = vecteur)

cook.johnson <- function(a, u) (sum(u^(-a))-19)^(-1/a)

##La fonction de répartition conjointe est donc
rep.conj <- function(a, x) cook.johnson(a, pbinom(x, size=1, prob=0.05))

##La fonction de densité conjointe est donc (on emploie équation 16+1 puisque les I sont i.d.) (a = param.dép., x = vecteur de 0 et de 1)
dens.conj <- function(a, x) {
  j <- sum(x)
  
  ToSum <- numeric(j+1)
  for (k in (0:j)){
    ToSum[k+1] <- (-1)^k*factorial(j)/(factorial(k)*factorial(j-k))*rep.conj(a,c(rep(1,j-k),rep(0,20+k-j)))
  }
  sum(ToSum)
}

fM <- sapply(0:20, dens.conj, a=1)*factorial(20)/(factorial(20-0:20)*factorial(0:20))
##La fonction génératrice de M ou M = somme(I) détermine le nombre de fois qu'il faut sommer (puisque le B_i sont iid)

pgf.M <- function(t, dens){
  sum(dens[1:21]*t^(0:20))
}

##Évaluation de la somme composée à l'aide de FFT (on a fgp de M et fdm de B)

#ffst <- function(densM){

#Discrétisation de B
i <- (0:10000)/100
densite.upperB <- sapply(i, function(x) pgamma(x+0.01, shape = 1, rate = 0.5)-
                           pgamma(x, shape = 1, rate = 0.5))

nfft <- 2^20

fB <- c(densite.upperB, rep(0, nfft-length(densite.upperB)))
ftB <- fft(fB)
ftS <- sapply(ftB, pgf.M, dens=fM)
fS <- Re(fft(ftS,inverse = TRUE))/nfft
fS
#}
# 
# fS <- ffst(fM)
sum(fS)

dat <- data.frame(fS, c(i))

library(ggplot2)
ggplot(data, aes(x = i, y = fS)) + geom_line()

