
######## definition des poids sigma, pi, omega et alpha
library(dplyr)
library(actuar)

#### fonction qui ajoute des zeros à un vecteur de longueur totale égale à n
fvecteur <- function(n,p){
  valeur <- vector()
  for (i in 1:length(p)){
    valeur[i]<- p[i]
  }
  for (i in ((length(p)+1):n)){
    valeur[i]<- 0
  }
  return(valeur)
}

## fonction qui calcule le poids sigma_2 pour lindice k
fsigma <- function(k,p1,p2){
  if (k==1){valeur <- 0}
  else {
    #valeur<- sum(sapply(1:k-1, function(i) p1[i]*p2[k-i]))
    valeur <- 0
    for (i in 1:(k-1)) {
      valeur <- valeur+p1[i]*p2[k-i]
    }
  }
  return(valeur)
}

#fsigma(7,c(2,3,4,4,7,8),c(1,7,9,3,4,2))

## fonction qui calcule le poids pi pour lindice k
fpi <- function(j,p,n){
  #(2^(-j+1))*sum(sapply(1:j,
  # function(k)(factorial(j-1)/(factorial(k-1)*factorial(j-k)))*p[k]*sum(sapply(j-k+1:n,function(i) fvecteur(n,p)[i]))))
  valeur <- 0
  for (k in 1:(j)){
    valeur1 <-0
    for (l in (j-k+1):n){
      valeur1 <- valeur1+fvecteur(n,p)[l]
    }
    valeur <- valeur+ (factorial(j-1)/(factorial(k-1)*factorial(j-k)))*fvecteur(n,p)[k]*valeur1
  }
  return((2^(-j+1))*valeur)
}
fpi(5,p1,40)
#fpi(4,p)

## fonction qui calcule le poids omega pour lindice k
fomega <- function(k,p,beta1,beta2,n){
  #sum(sapply(1:k, function(j) p[j]*(factorial(k-1)/(factorial(j-1)*factorial(k-j)))*(beta1/beta2)^j *(1- beta1/beta2)^(k-j)))
  valeur <- 0
  for (j in 1:k){
    valeur <- valeur+fvecteur(n,p)[j]*(factorial(k-1)/(factorial(j-1)*factorial(k-j)))*((beta1/beta2)^j)*((1- beta1/beta2)^(k-j))
  }
  return(valeur)
}

fomega(6,p1,beta1,beta2,10000)

#fomega(5,p,2,3)

## fonction qui calcule le poids alpha pour lindice k
falpha <- function(k,p){
  if (k==1){valeur <- 0}
  else {
    #valeur<- (k-1)*p[k-1]/sum(sapply(1:length(p), function(i) i*p[i]))
    valeur <- 0
    for (i in 1:length(p)){
      valeur <- valeur+(i*p[i])
    }
  }
  return(((k-1)*p[k-1])/valeur)
}

#falpha(6,p)




########### création des vecteurs des poids sigmma, pi, omega et alpha
# vecteur du poids sigma
vsigma <- function(p1,p2){
  valeur <- vector()
  n <- min(length(p1),length(p2))+1
  for (i in 1:n){
    valeur[i] <- fsigma(i,p1,p2)
  }
  return(valeur)
} 

#vsigma(p1,p2)

#vecteur du poids pi
vpi <- function(p,n){
  valeur <- vector()
  for (i in 1:n){
    valeur[i] <- fpi(i,p,n)
  }
  return(valeur)
} 

#vpi(p)

#vecteur du poids omega
vomega <- function(p,beta1,beta2,n){
  valeur <- vector()
  for (i in 1:n){
    valeur[i] <- fomega(i,p,beta1,beta2,n)
  }
  return(valeur)
} 

#vomega(p,1,2)

#vecteur du poids alpha
valpha <- function(p,n){
  valeur <- vector()
  valeur[1]<-0
  for (i in 2:n){
    valeur[i] <- falpha(i,fvecteur(n,p))
  }
  return(valeur)
}

#valpha(p)

### valeur de p2[j] : poids de la somme S2 à l'indice j
p_somme <- function(j,theta,beta1,beta2,p1,p2){
  (1+theta)*fsigma(j,vomega(p1,beta1,2*beta2,n),vomega(p2,beta1,2*beta2,n))
  -theta*fsigma(j,vomega(vpi(p1),2*beta1,2*beta2,n),vomega(p2,beta1,2*beta2,n))
  -theta*fsigma(j,vomega(p1,beta1,2*beta2,n),vpi(p2))
  +theta*fsigma(j,vomega(vpi(p1),2*beta1,2*beta2,n),vpi(p2))
}

vp_somme <- function(theta,beta1,beta2,p1,p2){
  valeur <- vector()
  n <- min(length(p1),length(p2))+1
  for (i in 1:n){
    valeur[i] <- p_somme(i,theta,beta1,beta2,p1,p2)
  }
  return(valeur)
}

################################################## Pour notre example

p1 <- c(0.6,0.4)
p2 <- c(0.3,0.5,0.2)
beta1 <- 0.1
beta2 <- 0.15
theta <- 0.5
n <- 50

### 1er terme du poids de p_2
v1 <- vomega(p1,beta1,2*beta2,n)
v1
v2 <- vomega(p2,beta2,2*beta2,n)
v2
fsigma(1,v1,v2) # pour j=1
fsigma(2,v1,v2) # pour j=2
fsigma(3,v1,v2) # pour j=3
fsigma(4,v1,v2) # pour j=4
fsigma(5,v1,v2) # pour j=5
fsigma(6,v1,v2) # pour j=6

### 2eme terme du poids de p_2
v3 <- vomega(vpi(p1,n),2*beta1,2*beta2,n)
v3
v4 <- vomega(p2,beta2,2*beta2,n)
v4
fsigma(1,v3,v4) # pour j=1
fsigma(2,v3,v4) # pour j=2
fsigma(3,v3,v4) # pour j=3
fsigma(4,v3,v4) # pour j=4
fsigma(5,v3,v4) # pour j=5
fsigma(6,v3,v4) # pour j=6

### 3eme terme du poids de p_2
v5 <- vomega(p1,beta1,2*beta2,n)
v6 <- vpi(p2,n)
fsigma(1,v5,v6) # pour j=1
fsigma(2,v5,v6) # pour j=2
fsigma(3,v5,v6) # pour j=3
fsigma(4,v5,v6) # pour j=4
fsigma(5,v5,v6) # pour j=5
fsigma(6,v5,v6) # pour j=6

### 4eme terme du poids de p_2
v7 <- vomega(vpi(p1,n),2*beta1,2*beta2,n)
v8 <- vpi(p2,n)
fsigma(1,v7,v8) # pour j=1
fsigma(2,v7,v8) # pour j=2
fsigma(3,v7,v8) # pour j=3
fsigma(4,v7,v8) # pour j=4
fsigma(5,v7,v8) # pour j=5
fsigma(6,v7,v8) # pour j=6

### calcul du poids p_2 
p_2_1 <- (1+theta)*fsigma(1,v1,v2) - theta*fsigma(1,v3,v4) - theta*fsigma(1,v5,v6) + theta*fsigma(1,v7,v8)
p_2_1 # pour j=1
p_2_2 <- (1+theta)*fsigma(2,v1,v2) - theta*fsigma(2,v3,v4) - theta*fsigma(2,v5,v6) + theta*fsigma(2,v7,v8)
p_2_2 # pour j=2
p_2_3 <- (1+theta)*fsigma(3,v1,v2) - theta*fsigma(3,v3,v4) - theta*fsigma(3,v5,v6) + theta*fsigma(3,v7,v8)
p_2_3 # pour j=3
p_2_4 <- (1+theta)*fsigma(4,v1,v2) - theta*fsigma(4,v3,v4) - theta*fsigma(4,v5,v6) + theta*fsigma(4,v7,v8)
p_2_4 # pour j=4
p_2_5 <- (1+theta)*fsigma(5,v1,v2) - theta*fsigma(5,v3,v4) - theta*fsigma(5,v5,v6) + theta*fsigma(5,v7,v8)
p_2_5 # pour j=5
p_2_6 <- (1+theta)*fsigma(6,v1,v2) - theta*fsigma(6,v3,v4) - theta*fsigma(6,v5,v6) + theta*fsigma(6,v7,v8)
p_2_6 # pour j=6

# le poids p(2) de S2
p_somme <- function(j,theta,beta1,beta2,p1,p2,n){
  (1+theta)*fsigma(j,v1,v2) - theta*fsigma(j,v3,v4) - theta*fsigma(j,v5,v6) + theta*fsigma(j,v7,v8)
}

p_somme(3,theta,beta1,beta2,p1,p2,n)

vp_somme <- function(theta,beta1,beta2,p1,p2,n){
  valeur <- vector()
  #n <- min(length(p1),length(p2))+1
  for (i in 1:n){
    valeur[i] <- p_somme(i,theta,beta1,beta2,p1,p2,n)
  }
  return(valeur)
}

vp_somme(theta,beta1,beta2,p1,p2,n)

# calcul de la fonction de répartition de S2
FS2 <- function(x,theta,beta1,beta2,p1,p2,n){
  valeur <-0
  for (i in 1:n){
    valeur <- valeur+p_somme(i,theta,beta1,beta2,p1,p2,n)*pgamma(x,i,2*beta2)
  }
  return(valeur)
}

FS2(1:10000,theta,beta1,beta2,p1,p2,n)

# calcul de la VaR de S2


VaR_S2 <- function(kappa,theta,beta1,beta2,p1,p2,n){
    f <- function(x) abs(FS2(x,theta,beta1,beta2,p1,p2,n)-kappa)
    res<-optimize(f, c(0:110))
    VaR<-res$minimum
    return(VaR)
}

p<-c(0.05,0.10,0.50,0.75,0.90,0.95,0.99,0.995,0.999)
##### VaR de S2
VaR_S2(0.05,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.10,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.50,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.75,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.90,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.95,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.99,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.995,theta,beta1,beta2,p1,p2,n) 
VaR_S2(0.999,theta,beta1,beta2,p1,p2,n) 

#### TVaR de S2
TVaR_S2 <- function(kappa,theta,beta1,beta2,p1,p2,n){
  valeur <-0
  for (i in 1:n){
    valeur <- valeur+p_somme(i,theta,beta1,beta2,p1,p2,n)*(i/(2*beta2))*(1-pgamma(VaR_S2(kappa,theta,beta1,beta2,p1,p2,n),i+1,2*beta2))
  }
  return(valeur/(1-kappa))
}

TVaR_S2(0.05,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.10,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.50,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.75,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.90,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.95,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.99,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.995,theta,beta1,beta2,p1,p2,n) 
TVaR_S2(0.999,theta,beta1,beta2,p1,p2,n) 

##################################################### Calcul de la contribution par TVaR

# esperance et Pi de X1 et X2

EX1 <- (p1[1]/beta1)+(p1[2]*2/beta1)
EX1
EX2 <- (p2[1]/beta2)+(p2[2]*2/beta2)+(p2[3]*3/beta2)
EX2

pi <- function(p,beta,n){
  valeur <- 0
  for (i in 1:n){
    valeur <- valeur+i*fpi(i,p,n)
  }
  return(valeur/(2*beta))
}

pi1 <- pi(p1,beta1,n)
pi1
pi2 <- pi(p2,beta2,n)
pi2

### Calcul du poids q2 de la contribution du risque X2 pour le risque total
a1 <- vomega(valpha(p2,n),beta2,2*beta2,n)
a1
a2 <- vomega(p1,beta1,2*beta2,n)
a2
a3 <- valpha(vpi(p2,n),n)
a3
a4 <- vomega(p1,beta1,2*beta2,n)
a4
a5 <- vomega(valpha(p2,n),beta2,2*beta2,n)
a5
a6 <- vomega(vpi(p1,n),2*beta1,2*beta2,n)
a6
a7 <- valpha(vpi(p2,n),n)
a7
a8 <- vomega(vpi(p1,n),2*beta1,2*beta2,n)
a8

q2 <- function(j,theta,beta1,beta2,p1,p2,n){
  (1+theta)*EX2*fsigma(j,a1,a2)-theta*pi2*fsigma(j,a3,a4)-theta*EX2*fsigma(j,a5,a6)+theta*pi2*fsigma(j,a7,a8)
}

vq2 <- function(theta,beta1,beta2,p1,p2,n){
  valeur <- vector()
  for (i in 1:n){
    valeur[i] <- q2(i,theta,beta1,beta2,p1,p2,n)
  }
  return(valeur)
}

vq2(theta,beta1,beta2,p1,p2,n)

### Calcul du poids q1 de la contribution du risque X1 pour le risque total
b1 <- vomega(valpha(p1,n),beta1,2*beta2,n)
b1
b2 <- vomega(p2,beta2,2*beta2,n)
b2
b3 <- vomega(valpha(vpi(p1,n),n),2*beta1,2*beta2,n)
b3
b4 <- vomega(p2,beta2,2*beta2,n)
b4
b5 <- vomega(valpha(p1,n),beta1,2*beta2,n)
b5
b6 <- vpi(p2,n)
b6
b7 <- vomega(valpha(vpi(p1,n),n),2*beta1,2*beta2,n)
b7
b8 <- vpi(p2,n)
b8

q1 <- function(j,theta,beta1,beta2,p1,p2,n){
  (1+theta)*EX1*fsigma(j,b1,b2)-theta*pi1*fsigma(j,b3,b4)-theta*EX1*fsigma(j,b5,b6)+theta*pi1*fsigma(j,b7,b8)
}

vq1 <- function(theta,beta1,beta2,p1,p2,n){
  valeur <- vector()
  for (i in 1:n){
    valeur[i] <- q1(i,theta,beta1,beta2,p1,p2,n)
  }
  return(valeur)
}

vq1(theta,beta1,beta2,p1,p2,n)

### TVaR de X1 et de X2
TVaR_X1 <- function(kappa,theta,beta1,beta2,p1,p2,n){
  valeur <-0
  for (i in 1:n){
    valeur <- valeur+q1(i,theta,beta1,beta2,p1,p2,n)*(1-pgamma(VaR_S2(kappa,theta,beta1,beta2,p1,p2,n),i,2*beta2))
  }
  return(valeur/(1-kappa))
}

TVaR_X2 <- function(kappa,theta,beta1,beta2,p1,p2,n){
  valeur <-0
  for (i in 1:n){
    valeur <- valeur+q2(i,theta,beta1,beta2,p1,p2,n)*(1-pgamma(VaR_S2(kappa,theta,beta1,beta2,p1,p2,n),i,2*beta2))
  }
  return(valeur/(1-kappa))
}

TVaR_X1(0.95,-1,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,-0.8,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,-0.6,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,-0.4,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,-0.2,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,0,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,0.2,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,0.4,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,0.6,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,0.8,beta1,beta2,p1,p2,n)
TVaR_X1(0.95,1,beta1,beta2,p1,p2,n)

##
TVaR_X2(0.95,-1,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,-0.8,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,-0.6,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,-0.4,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,-0.2,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,0,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,0.2,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,0.4,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,0.6,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,0.8,beta1,beta2,p1,p2,n)
TVaR_X2(0.95,1,beta1,beta2,p1,p2,n)

################################ Calcul de la contribution par covariance
  
#####
# rhoooo <- function(pi,betai,theta,beta1,beta2,p1,p2,n){
#   aa <-0
#   bb <-0
#   for (l in 1:n) {
#     aa <- aa+l*fomega(l,pi,betai,beta2,n)
#     bb <- bb+l*p_somme(l,theta,beta1,beta2,p1,p2,n)
#   }
#   cc <-0
#   for (l in 1:n){
#     for (m in 1:n){
#       cc <- cc+l*m*(fpi(l,vomega(p1,beta1,beta2,n),n)/2-fomega(l,p1,beta1,beta2,n))*(fpi(m,p2,n)/2 -fvecteur(n,p2)[m])
#     }
#   }
#   return((aa-aa^2+theta*cc)/(bb-bb^2))
# }
# 
# rhoooo(p1,beta1,-1,beta1,beta2,p1,p2,n)
#####
terme_a <- function(pi,betai,theta,beta1,beta2,p1,p2,n){
  aa <-0
  aaa <-0
  for (l in 1:n) {
    aa <- aa+l*(l+1)*fomega(l,pi,betai,beta2,n)
    aaa <-aaa+l*fomega(l,pi,betai,beta2,n)
  }
  return(aa -aaa^2)
}

terme1_a <- terme_a(p1,beta1,theta,beta1,beta2,p1,p2,n)
terme1_a
terme2_a <- terme_a(p2,beta2,theta,beta1,beta2,p1,p2,n)
terme2_a 

terme_b <- function(theta,beta1,beta2,p1,p2,n){
  bb <-0
  bbb<-0
  for (l in 1:n) {
    bb <- bb+ l*(l+1)*p_somme(l,theta,beta1,beta2,p1,p2,n)
    bbb <- bbb+ l*p_somme(l,theta,beta1,beta2,p1,p2,n)
  }
  return(bb -bbb^2)
}

termebb <-terme_b(theta,beta1,beta2,p1,p2,n)
termebb 

terme_c <- function(theta,beta1,beta2,p1,p2,n){
  cc <- 0
  for (l in 1:n){
    for (m in 1:n){
      cc <- cc+l*m*(fpi(l,vomega(p1,beta1,beta2,n),n)/2 -fomega(l,p1,beta1,beta2,n))*(fpi(m,p2,n)/2 -fvecteur(n,p2)[m])
    }
  }
  return(cc)
}

termec <- terme_c(theta,beta1,beta2,p1,p2,n)
termec #On a trouvé que terme_c = 0.809325

rhoo1 <- function(theta,beta1,beta2,p1,p2,n){
  a <- terme_a(p1,beta1,theta,beta1,beta2,p1,p2,n)
  b <- terme_b(theta,beta1,beta2,p1,p2,n)
  return(a/b + (theta*termec)/b)
}

rhoo1(-1,beta1,beta2,p1,p2,n)
rhoo1(-0.8,beta1,beta2,p1,p2,n)
rhoo1(-0.6,beta1,beta2,p1,p2,n)
rhoo1(-0.4,beta1,beta2,p1,p2,n)

rhoo2 <-   function(theta,beta1,beta2,p1,p2,n){
  a <- terme_a(p2,beta2,theta,beta1,beta2,p1,p2,n)
  b <- terme_b(theta,beta1,beta2,p1,p2,n)
  return(a/b + (theta*termec)/b)
}

rhoo2(-1,beta1,beta2,p1,p2,n)
rhoo2(-0.8,beta1,beta2,p1,p2,n)
rhoo2(-0.6,beta1,beta2,p1,p2,n)
rhoo2(-0.4,beta1,beta2,p1,p2,n)

c1k <- function(k,theta,beta1,beta2,p1,p2,n,kappa){
  fomega(k,p1,beta1,beta2,n)+2*rhoo1(theta,beta1,beta2,p1,p2,n)*p_somme(k,theta,beta1,beta2,p1,p2,n)*(((1-pgamma(VaR_S2(kappa,theta,beta1,beta2,p1,p2,n),k+1,2*beta2))/(1-kappa))-1)
}

c2k <- function(k,theta,beta1,beta2,p1,p2,n,kappa){
  fomega(k,p2,beta2,beta2,n)+2*rhoo2(theta,beta1,beta2,p1,p2,n)*p_somme(k,theta,beta1,beta2,p1,p2,n)*(((1-pgamma(VaR_S2(kappa,theta,beta1,beta2,p1,p2,n),k+1,2*beta2))/(1-kappa))-1)
}

C1k <- function(theta,beta1,beta2,p1,p2,n,kappa){
  valeur <-0
  for (k in 1:n){
    valeur <- valeur+c1k(k,theta,beta1,beta2,p1,p2,n,kappa)*(k/beta2)
  }
  return(valeur)
}

C2k <- function(theta,beta1,beta2,p1,p2,n,kappa){
  valeur <-0
  for (k in 1:n){
    valeur <- valeur+c2k(k,theta,beta1,beta2,p1,p2,n,kappa)*(k/beta2)
  }
  return(valeur)
}

#### pour le risque X1 à kappa=0.95
C1k(-1,beta1,beta2,p1,p2,n,0.95)
C1k(-0.8,beta1,beta2,p1,p2,n,0.95)
C1k(-0.6,beta1,beta2,p1,p2,n,0.95)
C1k(-0.4,beta1,beta2,p1,p2,n,0.95)
C1k(-0.2,beta1,beta2,p1,p2,n,0.95)
C1k(0,beta1,beta2,p1,p2,n,0.95)
C1k(0.2,beta1,beta2,p1,p2,n,0.95)
C1k(0.4,beta1,beta2,p1,p2,n,0.95)
C1k(0.6,beta1,beta2,p1,p2,n,0.95)
C1k(0.8,beta1,beta2,p1,p2,n,0.95)
C1k(1,beta1,beta2,p1,p2,n,0.95)

#### pour le risque X2 à kappa=0.95
C2k(-1,beta1,beta2,p1,p2,n,0.95)
C2k(-0.8,beta1,beta2,p1,p2,n,0.95)
C2k(-0.6,beta1,beta2,p1,p2,n,0.95)
C2k(-0.4,beta1,beta2,p1,p2,n,0.95)
C2k(-0.2,beta1,beta2,p1,p2,n,0.95)
C2k(0,beta1,beta2,p1,p2,n,0.95)
C2k(0.2,beta1,beta2,p1,p2,n,0.95)
C2k(0.4,beta1,beta2,p1,p2,n,0.95)
C2k(0.6,beta1,beta2,p1,p2,n,0.95)
C2k(0.8,beta1,beta2,p1,p2,n,0.95)
C2k(1,beta1,beta2,p1,p2,n,0.95)


##################################################### illustration numérique avec graphes

data1 <- data.frame(Theta=double(),TVaRS=double(),TVaR_X1=double(),C_X1=double(),TVaR_X2=double(),C_X2=double())


theta_valeurs <- c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
ii <-1
for(t in theta_valeurs ){
  data1[ii,1]<- t
  data1[ii,2]<-TVaR_S2(0.95,t,beta1,beta2,p1,p2,n)
  data1[ii,3]<-TVaR_X1(0.95,t,beta1,beta2,p1,p2,n)
  data1[ii,4]<-C1k(t,beta1,beta2,p1,p2,n,0.95)
  data1[ii,5]<-TVaR_X2(0.95,t,beta1,beta2,p1,p2,n)
  data1[ii,6]<-C2k(t,beta1,beta2,p1,p2,n,0.95)
  ii <- ii+1
}
data1

data2 <- data.frame(Kappa=double(),TVaRS=double(),TVaR_X1=double(),C_X1=double(),TVaR_X2=double(),C_X2=double())
###########
kappa_valeurs <- c(0.05,0.10,0.50,0.75,0.90,0.95,0.99,0.995,0.999)
jj <-1
for(k in kappa_valeurs ){
  data2[jj,1]<- k
  data2[jj,2]<-TVaR_S2(k,0.5,beta1,beta2,p1,p2,n)
  data2[jj,3]<-TVaR_X1(k,0.5,beta1,beta2,p1,p2,n)
  data2[jj,4]<-C1k(0.5,beta1,beta2,p1,p2,n,k)
  data2[jj,5]<-TVaR_X2(k,0.5,beta1,beta2,p1,p2,n)
  data2[jj,6]<-C2k(0.5,beta1,beta2,p1,p2,n,k)
  jj <- jj+1
}

data2

library(tidyr)
data1_long <- gather(data1, Allocations, Valeurs, TVaRS:C_X2, factor_key=TRUE)
data1_long

data2_long <- gather(data2, Allocations, Valeurs, TVaRS:C_X2, factor_key=TRUE)
data2_long

data1_long$Theta <- as.factor(data1_long$Theta)
data2_long$Kappa <- as.factor(data2_long$Kappa)
library(ggplot2)
########### graphe barplot de toutes les allocations lorsque theta varie et kappa=0.95
ggplot(data=data1_long, aes(x=Theta, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)+ 
  scale_fill_manual("Allocations", values = c("TVaRS" = "red", "TVaR_X1" = "orange", "C_X1" = "blue", "TVaR_X2" = "yellow", "C_X2" = "purple"))

########### graphe barplot de toutes les allocations lorsque kappa varie et theta=0.5
ggplot(data=data2_long, aes(x=Kappa, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)+ 
  scale_fill_manual("Allocations", values = c("TVaRS" = "red", "TVaR_X1" = "orange", "C_X1" = "blue", "TVaR_X2" = "yellow", "C_X2" = "purple"))

########### graphe barplot de l'allocation par TVaR lorsque theta varie et kappa=0.95
data3 <- data1[,c(1,2,3,5)]
data3_long <- gather(data3, Allocations, Valeurs, TVaRS:TVaR_X2, factor_key=TRUE)
data3_long

data3_long$Theta <- as.factor(data3_long$Theta)
ggplot(data=data3_long, aes(x=Theta, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)

########### graphe barplot de l'allocation par TVaR lorsque kappa varie et theta=0.5
data4 <- data2[,c(1,2,3,5)]
data4_long <- gather(data4, Allocations, Valeurs, TVaRS:TVaR_X2, factor_key=TRUE)
data4_long

data4_long$Kappa <- as.factor(data4_long$Kappa)
ggplot(data=data4_long, aes(x=Kappa, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)

########### graphe barplot de l'allocation par Covariance lorsque theta varie et kappa=0.95
data5 <- data1[,c(1,2,4,6)]
data5_long <- gather(data5, Allocations, Valeurs, TVaRS:C_X2, factor_key=TRUE)
data5_long

data5_long$Theta <- as.factor(data5_long$Theta)
ggplot(data=data5_long, aes(x=Theta, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)

########### graphe barplot de l'allocation par TVaR lorsque kappa varie et theta=0.5
data6 <- data2[,c(1,2,4,6)]
data6_long <- gather(data6, Allocations, Valeurs, TVaRS:C_X2, factor_key=TRUE)
data6_long

data6_long$Kappa <- as.factor(data6_long$Kappa)
ggplot(data=data6_long, aes(x=Kappa, y=Valeurs, fill=Allocations)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)



############################ cette partie est un supplément ###############################
###### calcul des tvar x1 x2 sans tenir compte de la dependance pour comparer avec les allocations quand il y a dependance

# calcul de la fonction de répartition de X1 et X2
FX1 <- function(x,beta1,p1){
  valeur <-0
  for (i in 1:length(p1)){
    valeur <- valeur+p1[i]*pgamma(x,i,beta1)
  }
  return(valeur)
}

FX1(1:10000,beta1,p1)

#####
FX2 <- function(x,beta2,p2){
  valeur <-0
  for (i in 1:length(p2)){
    valeur <- valeur+p2[i]*pgamma(x,i,beta2)
  }
  return(valeur)
}

FX2(1:10000,beta2,p2)

# calcul de la VaR de X1 et X2

VaRX1 <- function(kappa,beta1,p1){
  f <- function(x) abs(FX1(x,beta1,p1)-kappa)
  res<-optimize(f, c(0:110))
  VaR<-res$minimum
  return(VaR)
}

TVaRX1 <- function(kappa,beta1,p1){
  valeur <-0
  for (i in 1:length(p1)){
    valeur <- valeur+p1[i]*(i/beta1)*(1-pgamma(VaRX1(kappa,beta1,p1),i+1,beta1))
  }
  return(valeur/(1-kappa))
}

#####
VaRX2 <- function(kappa,beta2,p2){
  f <- function(x) abs(FX2(x,beta2,p2)-kappa)
  res<-optimize(f, c(0:110))
  VaR<-res$minimum
  return(VaR)
}

TVaRX2 <- function(kappa,beta2,p2){
  valeur <-0
  for (i in 1:length(p2)){
    valeur <- valeur+p2[i]*(i/beta2)*(1-pgamma(VaRX2(kappa,beta2,p2),i+1,beta2))
  }
  return(valeur/(1-kappa))
}


data7 <- data.frame(Kappa=double(),TVaRX1=double(),TVaRX2=double())

kappa_valeurs <- c(0.05,0.10,0.50,0.75,0.90,0.95,0.99,0.995,0.999)
ll <-1
for(k in kappa_valeurs ){
  data7[ll,1]<- k
  data7[ll,2]<-TVaRX1(k,beta1,p1)
  data7[ll,3]<-TVaRX2(k,beta2,p2)
  ll <- ll+1
}

data7

data7_long <- gather(data7, TVaR, Valeurs, TVaRX1:TVaRX2, factor_key=TRUE)
data7_long

################## graphe de TVar(X1) et TVar(X2)
data7_long$Kappa <- as.factor(data7_long$Kappa)
ggplot(data=data7_long, aes(x=Kappa, y=Valeurs, fill=TVaR)) +
  geom_bar(stat="identity", position=position_dodge(),width = 0.5)
