

###################################### Exemple 1 #############################################


###################### tracage de Fs

Fv = function(x) {
  return(pgamma(x,20,0.5) )
}

discr.B <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [2] = pgamma(1,1,0.5)
  for (i in (2:k)) 
  {
    f_g[i+1] = pgamma(i,1,0.5) - pgamma(i-1,1,0.5)
  }
  return(f_g)
}

panjer.binom<-function(nn,qq,ff,smax) {
  # Algorithme de Panjer
  # Cas Binomiale
  # Loi discrete pour B 
  aa<- -qq/(1-qq)
  bb<- -(nn+1)*aa 
  ll<-length(ff) 
  ffs<-(1-qq+qq*ff[1])^nn
  ff<-c(ff,rep(0,smax-ll+1))
  for (i in 1 :smax)
  {
    j<-i+1
    ffs<-c(ffs,(1/(1-aa*ff[1]))*sum(ff[2 :j]*ffs[i :1]*(bb*(1 :i)/i+aa))) }
  return(ffs)
}

Fu = function(x,q) {
  return(cumsum(panjer.binom(20,(1-0.95/(1-q)),discr.B(50),50))[x+1] )
}

Fs = function(x,q) {
  return( (1-q)*Fu(x,q) + q*Fv(x) )
}

x=0:50

plot(x,Fs(x,0),type = "l", lty = 3, ylim = 0:1, col= "red", xlab = "x", ylab = "Fs(x)")
lines(c(0,0,0),c(0,0.4,0.79))
lines(c(0,0),c(0.8,0.86),type = "l", lty = 5)
lines(x,Fs(x,0.015),type = "l", col= "blue")
lines(x,Fs(x,0.025),type = "l", lty = 2 , col= "green")
lines(x,Fs(x,0.045),type = "l", lty = 5, col= "black" )
legend(30, 0.5, legend = c("q0=0", "q0=0.015", "q0=0.025", "q0=0.045"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))


############################ calcul de pi(d)
o=1:50
FFs1 = Fs(o,0)
FFs2 = Fs(o,0.015)
FFs3 = Fs(o,0.025)
FFs4 = Fs(o,0.045)

p1=rep(0,50)
p2=rep(0,50)
p3=rep(0,50)
p4=rep(0,50)

pi1 = function() {
for (d in (1:50)) {
  s=0
  for (i in (d:50)) {
     s= s+ (1-FFs1[i])
 }
   p1[d]=s
}
  p1[1]=2
  return(p1)
}

pi2 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-FFs2[i])
    }
    p2[d]=s
  }
  p2[1]=2
  return(p2)
}
pi3 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-FFs3[i])
    }
    p3[d]=s
  }
  p3[1]=2
  return(p3)
}
pi4 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-FFs4[i])
    }
    p4[d]=s
  }
  p4[1]=2
  return(p4)
}

x=0:49

plot(x,pi1(),type = "l", lty = 3, ylim = c(0,2), col= "red", xlab = "u", ylab = "pi(u)")
lines(x,pi2(),type = "l", col= "blue")
lines(x,pi3(),type = "l", lty = 2, col= "green" )
lines(x,pi4(),type = "l", lty = 5, col= "black" )
legend(30, 1.5, legend = c("q0=0", "q0=0.015", "q0=0.025", "q0=0.045"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))

########################### calcul de la Tvar 

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FFs1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FFs1[i+1] - FFs1[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FFs2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FFs2[i+1] - FFs2[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FFs3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FFs3[i+1] - FFs3[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FFs4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FFs4[i+1] - FFs4[i]
  }
  return(f_g)
}
ind = seq(0,49, by=1)

VaR_k = function(f,k) 
{
  i= min(which(cumsum(f)>=0.95))
  return(i-1)  
}

TVaR <- function(f,  x, kappa){
  Varrrr <- VaR_k(f, kappa)
  x_zero <- x * (x > Varrrr) 
  esp <- sum(x_zero * f)
  x_var_id <- which(x == min(x[x >= Varrrr]))
  Fx <- sum(f[1:x_var_id])
  return((esp + Varrrr * (Fx - kappa))/(1 - kappa))
}

k=seq(0,1, by=0.1)

k = c(0,0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.80,0.85,0.9,0.95,0.995,1)
T1= rep(0,14)
T2= rep(0,14)
T3= rep(0,14)
T4= rep(0,14)

j=1
for (i in (k)) {
  T1[j] = TVaR(f.S1(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T2[j] = TVaR(f.S2(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T3[j] = TVaR(f.S3(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T4[j] = TVaR(f.S4(49),ind,i)
  j = j+1
}

plot(k,T1,type = "l", lty = 3, ylim = c(0,70), col= "red", xlab = "k", ylab = "TVaR_k(S)")
lines(k,T2,type = "l", lty = 1, col= "blue")
lines(k,T3,type = "l", lty = 2, col= "green" )
lines(k,T4,type = "l", lty = 5 , col= "black")
legend(0.2, 50, legend = c("q0=0", "q0=0.015", "q0=0.025", "q0=0.045"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))


####################################### fin exemple 1 #########################################

####################################### Exemple 2 #############################################

function(q0,qj,qi) {
  a=1-qi
  b=1-q0
  c=1-qj
  return(1-a/(b*c))
}


###################### tracage de Fs


###### q0=0


Fc1 = function(qj) {
  F1= rep(0,50)
  for (i in (1:50))
    F1[i] = qj*pgamma(i-1,5,0.5)  +(1-qj)*cumsum(panjer.binom(5,q_tildeij(0,0.006,0.02),discr.B(50),50))[i]
  return(F1)
}

Fc2 = function(qj) {
  F2= rep(0,50)
  for (i in (1:50))
    F2[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0,0.006,0.035),discr.B(50),50))[i]
  return(F2)
}

Fc3 = function(qj) {
  F3= rep(0,50)
  for (i in (1:50))
    F3[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0,0.006,0.05),discr.B(50),50))[i]
  return(F3)
}

Fc4 = function(qj) {
  F4= rep(0,50)
  for (i in (1:50))
    F4[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0,0.006,0.065),discr.B(50),50))[i]
  return(F4)
}

FF1 = Fc1(0.06)
FF2 = Fc2(0.06)
FF3 = Fc3(0.06)
FF4 = Fc4(0.06)

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF1[i+1] - FF1[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF4[i+1] - FF4[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF3[i+1] - FF3[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF2[i+1] - FF2[i]
  }
  return(f_g)
}

ff1 = f.S1(49)
ff2 = f.S2(49)
ff3 = f.S3(49)
ff4 = f.S4(49)

matf = rbind(ff1,ff2,ff3,ff4)

fft.nrisks<-function(matff,v.n,m=14) {
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx))) 
  fts<-(ftx)^v.n
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx))) 
    fts<-fts*(ftx^v.n)
  }
  ffs <- Re(fft(fts, TRUE))/aa 
  return(ffs)
}

Fuu = cumsum(fft.nrisks(matf,v.n = 1 ,m=14))

Fu = function(x) {
  return(Fuu[x+1] )
}


Fs = function(x,q0) {
  return( (1-q0)*Fu(x) + q0*Fv(x) )
}
Fs1 = Fs(x,0)


###### q0=0.005


Fc1 = function(qj) {
  F1= rep(0,50)
  for (i in (1:50))
    F1[i] = qj*pgamma(i-1,5,0.5)  +(1-qj)*cumsum(panjer.binom(5,q_tildeij(0.005,0.006,0.02),discr.B(50),50))[i]
  return(F1)
}

Fc2 = function(qj) {
  F2= rep(0,50)
  for (i in (1:50))
    F2[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.005,0.006,0.035),discr.B(50),50))[i]
  return(F2)
}

Fc3 = function(qj) {
  F3= rep(0,50)
  for (i in (1:50))
    F3[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.005,0.006,0.05),discr.B(50),50))[i]
  return(F3)
}

Fc4 = function(qj) {
  F4= rep(0,50)
  for (i in (1:50))
    F4[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.005,0.006,0.065),discr.B(50),50))[i]
  return(F4)
}

FF1 = Fc1(0.06)
FF2 = Fc2(0.06)
FF3 = Fc3(0.06)
FF4 = Fc4(0.06)

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF1[i+1] - FF1[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF4[i+1] - FF4[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF3[i+1] - FF3[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF2[i+1] - FF2[i]
  }
  return(f_g)
}

ff1 = f.S1(49)
ff2 = f.S2(49)
ff3 = f.S3(49)
ff4 = f.S4(49)

matf = rbind(ff1,ff2,ff3,ff4)

fft.nrisks<-function(matff,v.n,m=14) {
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx))) 
  fts<-(ftx)^v.n
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx))) 
    fts<-fts*(ftx^v.n)
  }
  ffs <- Re(fft(fts, TRUE))/aa 
  return(ffs)
}

Fuu = cumsum(fft.nrisks(matf,v.n = 1 ,m=14))

Fu = function(x) {
  return(Fuu[x+1] )
}


Fs = function(x,q0) {
  return( (1-q0)*Fu(x) + q0*Fv(x) )
}
Fs2 = Fs(x,0.005)


#### q0 = 0.01

Fc1 = function(qj) {
  F1= rep(0,50)
  for (i in (1:50))
    F1[i] = qj*pgamma(i-1,5,0.5)  +(1-qj)*cumsum(panjer.binom(5,q_tildeij(0.01,0.006,0.02),discr.B(50),50))[i]
  return(F1)
}

Fc2 = function(qj) {
  F2= rep(0,50)
  for (i in (1:50))
    F2[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.01,0.006,0.035),discr.B(50),50))[i]
  return(F2)
}

Fc3 = function(qj) {
  F3= rep(0,50)
  for (i in (1:50))
    F3[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.01,0.006,0.05),discr.B(50),50))[i]
  return(F3)
}

Fc4 = function(qj) {
  F4= rep(0,50)
  for (i in (1:50))
    F4[i] = qj*pgamma(i-1,5,0.5)  + (1-qj)*cumsum(panjer.binom(5,q_tildeij(0.01,0.006,0.065),discr.B(50),50))[i]
  return(F4)
}

FF1 = Fc1(0.06)
FF2 = Fc2(0.06)
FF3 = Fc3(0.06)
FF4 = Fc4(0.06)

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF1[i+1] - FF1[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF4[i+1] - FF4[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF3[i+1] - FF3[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF2[i+1] - FF2[i]
  }
  return(f_g)
}

ff1 = f.S1(49)
ff2 = f.S2(49)
ff3 = f.S3(49)
ff4 = f.S4(49)

matf = rbind(ff1,ff2,ff3,ff4)

fft.nrisks<-function(matff,v.n,m=14) {
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx))) 
  fts<-(ftx)^v.n
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx))) 
    fts<-fts*(ftx^v.n)
  }
  ffs <- Re(fft(fts, TRUE))/aa 
  return(ffs)
}

Fuu = cumsum(fft.nrisks(matf,v.n = 1 ,m=14))

Fu = function(x) {
  return(Fuu[x+1] )
}


Fs = function(x,q0) {
  return( (1-q0)*Fu(x) + q0*Fv(x) )
}
Fs3 = Fs(x,0.01)

##### independance

Fc1 = function() {
  F1= rep(0,50)
  for (i in (1:50))
    F1[i] = cumsum(panjer.binom(5,q_tildeij(0,0,0.02),discr.B(50),50))[i]
  return(F1)
}

Fc2 = function(qj) {
  F2= rep(0,50)
  for (i in (1:50))
    F2[i] = cumsum(panjer.binom(5,q_tildeij(0,0,0.035),discr.B(50),50))[i]
  return(F2)
}

Fc3 = function(qj) {
  F3= rep(0,50)
  for (i in (1:50))
    F3[i] = cumsum(panjer.binom(5,q_tildeij(0,0,0.05),discr.B(50),50))[i]
  return(F3)
}

Fc4 = function(qj) {
  F4= rep(0,50)
  for (i in (1:50))
    F4[i] = cumsum(panjer.binom(5,q_tildeij(0,0,0.065),discr.B(50),50))[i]
  return(F4)
}

FF1 = Fc1(0)
FF2 = Fc2(0)
FF3 = Fc3(0)
FF4 = Fc4(0)

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF1[i+1] - FF1[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF4[i+1] - FF4[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF3[i+1] - FF3[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = FF2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = FF2[i+1] - FF2[i]
  }
  return(f_g)
}

ff1 = f.S1(49)
ff2 = f.S2(49)
ff3 = f.S3(49)
ff4 = f.S4(49)

matf = rbind(ff1,ff2,ff3,ff4)

fft.nrisks<-function(matff,v.n,m=14) {
  aa <- 2^m
  nbrisks<-dim(matff)[1]
  fx<-matff[1,]
  nx <- length(fx)
  ftx <- fft(c(fx, rep(0, aa - nx))) 
  fts<-(ftx)^v.n
  for (i in 2:nbrisks)
  {
    fx<-matff[i,]
    nx <- length(fx)
    ftx <- fft(c(fx, rep(0, aa - nx))) 
    fts<-fts*(ftx^v.n)
  }
  ffs <- Re(fft(fts, TRUE))/aa 
  return(ffs)
}

Fuu = cumsum(fft.nrisks(matf,v.n = 1 ,m=14))

Fu = function(x) {
  return(Fuu[x+1] )
}


Fs = function(x,q0) {
  return( (1-q0)*Fu(x) + q0*Fv(x) )
}
Fs4 = Fs(x,0)



##############

x=0:49
plot(x,Fs1,type = "l", lty = 3, ylim = 0:1, col= "red", xlab = "x", ylab = "Fs(x)")
lines(c(0,0),c(0,0.38),type = "l", lty = 3, col= "red" )
lines(c(0,0),c(0.38,0.4),type = "l", lty = 1, , col= "blue")
lines(c(0,0),c(0.4,0.44),type = "l", lty = 2, , col= "green")
lines(x,Fs2,type = "l", col= "blue")
lines(x,Fs3,type = "l", lty = 2 , col= "green")
lines(x,Fs4,type = "l", lty = 5, col= "black" )
legend(30, 0.5, legend = c("q0=0", "q0=0.005", "q0=0.01", "independance"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))

##################### stop loss 

p1=rep(0,50)
p2=rep(0,50)
p3=rep(0,50)
p4=rep(0,50)

pi1 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-Fs1[i])
    }
    p1[d]=s
  }
  return(p1)
}

pi2 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-Fs2[i])
    }
    p2[d]=s
  }
  return(p2)
}
pi3 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-Fs3[i])
    }
    p3[d]=s
  }
  return(p3)
}
pi4 = function() {
  for (d in (1:50)) {
    s=0
    for (i in (d:50)) {
      s= s+ (1-Fs4[i])
    }
    p4[d]=s
  }
  return(p4)
}

x=0:49

plot(x,pi1(),type = "l", lty = 3, ylim = c(0,1.7), col= "red", xlab = "u", ylab = "pi(u)")
lines(x,pi2(),type = "l", col= "blue")
lines(x,pi3(),type = "l", lty = 2, col= "green" )
lines(x,pi4(),type = "l", lty = 5, col= "black" )
legend(30, 1.5, legend = c("q0=0", "q0=0.005", "q0=0.01", "independance"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))

##### calcul de TVaR

f.S1 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = Fs1[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = Fs1[i+1] - Fs1[i]
  }
  return(f_g)
}

f.S2 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = Fs2[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = Fs2[i+1] - Fs2[i]
  }
  return(f_g)
}

f.S3 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = Fs3[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = Fs3[i+1] - Fs3[i]
  }
  return(f_g)
}

f.S4 <- function(k) 
{
  f_g = rep(0,k+1)
  f_g [1] = Fs4[1]
  for (i in (1:k)) 
  {
    f_g[i+1] = Fs4[i+1] - Fs4[i]
  }
  return(f_g)
}

ind = seq(0,49, by=1)

VaR_k = function(f,k) 
{
  i= min(which(cumsum(f)>=0.95))
  return(i-1)  
}

TVaR <- function(f,  x, kappa){
  Varrrr <- VaR_k(f, kappa)
  x_zero <- x * (x > Varrrr) 
  esp <- sum(x_zero * f)
  x_var_id <- which(x == min(x[x >= Varrrr]))
  Fx <- sum(f[1:x_var_id])
  return((esp + Varrrr * (Fx - kappa))/(1 - kappa))
}

k=seq(0,1, by=0.1)

k = c(0,0.1, 0.2, 0.3,0.4,0.5,0.6,0.7,0.80,0.85,0.9,0.95,0.995,1)
T1= rep(0,14)
T2= rep(0,14)
T3= rep(0,14)
T4= rep(0,14)

j=1
for (i in (k)) {
  T1[j] = TVaR(f.S1(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T2[j] = TVaR(f.S2(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T3[j] = TVaR(f.S3(49),ind,i)
  j = j+1
}

j=1
for (i in (k)) {
  T4[j] = TVaR(f.S4(49),ind,i)
  j = j+1
}

plot(k,T1,type = "l", lty = 3, ylim = c(0,70), col= "red", xlab = "k", ylab = "TVaR_k(S)")
lines(k,T2,type = "l", lty = 1, col= "blue")
lines(k,T3,type = "l", lty = 2, col= "green" )
lines(k,T4,type = "l", lty = 5 , col= "black")
legend(0.2, 50, legend = c("q0=0", "q0=0.005", "q0=0.01", "independance"),col = c("red","blue","green","black")  , lty = c(3,1,2,5))
