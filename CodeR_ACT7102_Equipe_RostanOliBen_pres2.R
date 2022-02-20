##=========================================================================
## ACT 7102 
##
## Modèles avancés en théorie du risque - session d'Hiver 2022 
##
## Lien de l'article :  https://doi.org/10.1016/j.insmatheco.2011.11.006 
##
## Code R associé aux exemples 12 et 13 de l'article 
##=========================================================================

###=======================================
### Exemple 13
###=======================================

## Trouver l'espérance et la variance des Xi

vect.lambda <- c(0.1, 0.2) 
beta <- 0.1


k1 <- c(0, 0.7, 0.2, 0.1) 
k2 <- c(0, 0.1, 0.4, 0.5) 


phi.Ki <- function(k.i, m = 14) {
  kmax <- 2^m
  ki.long <- c(k.i, rep(0, kmax - length(k.i)))
  fft(ki.long)
}


phi.Mstar <- function(lam, k.i) exp(lam * (phi.Ki(k.i) - 1)) 

f.Mstar <- function(lam, k.i, m = 14) {
  Re(fft(phi.Mstar(lam, k.i), inverse = T))/(2^m)
}

## Espérance de Xi
Esp.X <- function(lam, k.i) {
  fm <- f.Mstar(lam, k.i)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) / beta)
}

## Espérance de Xi^2  
Esp.X2 <- function(lam, k.i) {
  fm <- f.Mstar(lam, k.i)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) * (2 : len) / beta)
}

## Variance de Xi
var.X <- function(lam, k.i) Esp.X2(lam, k.i) - (Esp.X(lam, k.i))^2

## Applications
Esp.X(vect.lambda[1], k1)
Esp.X(vect.lambda[2], k2)

## Je ne comprends pas pourquoi ça ne donne pas 
var.X(vect.lambda[1], k1)
var.X(vect.lambda[2], k2)

dens.Mstar1 <- f.Mstar(vect.lambda[1], k1)
sum(dens.Mstar1)
 


