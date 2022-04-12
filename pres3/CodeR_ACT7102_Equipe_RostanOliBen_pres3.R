##=========================================================================
## ACT 7102 
##
## Modèles avancés en théorie du risque - session d'Hiver 2022 
##
## Lien de l'article :  https://doi.org/10.1016/j.insmatheco.2007.10.005  
##=========================================================================

###
### Liste des copules choisies :
###

## Copule bivariée de Ali-Mikhail-Haq (AMH)
F.AMH <- function(alpha, u1, u2) {
  ## alpha \in [-1, 1]
  u1 * u2 / (1 - alpha * (1 - u1) * (1 - u2))
}

f.AMH <- function(alpha, u1, u2){
  ## alpha \in [-1, 1]
  a1 <- F.AMH(alpha, u1, u2)
  num <- 1 - alpha + 2 * alpha * a1
  denom <- {1 - alpha * (1 - u1) * (1 - u2) }^2
  
  num / denom
}

## Copule bivariée de Clayton
F.clayton <- function(alpha, u1, u2) {
  ## alpha > 0, (u1, u2) \in [0, 1]^2
  {u1^(-alpha) + u2^(-alpha) - 1}^(-1/alpha)
}

f.clayton <- function(alpha, u1, u2) {
  ## alpha > 0, (u1, u2) \in [0, 1]^2
  a1 <- (1 + alpha)/{(u1 * u2)^(1 + alpha)}
  a2 <- {u1^(-alpha) + u2^(-alpha) - 1}^(-2 -1/alpha)
  
  a1 * a2
}

## Copule bivariée de Frank
F.frank <- function(alpha, u1, u2) {
  ## alpha != 0, (u1, u2) \in [0, 1]^2
  
  f1 <- function(alph = alpha, u) exp(-alph * u) - 1
  
  frac1 <- f1(u1) * f1(u2) / f1(1)
  (- 1 / alpha) * log(1 + frac1)
}

## Utilisation du package copula
install.packages("copula")
library(copula)
library(actuar)


nsim <- 1e5

### copules

table.simul <- function(family, theta, nsim = nsim, seed = 2022) {
  ## family = "clayton", "frank", "amh", "gumbel", "normal"
  set.seed(seed)
  if (family == "normal") {
    cop
  }
}

set.seed(1)  
rCopula(3, normalCopula(0.5, dim = 2)) 

set.seed(1)
cop1 <- archmCopula(family = "clayton", param = 0.5, dim = 2)
rCopula(3, cop1)

set.seed(1) 
rCopula(3, claytonCopula(param = 0.5, dim = 2))


cop.Clayton <- claytonCopula(param = 0.5, dim = 2)
rCopula(nsim, cop.Clayton)


frankCopula(param = 0.5, dim = 2)
gumbelCopula(param = 0.5, dim = 2)
amhCopula(param = 0.5, dim = 2)

##empCopula

