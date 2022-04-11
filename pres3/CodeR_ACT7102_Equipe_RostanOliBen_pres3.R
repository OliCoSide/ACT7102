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
