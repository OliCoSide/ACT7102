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
cop1 <- normalCopula(param = 0.7, dim = 2)
u_simulated <- rCopula(100, cop1) 

## répartition de mix erl
p_dist_oli <- function(x){
  sum(dnbinom(0:50, size = 4, mu = 5) *
                                       pgamma(x, 0:50 + 1, 0.1))
}

VaR_MixErl <- function(kap){optimize(function(x){log(abs(kap - p_dist_oli(x)))},
                       interval = c(0, 350))$minimum}

## sim x1
sim_x1 <- sapply(u_simulated[, 1], function(u) VaR_MixErl(u))

## diagnostic rapide
mean(sim_x1)
diff_u1 <- u_simulated[, 1] - sapply(sim_x1, function(x) p_dist_oli(x))
hist(diff_u1)
summary(diff_u1)
plot(u_simulated[, 1] ~ sapply(sim_x1, function(x) p_dist_oli(x)))


## sim x2
library(tidyverse)
library(actuar)
sim_x2 <- qgenpareto(u_simulated[, 2], shape1 = 3,
           shape2 = 5,
           scale = 50) 

sims  <- data.frame(u_simulated, x1 = sim_x1, x2 = sim_x2)
names(sims)[1:2] <- c("u1", "u2")
names(sims)
sims %>% ggplot(aes(x = u1, y = u2)) + geom_point()
sims %>% ggplot(aes(x = x1, y = x2)) + geom_point()

sims$u1_hat <- pobs(sims$x1)
sims$u2_hat <- pobs(sims$x2)

library(latex2exp)
g <- rbind(sims %>% 
  select(u1, u2) %>% 
    mutate(hat = "Simulés"),
  sims %>%
    select(u1_hat, u2_hat) %>%
    rename("u1" = "u1_hat",
           "u2" = "u2_hat") %>% 
    mutate(hat = "Empiriques")
) %>% 
  mutate(paired = rep(1:nrow(sims), 2)) %>% 
  ggplot(aes(x = u1, y = u2, color = hat)) + 
  geom_point(alpha = 0.4, size = 1.5) + 
  theme_bw()+
  geom_line(aes(group = paired), color = "black", alpha = 0.15, size = 0.5) + 
  scale_color_brewer(palette = "Dark2") +
  labs(title = TeX("Relation entre les u simulés et les quantiles empiriques $\\hat{u}$"),
       x = TeX("$u_1$"),
       y = TeX("$u_2$")) + 
  theme(legend.title = element_blank())

# ggsave("fig/comp_u1u2.png",
#        g, 
#        width = 15,
#        height = 11,
#        units = "cm")



mat1 <- matrix(c(sims$u1_hat, sims$u2_hat),
               byrow = FALSE, ncol = 2)

## Inverting Kendall's tau
fit.alph1 <- fitCopula(normalCopula(), mat1, method="mpl")
fit.alph1
confint(fit.alph1)

fit.alph2 <- fitCopula(frankCopula(), mat1, method="mpl")
fit.alph2
confint(fit.alph2)

simuls_c_thetan1<-rCopula(1e3, copula = frankCopula(coef(fit.alph2), dim = 2))
plot(simuls_c_thetan1)

empcop1 <- function(u) C.n(X= mat1, u)
empcop1_est <- function(u) C.n(X= simuls_c_thetan1, u)
dist_cop1 <- function(u) sqrt(nrow(sims)) *(empcop1(u) - empcop1_est(u))

(S_n1 <- sum(sapply(1:nrow(sims), function(i) dist_cop1(mat1[i, ])^2)))
(T_n1 <- max(sapply(1:nrow(sims), function(i) abs(dist_cop1(mat1[i, ])))))
mat1[86,]


Sn <- function(sims, copula_to_test = NULL, N = 10 * nrow(sims)){
  simuls2 <- pobs(sims)
  
  empcop1 <- function(u) C.n(X= simuls2, u) # copule empirique
  fit.alph99 <- fitCopula(copula_to_test, simuls2, method="mpl") # meilleur param
  
  simuls_c_thetan1<-rCopula(N, copula = copula_to_test(coef(fit.alph99), dim = 2)) # simuls selon meilleur param
  empcop1_est <- function(u) C.n(X= simuls_c_thetan1, u) # approx meilleur copule frank
  
  S_star_nk <- sum(sapply(1:nrow(simuls2), function(i){ # dist
    (sqrt(nrow(sims)) *(empcop1(simuls2[i, ]) - empcop1_est(simuls2[i, ])))^2 
  }))
  
  return(S_star_nk)
}

N <- 1e4
Sn_star <- sapply(seq(N), function(i){
  ## n simuls
  simuls<-rCopula(1e2, copula = frankCopula(coef(fit.alph2), dim = 2))
  Sn(simuls, copula_to_test = frankCopula())
})

(approx_p_value <- mean(Sn_star > Sn1))


### TEST AVEC K 

best_param <- coef(fit.alph2)
simuls_c_thetan1<-rCopula(1e3, copula = frankCopula(best_param, dim = 2))
empcop1 <- function(u) C.n(X= mat1, u)
empcop1_est <- function(u) C.n(X= simuls_c_thetan1, u)
V_i <- function(i) mean(sapply(1:nrow(simuls_c_thetan1), function(j) prod(simuls_c_thetan1[j,] < simuls_c_thetan1[i,])))
dist_cop2 <- function(u) sqrt(nrow(sims)) *(empcop1(u) - empcop1_est(u))

(S_nK <- sum(sapply(1:nrow(sims), function(i) dist_cop2(simuls_c_thetan1[i, ])^2)))
(T_nK <- max(sapply(1:nrow(sims), function(i) abs(dist_cop2(simuls_c_thetan1[i, ])))))


## visualisation : idées
## https://r-graph-gallery.com/3d-surface-plot.html
## 


sapply(seq(0.01, 0.99, length.out = 99), function(u1){
  sapply(seq(0.01, 0.99, length.out = 99), function(u2){
    empcop1(c(u1, u2))
  })
})

simuls_to_graph <- rCopula(1e3, normalCopula(-0.9))
cop_to_graph <- function(u) C.n(X= simuls_to_graph, u)
to_show_in_3d <-sapply(seq(0.01, 0.99, length.out = 99), function(u1){
  sapply(seq(0.01, 0.99, length.out = 99), function(u2){
    dCopula(c(u1, u2), normalCopula(-0.9))
  })
})
plot_ly(z = to_show_in_3d, type = "surface")
