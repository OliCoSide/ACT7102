##=========================================================================
## ACT 7102 
##
## Modèles avancés en théorie du risque - session d'Hiver 2022 
##
## Lien de l'article :  https://doi.org/10.1016/j.insmatheco.2007.10.005  
##=========================================================================
 
## Utilisation du package copula
# install.packages("copula")
library(copula)
library(tidyverse)
library(actuar)
library(latex2exp)
library(xtable)

## Simulation de U à partir d'une copule

u_simulated <- function(param, real_copula = NULL, nsim = 100) {
  # real_copula = frankCopula, normalCopula, gumbelCopula, amhCopula, claytonCopula 
  rCopula(nsim, copula = real_copula(param = param, dim = 2))
}
#u_simulated(0.5, real_copula = claytonCopula)

## Marginale F1 : répartition de MixErl
p_dist_oli <- function(x){
  poids <- dnbinom(0:50, size = 4, mu = 5)
  sum(poids * pgamma(x, 0:50 + 1, 0.1))
}

VaR_MixErl <- function(kap){
  optimize(function(x) {
    log(abs(kap - p_dist_oli(x)))
            }, interval = c(0, 350))$minimum
}

## Observations X1
sim_x1 <- function(param, real_copula) {
  # real_copula = frankCopula, normalCopula, gumbelCopula, amhCopula, claytonCopula
  u_simul <- u_simulated(param, real_copula = real_copula)
  sapply(u_simul[, 1], function(u) VaR_MixErl(u))
}
#sim_x1(0.5, claytonCopula)  
 
## Marginale F2 : Pareto généralisée
## Simulation de X2
sim_x2 <- function(param, real_copula) {
  u_simul <- u_simulated(param, real_copula = real_copula)
  qgenpareto(u_simul[, 2], shape1 = 3, shape2 = 5, scale = 50)
}
#sim_x2(0.5, claytonCopula)   
   
sims  <- function(param, real_copula) {
  
  ## Simulation de U à partir d'une copule
  u_simul <- u_simulated(param, real_copula = real_copula)
  
  ## Observations X1 et X2 à partir des marginales F1 et F2
  simul_x1 <- sapply(u_simul[, 1], function(u) VaR_MixErl(u))
  simul_x2 <- qgenpareto(u_simul[, 2], shape1 = 3, shape2 = 5, scale = 50)
  
  ## pseudo-observations U1 et U2 (chapeau)
  u1_hat <- pobs(simul_x1)
  u2_hat <- pobs(simul_x2)
 
  to_data <- data.frame(u_simul, x1 = simul_x1, x2 = simul_x2, 
                        u1_hat = u1_hat, u2_hat = u2_hat)
  names(to_data)[1:2] <- c("u1", "u2")
  to_data
}
#names(sims(0.5, claytonCopula))

#sims(0.5, claytonCopula) %>% ggplot(aes(x = u1, y = u2)) + geom_point()
#sims(0.5, normalCopula) %>% ggplot(aes(x = x1, y = x2)) + geom_point()
 

g <- function(param, real_copula) {
  
  to_sims <- sims(param, real_copula)
  
  g1 <- rbind(to_sims %>% 
               select(u1, u2) %>% 
               mutate(hat = "Simulés"),
             to_sims %>%
               select(u1_hat, u2_hat) %>%
               rename("u1" = "u1_hat",
                      "u2" = "u2_hat") %>% 
               mutate(hat = "Empiriques")
  ) %>% 
    mutate(paired = rep(1:nrow(to_sims), 2)) %>% 
    ggplot(aes(x = u1, y = u2, color = hat)) + 
    geom_point(alpha = 0.4, size = 1.5) + 
    theme_bw()+
    geom_line(aes(group = paired), color = "black", alpha = 0.15, size = 0.5) + 
    scale_color_brewer(palette = "Dark2") +
    labs(title = TeX("Relation entre les u simulés et les quantiles empiriques $\\hat{u}$"),
         #subtitle = paste0('Copule ', real_copula()),
         x = TeX("$u_1$"),
         y = TeX("$u_2$")) + 
    theme(legend.title = element_blank())
  
  return(g1)
  
}

## family = "clayton", "frank", "amh", "gumbel", "normal"
g(param = 0.5, amhCopula) 
g(param = 0.5, frankCopula)

fit.alph <- function(param, real_copula, copula_to_test) {
  ## real_copula = frankCopula, normalCopula, gumbelCopula, amhCopula, claytonCopula
  ## copula_to_test = frankCopula(), normalCopula(), gumbelCopula(), amhCopula(), claytonCopula()
  
  to_sims <- sims(param = param, real_copula = real_copula)
  to_mat1 <- matrix(c(to_sims$u1_hat, to_sims$u2_hat),
                    byrow = FALSE, ncol = 2)
 
  fitCopula(copula_to_test, to_mat1, method = "mpl")
}

(fit.alph1 <- fit.alph(0.5, claytonCopula, claytonCopula()))
confint(fit.alph1, level = 0.98)

simuls_c_thetan1 <- function(param, real_copula, copula_to_test = NULL, m_sim = 1e3) {
  ## real_copula, copula_to_test = frankCopula, normalCopula, gumbelCopula, amhCopula, claytonCopula 
  
  fit.alph1 <- fit.alph(param, real_copula, copula_to_test())
  coef <- coef(fit.alph1)
  
  u_simulated(param = coef, real_copula = real_copula, nsim = m_sim)
  rCopula(m_sim, copula = copula_to_test(coef, dim = 2))
}
  
#rCopula(1e3, copula = frankCopula(coef(fit.alph2), dim = 2))
#simuls_c_thetan1(0.5, claytonCopula, normalCopula)
 

Sn_Tn_K <- function(param, real_copula, copula_to_test = NULL) {
  # real_copula, copula_to_test = frankCopula, normalCopula, gumbelCopula, amhCopula, claytonCopula
  to_sims1 <- sims(param = param, real_copula = real_copula)
  to_sims <- matrix(c(to_sims1$u1_hat, to_sims1$u2_hat),
                    byrow = FALSE, ncol = 2)
  
  simuls2 <- pobs(to_sims)
  
  M <- 10 * nrow(to_sims)
  
  empcop1 <- function(u) C.n(X = simuls2, u) # copule empirique
  fit.alph99 <- fitCopula(copula_to_test(dim = 2), simuls2, method = "mpl") # meilleur param

  simuls_c_thetan <- rCopula(M, copula = copula_to_test(coef(fit.alph99), dim = 2)) # simuls selon meilleur param
  empcop1_est <- function(u) C.n(X = simuls_c_thetan, u) # approx meilleur copule  

  ##
  ## TESTS BASÉS SUR LA COPULE EMPIRIQUE
  ##
  
  dist_cop1 <- function(u) sqrt(nrow(to_sims)) *(empcop1(u) - empcop1_est(u))
  
  S_star_nk <- sum(sapply(1:nrow(simuls2), function(i) { # dist
    dist_cop1(simuls2[i, ])^2
  }))

  T_star_nk <- max(sapply(1:nrow(simuls2), function(i) { # dist
    abs(dist_cop1(simuls2[i, ]))
  }))
 
  ##
  ## TESTS BASÉS SUR LA TRANSFORMATION DE KENDALL
  ##
  
  V_i <- function(i) {
    mean(sapply(1:nrow(simuls_c_thetan), function(j) {
      prod(simuls_c_thetan[j,] <= simuls_c_thetan[i,])
    }))
  }
  
  to_V <- sapply(1:M, function(i) V_i(i))
  
  K_theta_n <- rank(to_V)/M
  K_n <- rank(to_V[1:nrow(simuls2)])/nrow(simuls2)
  
  dist_K <- sqrt(nrow(simuls2) / M) * (K_n - K_theta_n)
  
  Sn_K <- sum(dist_K^2)
  Tn_K <- max(abs(dist_K))
    
  return(list(S_star_nk = S_star_nk,
              T_star_nk = T_star_nk,
              Sn_K = Sn_K,
              Tn_K = Tn_K)) 
}

Sn_Tn_K(0.5, normalCopula, frankCopula)

## Ensuite on répète l'expérience N fois 
Sn_Tn_K_star <- function(param, real_copula, copula_to_test, N = 150) {
  
  sapply(seq(N), function(i) {
    Sn_Tn_K(param, real_copula, copula_to_test)
  })
  
}
#cl <- Sn_Tn_K_star(0.5, frankCopula, amhCopula)  

p_values <- function(param, real_copula, copula_to_test, N = 150) {
  values <- Sn_Tn_K_star(param, real_copula, copula_to_test, N = N)
  
  p_value_Sn <- mean(unlist(values[1, ]) > values[1, 1])
  p_value_Tn <- mean(unlist(values[2, ]) > values[2, 1])
  
  p_value_Sn_K <- mean(unlist(values[3, ]) > values[3, 1])
  p_value_Tn_K <- mean(unlist(values[4, ]) > values[4, 1])
  
  return(list(p_Sn = p_value_Sn, 
              p_Tn = p_value_Tn,
              p_Sn_K = p_value_Sn_K, 
              p_Tn_K = p_value_Tn_K))
}

## copule sous H0 : frankCopula
p11 <- p_values(0.2, frankCopula, frankCopula)
p12 <- p_values(0.2, frankCopula, gumbelCopula)
p13 <- p_values(0.2, frankCopula, amhCopula)
p14 <- p_values(0.2, frankCopula, normalCopula)
# p15 <- p_values(0.5, frankCopula, claytonCopula)

## copule sous H0 : normalCopula
p21 <- p_values(0.2, normalCopula, gumbelCopula)
p22 <- p_values(0.2, normalCopula, frankCopula)
#p23 <- p_values(0.5, normalCopula, amhCopula)
p24 <- p_values(0.2, normalCopula, normalCopula)
#p25 <- p_values(0.5, gumbelCopula, claytonCopula)
 

data1 <- data.frame(Copule_H0 = c('Frank', ' ', ' ', ' '),
                    Copule_test = c('Frank', 'Gumbel', 'AMH', 'Normale'),
                    Sn = c(p11$p_Sn, p12$p_Sn, p13$p_Sn, p14$p_Sn),
                    Tn = c(p11$p_Tn, p12$p_Tn, p13$p_Tn, p14$p_Tn),
                    Sn_K = c(p11$p_Sn_K, p12$p_Sn_K, p13$p_Sn_K, p14$p_Sn_K),
                    Tn_K = c(p11$p_Tn_K, p12$p_Tn_K, p13$p_Tn_K, p14$p_Tn_K)
                    )

data2 <- data.frame(Copule_H0 = c('Normal', ' ', ' '),
                    Copule_test = c('Gumbel', 'Frank', 'Normale'),
                    Sn = c(p21$p_Sn, p22$p_Sn, p24$p_Sn),
                    Tn = c(p21$p_Tn, p22$p_Tn, p24$p_Tn),
                    Sn_K = c(p21$p_Sn_K, p22$p_Sn_K, p24$p_Sn_K),
                    Tn_K = c(p21$p_Tn_K, p22$p_Tn_K, p24$p_Tn_K)
                    )

xtable(data1)
xtable(data2)

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
