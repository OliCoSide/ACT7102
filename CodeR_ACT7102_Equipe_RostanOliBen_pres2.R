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

phi.Mstar <- function(lam, k.i, m = 14) exp(lam * (phi.Ki(k.i, m) - 1)) 

f.Mstar <- function(lam, k.i, m = 14) {
  Re(fft(phi.Mstar(lam, k.i, m), inverse = T))/(2^m)
}

## Espérance de Xi
Esp.X <- function(lam, k.i, analytic = FALSE) {
  if(analytic){ # version analytique
    return(
      sum(lam * # Espérance de M
            (k.i * (seq(k.i) - 1)/beta) # espérance de B (pondération des gammas)
      )
    )
  }
  fm <- f.Mstar(lam, k.i)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) / beta)
}

## Espérance de Xi^2  
Esp.X2 <- function(lam, k.i) {
  fm <- f.Mstar(lam, k.i, m = 6)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) * (2 : len) / beta^2)
}

## Variance de Xi
var.X <- function(lam, k.i, analytic = FALSE){
  if(analytic){
    return(
      sum(lam * k.i * (seq(k.i) - 1) * (seq(k.i)) / beta^2)
    )
  }
  Esp.X2(lam, k.i) - (Esp.X(lam, k.i))^2
}

## Applications ---- [OC] analytical form 
Esp.X(vect.lambda[1], k1)
Esp.X(vect.lambda[1], k1, analytic = TRUE)
Esp.X(vect.lambda[2], k2)
Esp.X(vect.lambda[2], k2, analytic = TRUE)

## Je ne comprends pas pourquoi ça ne donne pas ---- [OC] fixed + analytical form 
var.X(vect.lambda[1], k1)
var.X(vect.lambda[1], k1, analytic = TRUE)
var.X(vect.lambda[2], k2)
var.X(vect.lambda[2], k2, analytic = TRUE)

## Moment d'ordre k de Mstar
Moments <- function(lam, k.i, ord) {
  fm <- f.Mstar(lam, k.i, m = 14)
  len <- length(fm)
  
  sum(fm * (0 : (len - 1))^ord)
}

## Variance de Xi : ça marche :)
var.Xi <- function(lam, k.i) { 
  mum <- Moments(lam, k.i, ord = 1) +
    Moments(lam, k.i, ord = 2) -
    (Moments(lam, k.i, ord = 1))^2
  mum / beta^2
}

var.Xi(vect.lambda[1], k1)
var.Xi(vect.lambda[2], k2)
 
dens.Mstar1 <- f.Mstar(vect.lambda[1], k1)
sum(dens.Mstar1 * (0 : (length(dens.Mstar1) - 1)))/beta

sum(dens.Mstar1)
 
## Moment d'ordre 'ord' de Xi
Esp.M <- function(lam, k.i, ord) {
  fm <- f.Mstar(lam, k.i, m = 6)
  
  k_table <- sapply(1:ord, function(i) {
    ((1 + i):(length(fm) - 1 + i)) - 1
  })
  
  k_factor <- unlist(apply(k_table, 1, prod))
  
  sum(fm[-1] * k_factor) / beta^ord
}

## Fonction de répartition de Xi
F.Xi <- function(lam, k.i, x, m = 6) {
  fm <- f.Mstar(lam, k.i, m = 6)
  len <- length(fm)
  
  fm[1] + sum(fm[-1] * pgamma(x, 1 : (len - 1), beta)) 
}

## VaR de Xi
VaR.Xi <- function(lam, k.i, kappa) {
  fm <- f.Mstar(lam, k.i, m = 14)
  if (fm[1] > kappa) return(0)
  optimise(function(x) abs(F.Xi(lam, k.i, x) - kappa), c(0, 200))$minimum
}

VaR.Xi(vect.lambda[1], k1, 0.995)
VaR.Xi(vect.lambda[2], k2, 0.995) # ça ne marche pas

F.Xi(vect.lambda[2], k2, 59.9921) # Résultat de l'article (devrait être 0.995)
F.Xi(vect.lambda[2], k2, VaR.Xi(vect.lambda[2], k2, 0.995)) # Notre résultat 

## TVaR de Xi
TVaR.Xi <- function(lam, k.i, kappa) {
  fm <- f.Mstar(lam, k.i, m = 14)
  vv <- VaR.Xi(lam, k.i, kappa)
  
  pgamma.bar <- function(k) 1 - pgamma(vv, k + 1, beta)
  
  k.sum <- sapply(1:(length(fm) - 1), function(k) {
    fm[k + 1] * k * pgamma.bar(k) / beta
  })
  
  sum(k.sum) / (1 - kappa)
}

TVaR.Xi(vect.lambda[1], k1, 0.995)
TVaR.Xi(vect.lambda[2], k2, 0.995)

## Vérifier que Var(S) = 1050 + 21 440 * alpha_0
sum(sapply(1:10, function(i) {
  sapply(1:10, function(j) (14 * (i <= 5) + 24 * (i > 5))*(14 * (j <= 5) + 24 * (j > 5)))
}))

2 * sum(sapply(1:9, function(i) {
  u <- i + 1
  sapply(u:10, function(j) (14 * (i <= 5) + 24 * (i > 5))*(14 * (j <= 5) + 24 * (j > 5)))
}))

##=====================================
## Distribution S = X1 + ... + X10
n <- 10

lambda.s <- function(alpha_0, remove_i = NULL){
  lam_s <- 5 * sum(vect.lambda) - (n - 1) * alpha_0
  
  if(!is.null(remove_i)){
    to_remove <- vect.lambda[1] - alpha_0 # On mets pour i = 1..5
    if(remove_i >= 6) to_remove <- vect.lambda[2] - alpha_0# on corrige pour i = 6...10
  }
} 

## Trouver v_k tel que P(K1 + ... + Kn = k) = v_k
v_k_fun <- function(m = 2^8, remove_i = NULL,
                    convol = 1) {
  
  ind_15 <- 0   # par défaut, on ne retire pas de i 
  ind_610 <- 0  
  
  ## On vérifie si on doit retirer des i
  if(!is.null(remove_i)){
    if(remove_i <= 5) ind_15 <- 1 # i = 1...5
    if(remove_i >= 6) ind_610 <- 1 # i = 6...10
  }
  
  ## On fait ce qu'on a à faire
  phi_sum_k <- phi.Ki(k1, m)^(5 - ind_15) * phi.Ki(k2, m)^(5 - ind_610)
  Re(fft(phi_sum_k^convol, inverse = T))/(m)
}

## ??? 
k.names <- sapply(1:3, function(p) paste0("k_", p)) 
l_names <- sapply(1:10, function(p) paste0("l_", p))


## Ok
kl_table <- sapply(1:3, function(k) {
  u1 <- c(0.7, 0.2, 0.1) 
  u2 <- c(0.1, 0.4, 0.5) 
  sapply(1:10, function(l) u1[k] * (l <= 5) + u2[k] * (l > 5))
})

## ok, poids
colnames(kl_table) <- k.names


kl_table2 <- as.data.frame(kl_table)
kl_table2$lambda_l <- c(rep(0.1, 5), rep(0.2, 5))

## à revoir
zeta <- function(k, i, convol = 1,
                 maxval = 2^8,
                 full_vec = FALSE){
  
  ## Pour i = 1...5
  zet <- c(k1, rep(0, maxval - length(k1)))
  ## Pour i = 6..10 on corrige
  if(i %in% 6:10) zet <- c(k2, rep(0, maxval - length(k2)))
  
  ## On convolue au besoin
  if(convol > 1){
    zet <- Re(fft(fft(zet)^convol, inverse = TRUE))/maxval
  }
  
  ## Si on veut le vecteur entier, on retourne le vecteur entier
  if(full_vec == TRUE) return(zet)
  
  ## sinon, seulement la valeur demandée
  zet[k + 1]
}

## trouver les thau_k
thau_vec <- function(alpha_0, table, maxval = 2^8,
                     convol = 1, remove_i = NULL) {
  ## En réalité thau_nk calcule la valeur de thau_0 et les 
  ## valeurs à partir de thau_4
  
  v_k <- v_k_fun(remove_i = remove_i, maxval = maxval)          # nu_k
  ls <- lambda.s(alpha_0, remove_i = remove_i)   # lambda_s
  

  
  ## Si on veut retirer un des "i" (\nu^{(-i)})
  if(!is.null(remove_i)){
   table <- table[-1 * remove_i,] # onretire la colonne que l'on ne veut plus
  }
  
  ## zeta_vec (on ajoute beaucoup de colonnes de 0)
  zeta <- cbind(rep(0, nrow(table)),
                table[, -4],
                matrix(rep(0, maxval * nrow(table)),
                       nrow = nrow(table)))
  
  ## On retourne un vecteur complet de tau (peut être très grand)
  full_vec <- sapply(1:maxval, function(k){
    alpha_0/ls * v_k[k + 1] + 
      sum((table$lambda_l - alpha0)/ls * zeta[, k + 1])
  })
  
  ## Si on doit convoluer, on convolue
  if(convol > 1){
    full_vec <- Re(fft(fft(full_vec)^convol, inverse = TRUE))/maxval
  }
  
  return(full_vec)
}

# ## Poids associés à la distribution D de mélange d'erlangs
# vect.thau <- function(alpha_0) {
#   v0 <- thau_nk(alpha_0, 0)
#   v1 <- thau_123(alpha_0, 1)
#   v2 <- thau_123(alpha_0, 2)
#   v3 <- thau_123(alpha_0, 3)
#   
#   vk <- sapply(4:(length(v_k) - 1), function(k) thau_nk(alpha_0, k))
#   
#   c(v0, v1, v2, v3, vk)
# } 

# sum(vect.thau(0.5) > 1)

## Poids associés à la distribution S de mélange d'erlangs
coef.vs <- function(alpha_0) {
  ft <- fft(vect.thau(alpha_0))
  ls <- lambda.s(alpha_0)
  exp.s <- exp(ls * (ft - 1))

  uu <- Re(fft(exp.s, inverse = T)) / (length(exp.s))
  round(uu, 8)[1:100]
}

# sum(coef.vs(0.025))
## Fonction de répartition de S
F.S <- function(x, alpha_0) {
  pk <- coef.vs(alpha_0)
  
  sum(pk[-1] * pgamma(x, 1:(length(pk) - 1), beta)) + pk[1]
}

##=========================== VaR
##===


VaR.S <- function(alpha_0, kappa) {
  fm <- coef.vs(alpha_0)
  if (fm[1] > kappa) return(0)
  
  optimise(function(x) abs(F.S(x, alpha_0) - kappa), c(0, 1000))$minimum
}

VaR.S(0, 0.995)
VaR.S(0.05, 0.995)
VaR.S(0.09, 0.995) 

#seq.alpha_0 <- seq(0, 0.1, by = 0.005)
#seq.kappa <- seq(0, 0.99, by = 0.01)

seq.alpha_0 <- seq(0, 0.1, by = 0.01)
seq.kappa <- c(seq(0, 0.85, by = 0.05), 0.95, 0.99, 0.995, 0.999)

data.var.S <- sapply(seq.alpha_0, function(al) {
  sapply(seq.kappa, function(ka) VaR.S(al, ka))
})

names(seq.alpha_0) <- sapply(seq.alpha_0, function(p) paste0("alpha_", p))
names(seq.kappa) <- sapply(seq.kappa, function(p) paste0("kappa_", p))

colnames(data.var.S) <- names(seq.alpha_0)
rownames(data.var.S) <- names(seq.kappa) 

data.var.S2 <- as.data.frame(data.var.S) 
data.var.S2$kappa <- seq.kappa 

data.var.S2_long <- reshape2::melt(data.var.S2, id.vars = "kappa")

data.var.S2_long$alpha <- as.numeric(str_replace(data.var.S2_long$variable, "alpha_", ""))

col4 <- hcl.colors(8, "Batlow")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark",  "Spectral"

ggsave("graph2_Var.S.png",
       data.var.S2_long %>% 
         ggplot(aes(x = kappa,
                    color = alpha,
                    group = factor(alpha))) + 
         geom_line(aes(y = value), lwd = 2, alpha = 0.8) +
         scale_colour_gradient(name = TeX("Valeur de $\\alpha_0$"), 
                               low = tail(col4, 1),
                               high = col4[1],
                               trans = "exp") + 
         theme_bw() + 
         labs(x = TeX("kappa"),
              y = TeX("$VaR_{\\kappa}(S)$"),
              title = TeX("$VaR_{\\kappa}(S)$ selon la valeur de $\\alpha_0$"),
              subtitle = TeX("n = 10, $\\lambda_i = 0.1 * (i < 6) + 0.2 * (i > 5)$, $\\beta = 0.1$")) + 
         scale_y_continuous(labels = scales::dollar) + 
         scale_x_continuous(labels = scales::dollar)
)
   

## TVaR de S
TVaR.S <- function(alpha_0, kappa) {
  fm <- coef.vs(alpha_0)
  vv <- VaR.S(alpha_0, kappa)
  
  pgamma.bar <- function(k) 1 - pgamma(vv, k + 1, beta)
  
  k.sum <- sapply(1:(length(fm) - 1), function(k) {
    fm[k + 1] * k * pgamma.bar(k) / beta
  })
  
  sum(k.sum) / (1 - kappa)
}

TVaR.S(0, 0.995)
TVaR.S(0.05, 0.995)
TVaR.S(0.09, 0.995)
 
data.var.S2_long$TVaR.S <- sapply(1:nrow(data.var.S2_long), function(k) {
                              TVaR.S(data.var.S2_long[k, 4], data.var.S2_long[k, 1])
                            })

ggsave("graph_TVar.S.png",
       data.var.S2_long %>% 
         ggplot(aes(x = kappa,
                    color = alpha,
                    group = factor(alpha))) + 
         geom_line(aes(y = TVaR.S), lwd = 2, alpha = 0.8) +
         scale_colour_gradient(name = TeX("Valeur de $\\alpha_0$"), 
                               low = tail(col4, 1),
                               high = col4[1],
                               trans = "exp") + 
         theme_bw() + 
         labs(x = TeX("kappa"),
              y = TeX("$TVaR_{\\kappa}(S)$"),
              title = TeX("$TVaR_{\\kappa}(S)$ selon la valeur de $\\alpha_0$"),
              subtitle = TeX("n = 10, $\\lambda_i = 0.1 * (i < 6) + 0.2 * (i > 5)$, $\\beta = 0.1$")) + 
         scale_y_continuous(labels = scales::dollar) + 
         scale_x_continuous(labels = scales::dollar)
)

##============== Calcul des contributions =====================
# kl_table2

## densité de J_(i)
dens.J.i <- function(i, alpha_0, x) {
  lam.i <- kl_table2[i, 4]
  alph <- lam.i - alpha_0
  
  dpois(x, alph)
}

## densité de J_(-i)
dens.J.im <- function(i, alpha_0, x) {
  lam.i <- sum(kl_table2[-i, 4])
  alph <- lam.i - alpha_0
  
  dpois(x, alph)
}
 
## densité de (M_i, N_-i) : pas complet
dens.MN.i <- function(i, ki, ni_moins, alpha_0) {
  to.sum <- sapply(0:min(ki, ni_moins), function(j) {
    d1 <- dens.J.i(i, alpha_0, ki - j)
    d2 <- dens.J.im(i, alpha_0, ni_moins - j)
    
    d1 * d2 * dpois(j, alpha_0)
    
    })
  sum(to.sum)
}

## Coefficients de C_(-i)
v_k_im <- function(im, m = 14) {
  
  if (im <= 5) {
    phi_sum_k <- (phi.Ki(k1, m))^4 * (phi.Ki(k2, m))^5
  } else {
    phi_sum_k <- (phi.Ki(k1, m))^5 * (phi.Ki(k2, m))^4
  }
  
  uu <- Re(fft(phi_sum_k, inverse = T))/(2^m)
  round(uu, 9)[1:100]
}

 

# ## Coeficients thau_im 
#  
# lambda.s.im <- function(alpha_0, im) sum(kl_table2[-im, 4]) - (n - 2) * alpha_0
# 
# thau_nk_im <- function(alpha_0, k, im) {
#   ## En réalité thau_nk calcule la valeur de thau_0 et les 
#   ## valeurs à partir de thau_4
#   
#   (alpha_0 / lambda.s.im(alpha_0, im)) * v_k_im(im)[k + 1]
# }
# 
# thau_nk_im(alpha_0, 9, 1) 
# 
# thau_123_im <- function(alpha_0, k, im) {
#   ls.im <- lambda.s.im(alpha_0, im)
#   
#   #sum1 <- sum(sapply(1:10, function(l) {
#   #  if (l == im) next
#     
#   #  lam.l <- kl_table2[l, 4]
#   #  kl_table2[l, k] * (lam.l - alpha_0) / lambda.s.im(alpha_0, l)
#   #}))
#   
#   ss1 <- 0
#   
#   for (l in 1:10) {
#     if (l == im) next
#     
#     lam.l <- kl_table2[l, 4]
#     ss1 <- ss1 + kl_table2[l, k] * (lam.l - alpha_0) / lambda.s.im(alpha_0, l)
#   }
#   
#   ss1 + (alpha_0 * v_k_im(im)[k + 1]) / ls.im
# }

# ## Poids associés à la distribution D_(-i) de mélange d'erlangs
# vect.thau.im <- function(alpha_0, im) {
#   v0 <- thau_nk_im(alpha_0, 0, im)
#   v1 <- thau_123_im(alpha_0, 1, im)
#   v2 <- thau_123_im(alpha_0, 2, im)
#   v3 <- thau_123_im(alpha_0, 3, im)
#   
#   vk <- sapply(4:(length(v_k_im(im)) - 1), function(k) thau_nk_im(alpha_0, k, im))
#   
#   c(v0, v1, v2, v3, vk)
# } 
# 
# vuu <- vect.thau.im(alpha_0, 6)
# sum(vuu)

## PAr Oli
TVaR_kap_Xi_S <- function(kap, i, alpha0, maxval_k = 3,
                          maxval_gen = 2^8){
  
  VaR <-VaR.S(alpha0, kappa = kap)
  
  terms <- sapply(0:maxal_k, function(ki){ ## sum ki = 0,  ... infty
    
    zet_i_ki <- zeta(k = 3, ## dummy value when full_vec = TRUE
                     i = i, maxval = maxval_gen,
                     convol = ki, full_vec = TRUE)
    
    sum(sapply(0:maxval_k, function(ni){ ## sum n-i = 0 , ... infty
      dens.MN.i(i, ki, ni, alpha_0 = alpha0) * # P(M_i = ki, N_-i = ni)
        sum(sapply(1:maxval_gen, function(k){ # sum k = 1 ... infty
          sum(sapply(1:k, function(l){
            
            v <- v_k_fun(m = maxval_gen)
            tau <- tau_vec(alpha0, kl_table2,
                           maxval = maxval_gen,
                           convol = ni - l,
                           remove_i = i)
            convol_v_tau <- Re(fft(fft(v)*fft(tau),inverse = TRUE))/maxval_gen
            return(
              zet_i_ki[l + 1] * convol_v_tau[k - l + 1] * l/beta *
                pgamma(VaR,
                       k + 1, beta,
                       lower.tail = FALSE))
          }))
        }))
    }))
  })
  
  return(
    1/(1 - kap) * sum(terms)
  )
}
