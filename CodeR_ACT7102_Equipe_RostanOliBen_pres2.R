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
  uu <- Re(fft(phi.Mstar(lam, k.i, m), inverse = T))/(2^m)
  
  round(uu, 8)[1:100]
}

## Espérance de Xi
Esp.X <- function(lam, k.i) {
  fm <- f.Mstar(lam, k.i)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) / beta)
}

## Espérance de Xi^2  
Esp.X2 <- function(lam, k.i) {
  fm <- f.Mstar(lam, k.i, m = 6)
  len <- length(fm)
  sum( fm[-1] * (1:(len - 1)) * (2 : len) / beta)
}

## Variance de Xi
var.X <- function(lam, k.i) Esp.X2(lam, k.i) - (Esp.X(lam, k.i))^2

## Applications
Esp.X(vect.lambda[1], k1)
Esp.X(vect.lambda[2], k2)

## Je ne comprends pas pourquoi ça ne donne pas 
#var.X(vect.lambda[1], k1)
#var.X(vect.lambda[2], k2)

## Moment d'ordre k de Mstar
Moments <- function(lam, k.i, ord) {
  fm <- f.Mstar(lam, k.i, m = 14)
  len <- length(fm)
  
  sum(fm * (0 : (len - 1))^ord)
}

## Variance de Xi : ça marche :)
var.Xi <- function(lam, k.i) { 
  mum <- Moments(lam, k.i, ord = 1) + Moments(lam, k.i, ord = 2) - (Moments(lam, k.i, ord = 1))^2
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
F.Xi <- function(lam, k.i, x) {
  fm <- f.Mstar(lam, k.i, m = 14)
  len <- length(fm)
  
  sum(fm[-1] * pgamma(x, 1 : (len - 1), beta)) + fm[1]
}

## VaR de Xi
VaR.Xi <- function(lam, k.i, kappa) {
  fm <- f.Mstar(lam, k.i, m = 14)
  if (fm[1] > kappa) return(0)
  
  optimise(function(x) abs(F.Xi(lam, k.i, x) - kappa), c(0, 200))$minimum
}

VaR.Xi(vect.lambda[1], k1, 0.995)
VaR.Xi(vect.lambda[2], k2, 0.995) # ça ne marche pas

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

lambda.s <- function(alpha_0) 5 * sum(vect.lambda) - (n - 1) * alpha_0

## Trouver v_k tel que P(K1 + ... + Kn = k) = v_k
 
v_k_fun <- function(m = 14) {
  phi_sum_k <- (phi.Ki(k1, m) * phi.Ki(k2, m))^5
  
  uu <- Re(fft(phi_sum_k, inverse = T))/(2^m)
  round(uu, 8)[1:100]
}
v_k <- v_k_fun()

k.names <- sapply(1:3, function(p) paste0("k_", p)) 
l_names <- sapply(1:10, function(p) paste0("l_", p))

kl_table <- sapply(1:3, function(k) {
  u1 <- c(0.7, 0.2, 0.1) 
  u2 <- c(0.1, 0.4, 0.5) 
  sapply(1:10, function(l) u1[k] * (l <= 5) + u2[k] * (l > 5))
})

colnames(kl_table) <- k.names

kl_table2 <- as.data.frame(kl_table)
kl_table2$lambda_l <- c(rep(0.1, 5), rep(0.2, 5))

## trouver les thau_k
thau_nk <- function(alpha_0, k) {
  ## En réalité thau_nk calcule la valeur de thau_0 et les 
  ## valeurs à partir de thau_4
  
  (alpha_0 / lambda.s(alpha_0)) * v_k[k + 1]
}

thau_123 <- function(alpha_0, k) {
  ls <- lambda.s(alpha_0)
  
  sum1 <- sum(sapply(1:10, function(l) {
                    lam.l <- kl_table2[l, 4]
                    kl_table2[l, k] * (lam.l - alpha_0) / ls
          }))
  
  sum1 + (alpha_0 * v_k[k + 1]) / ls
}

## Poids associés à la distribution D de mélange d'erlangs
vect.thau <- function(alpha_0) {
  v0 <- thau_nk(alpha_0, 0)
  v1 <- thau_123(alpha_0, 1)
  v2 <- thau_123(alpha_0, 2)
  v3 <- thau_123(alpha_0, 3)
  
  vk <- sapply(4:(length(v_k) - 1), function(k) thau_nk(alpha_0, k))
  
  c(v0, v1, v2, v3, vk)
} 

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
 
## densité de (M_i, N_-i) :  
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

 

## Coeficients thau_im 
 
lambda.s.im <- function(alpha_0, im) sum(kl_table2[-im, 4]) - (n - 2) * alpha_0

thau_nk_im <- function(alpha_0, k, im) {
  ## En réalité thau_nk calcule la valeur de thau_0 et les 
  ## valeurs à partir de thau_4
  
  (alpha_0 / lambda.s.im(alpha_0, im)) * v_k_im(im)[k + 1]
}

thau_nk_im(alpha_0, 9, 1) 

thau_123_im <- function(alpha_0, k, im) {
  ls.im <- lambda.s.im(alpha_0, im)
  
  #sum1 <- sum(sapply(1:10, function(l) {
  #  if (l == im) next
    
  #  lam.l <- kl_table2[l, 4]
  #  kl_table2[l, k] * (lam.l - alpha_0) / lambda.s.im(alpha_0, l)
  #}))
  
  ss1 <- 0
  
  for (l in 1:10) {
    if (l == im) next
    
    lam.l <- kl_table2[l, 4]
    ss1 <- ss1 + kl_table2[l, k] * (lam.l - alpha_0) / lambda.s.im(alpha_0, l)
  }
  
  ss1 + (alpha_0 * v_k_im(im)[k + 1]) / ls.im
}

## Poids associés à la distribution D_(-i) de mélange d'erlangs
vect.thau.im <- function(alpha_0, im) {
  v0 <- thau_nk_im(alpha_0, 0, im)
  v1 <- thau_123_im(alpha_0, 1, im)
  v2 <- thau_123_im(alpha_0, 2, im)
  v3 <- thau_123_im(alpha_0, 3, im)
  
  vk <- sapply(4:(length(v_k_im(im)) - 1), function(k) thau_nk_im(alpha_0, k, im))
  
  uu <- c(v0, v1, v2, v3, vk)
  uu / sum(uu)
} 

#vuu <- vect.thau.im(alpha_0, 6)
#sum(vuu)

## Coder la formule du Box II
#weight.Bi(ki, i)[l + 1] weight.DC

v_thau <- function(alpha_0, i, l, m, ni) {
  
  to.sum <- sapply(0:m, function(u) {
                    vu <- weight.Ci(i, l)[u + 1]
                    thau <- weight.Di(alpha_0, i, ni, l)[m - u + 1]
    
                    vu * thau
              })
  
  sum(to.sum)
}

weight.box2 <- function(alpha_0, kappa, ki, i, ni) {
  
  w1 <- weight.Bi(ki, i)
  len <- length(w1)
  vect.k <- 1:(len - 1)
  ## data.var.S2_long[kl, 3] = Valeur de VaR.S dans le data.frame data.var.S2_long
  
  wl <- which(data.var.S2_long$kappa == kappa & data.var.S2_long$alpha == alpha_0)
  
  sk <- sapply(vect.k, function(k) {
    vect.l <- 1:k
    pbar <- 1 - pgamma(data.var.S2_long[wl, 3], k + 1, beta)
    
    sl <- sapply(vect.l, function(l) {
      ww1 <- w1[l + 1]
      ww2 <- v_thau(alpha_0, i, l, k - l, ni)
      
      ww1 * ww2 * (l / beta) * pbar 
      
    })
    
    sum(sl)
  })
  
  sum(sk)
}

##
TVaR.box2 <- function(i, alpha_0, kappa, kimax = 3, nimax = 3) {
  vect.ki <- 1:kimax
  vect.ni <- 1:nimax
  
  s2 <- sapply(vect.ki, function(ki) {
    
    s1 <- sapply(vect.ni, function(ni) {
      
      dd1 <- dens.MN.i(i, ki, ni, alpha_0)
      dd2 <- weight.box2(alpha_0, kappa, ki, i, ni)
      
      dd1 * dd2
    })
    
    sum(s1)
    
  })
  
  sum(s2) / (1 - kappa)
}

TVaR.box2(1, alpha_0, kappa, kimax = 2, nimax = 2)

## Poids associé à  sum_{m = 1}^{ni - j} D_{i, m} + sum_{r = 1}^{j} C_{i, r}
weight.DC <- function(alpha_0, im, ni, j) { 
  p2 <- weight.Di(alpha_0, im, ni, j) 
  p3 <- weight.Ci(im, j)
  
  phip <-  fft(p2) * fft(p3)
  
  uu <- Re(fft(phip, inverse = T)) / length(phip)
  round(uu, 8)[1:100]
}


## Poids associés à sum_{l = 1}^ki B_{i,l}
weight.Bi <- function(ki, i) {
  
  if(i <= 5) {
    kk <- k1
  } else {
    kk <- k2
  }
  phi.kk <- (phi.Ki(kk, m = 8))^ki
  
  Re(fft(phi.kk, inverse = T)) / length(phi.kk)
}

#sum(weight.Bi(6, 3))

## Poids associés à sum_{m = 1}^{ni - j} D_{i, m}
weight.Di <- function(alpha_0, im, ni, j) {
  kk <- vect.thau.im(alpha_0, im)
  
  phi.kk <- (phi.Ki(kk, m = 14))^(ni - j)
  
  uu <- Re(fft(phi.kk, inverse = T)) / length(phi.kk)
  round(uu, 8)[1:100]
}

# length(weight.Di(alpha_0, 6, 3, 1))

## Poids associés à sum_{r = 1}^{j} C_{i, r}
weight.Ci <- function(im, j) {
  kk <- v_k_im(im)
  
  phi.kk <- (phi.Ki(kk, m = 14))^j
  
  uu <- Re(fft(phi.kk, inverse = T)) / length(phi.kk)
  round(uu, 8)[1:100]
}

#length(weight.Ci(4, 2))

## Poids associé à sum_{l = 1}^ki B_{i,l} + sum_{m = 1}^{ni - j} D_{i, m} + sum_{r = 1}^{j} C_{i, r}
weight.BDC <- function() {
  p1 <- weight.Bi(ki, i) 
  p2 <- weight.Di(alpha_0, im, ni, j) 
  p3 <- weight.Ci(im, j)
  
  phip <- fft(p1) * fft(p2) * fft(p3)
  
  uu <- Re(fft(phip, inverse = T)) / length(phip)
  round(uu, 8)[1:100]
}

## Proposition 8
qm10 <- function(vect.m, vect.lam, alpha_0) {
  
  sp1 <- sapply(0:min(vect.m), function(j) {
          dp1 <- dpois(j, alpha_0)
    
          pr1 <- sapply(1:length(vect.m), function(i) {
      
                    dpois(vect.m[i] - j, vect.lam[i] - alpha_0)
                })
    
          dp1 * prod(pr1)
  })
  
  sum(sp1)
}
 

psi.ji <- function(vect.j, vect.m, mi, i, m.max = 3) {
  pv <- weight.Bi(mi, i)
  
  sapply(1:length(vect.m), function(m.i) {
    valeurs <- m_values[m.i, ]
    
  })
  
} 

m1 <- 1:6
m2 <- 1:6
m3 <- 1:6

m_values <- expand.grid(m1, m2, m3)

m_values[33,]
##===================================================================
## EXEMPLE 12
##===================================================================

## Espérence de Xi
mean.X.Gam <- function(lam, gam) {
   lam * gam * 1000
}
mean.X.Gam(0.003, 2)
mean.X.Gam(0.004, 1)

## Variance de Xi
var.X.Gam <- function(lam, gam) lam * gam * (1 + gam) * 1e6

var.X.Gam(0.003, 2)
var.X.Gam(0.004, 1)

## VaR de Xi
cdf.X.gam <- function(x, lam, gam, kmax = 100) {
  dpois(0, lam) + sum(dpois(1:kmax, lam) * pgamma(x, (1:kmax) * gam, 1e-3))
}

# cdf.X.gam(x = 100, lam = 0.003, gam = 2, kmax = 100)

VaR.X.Gam <- function(kappa, lam, gam) {
  
  optimise(function(x) abs(cdf.X.gam(x, lam, gam) - kappa), c(0, 5000))$minimum
}

VaR.X.Gam(kappa = 0.9995 , lam = 0.003, gam = 2)
VaR.X.Gam(kappa = 0.9995 , lam = 0.004, gam = 1)

## TVaR de Xi 
TVaR.X.Gam <- function(kappa, lam, gam, kmax = 100) {
  varf <- VaR.X.Gam(kappa, lam, gam)
  
  ss1 <- sapply(1:kmax, function(k) {
    dpois(k, lam) * (k * gam * 1e3) * (1 - pgamma(varf, gam * k + 1, 1e-3))
  })
  
  sum(ss1) / (1 - kappa)
}

TVaR.X.Gam(kappa = 0.9995 , lam = 0.003, gam = 2)
TVaR.X.Gam(kappa = 0.9995 , lam = 0.004, gam = 1)

## VaR de S
m1 <- 1:2

list.m = list()
for (i in 1:10) {
  list.m[[i]] <- m1
}

m_values <- expand.grid(list.m)
m_values$qmn <- sapply(1:nrow(m_values), function(j) {
  qm10(vect.m = as.numeric(m_values[j, ]), vect.lam, alpha_0)
})
   

cdf_S <- function(x, n1, n2, alpha_0) {
  vect.lam <- c(rep(0.003, n1), rep(0.004, n2)) ## paramètres de la Poisson
  vect.gam <- c(rep(2, n1), rep(1, n2)) ## Bi ~ Gamma(vect.gam[i], 1e-3)
  
  n <- n1 + n2
  nn1 <- 0:1
  
  list.m = list()
  for (i in 1:n) {
    list.m[[i]] <- nn1
  }
  
  m_values <- expand.grid(list.m)
  m_values$qmn <- sapply(1:nrow(m_values), function(j) {
    qm10(vect.m = as.numeric(m_values[j, ])[1:n], vect.lam, alpha_0)
  })
  
  m_values$mi.al <- sapply(1:nrow(m_values), function(j) {
    sum(as.numeric(m_values[j, ])[1:n] * vect.gam)
  })
  
  m_values$pgam <- sapply(1:nrow(m_values), function(j) {
    pgamma(x, m_values$mi.al, 1e-3)
  })
  
  # q0 <- exp((n - 1) * alpha_0 - sum(vect.lam))
  sum(m_values$qmn[-1] * m_values$pgam[-1]) + m_values$qmn[1]
}

# cdf_S(4, 4, 5, alpha_0 = 0.001)

VaR.S.Gam <- function(kappa, n1, n2, alpha_0) {
  
  optimise(function(x) abs(cdf_S(x, n1, n2, alpha_0) - kappa), c(0, 10000))$minimum
}

memory.limit(size = 35000) 
VaR.S.Gam(0.995, 10, 10, 0)
