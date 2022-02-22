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
var.X(vect.lambda[1], k1)
var.X(vect.lambda[2], k2)

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

VaR.Xi <- function(lam, k.i, kappa) {
  fm <- f.Mstar(lam, k.i, m = 14)
  if (fm[1] > kappa) return(0)
  
  optimise(function(x) abs(F.Xi(lam, k.i, x) - kappa), c(0, 200))$minimum
}

VaR.Xi(vect.lambda[1], k1, 0.995)
VaR.Xi(vect.lambda[2], k2, 0.995) # ça ne marche pas

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
  
  Re(fft(phi_sum_k, inverse = T))/(2^m)
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

vect.thau <- function(alpha_0) {
  v0 <- thau_nk(alpha_0, 0)
  v1 <- thau_123(alpha_0, 1)
  v2 <- thau_123(alpha_0, 2)
  v3 <- thau_123(alpha_0, 3)
  
  vk <- sapply(4:(length(v_k) - 1), function(k) thau_nk(alpha_0, k))
  
  c(v0, v1, v2, v3, vk)
} 

sum(vect.thau(0.5) > 1)

coef.vs <- function(alpha_0) {
  ft <- fft(vect.thau(alpha_0))
  ls <- lambda.s(alpha_0)
  exp.s <- exp(ls * (ft - 1))
  
  Re(fft(exp.s, inverse = T)) / (length(exp.s))
}

# coef.vs(0.05)
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
VaR.S(0.1, 0.999) 

seq.alpha_0 <- seq(0, 0.1, by = 0.005)
seq.kappa <- seq(0, 1, by = 0.01)

data.var.S <- sapply(seq.alpha_0, function(al) {
  sapply(seq.kappa, function(ka) VaR.S(al, ka))
})

names(seq.alpha_0) <- sapply(seq.alpha_0, function(p) paste0("alpha_", p))
names(seq.kappa) <- sapply(seq.kappa, function(p) paste0("kappa_", p))

colnames(data.var.S) <- names(seq.alpha_0)
rownames(data.var.S) <- names(seq.kappa)

#data.pi[c(1:6, 37:41), 1:5]
# data.var <- t(data.var)

data.var.S2 <- as.data.frame(data.var.S)
#head(data.var2)
#colnames(data.var2) <- names(seq.alpha)
data.var.S2$kappa <- seq.kappa 

data.var.S2_long <- reshape2::melt(data.var.S2, id.vars = "kappa")

data.var.S2_long$alpha <- as.numeric(str_replace(data.var.S2_long$variable, "alpha_", ""))

col4 <- hcl.colors(8, "Batlow")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark",  "Spectral"

ggsave("graph_Var.S.png",
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

## TVaR de S : Revoir...petit problème dessus
TVaR.S <- sapply(seq.alpha_0, function(al) {
  ps <- coef.vs(al)
  
  sapply(seq.kappa, function(ka) {
     var.k <- data.var.S2[ka, al] 
     
     pgam <- sapply(1:(length(ps) - 1), function(k) 
       (1 - pgamma(var.k, k + 1, beta)) * ps[k + 1] * k / beta)
     
     sum(pgam)
  })
  
})

data.var.S2_long$len <- length(coef.vs(0))

data.var.S2_long$TVaR.S <- (1 - pgamma(data.var.S2_long$value, 2:length(coef.vs(al)), beta)) *
                            ps[k + 1] * k / beta)