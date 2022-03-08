##===================================================================
## EXEMPLE 12
##===================================================================
library(latex2exp)
library(stringr)
library(ggplot2)

library(xtable)

## Initialisation des paramètres

vect.lambda.ex12 <- c(0.003, 0.004) ## Paramètres lambda de la loi Poisson
beta.ex12 <- 1e-3 					## Paramètre beta de la loi Gamma


k1.ex12 <- c(0, 0, 1) 				## Poids associés à i = 1,...,n1	
k2.ex12 <- c(0, 1, 0) 				## poids associés à i = n1+1,..., n1+n2

## ================== Moments de X_i =========================================
## Esperence de Xi
mean.X.Gam <- function(lam, gam) {
   lam * gam * 1000
}
mean1 <- mean.X.Gam(0.003, 2)
mean2 <- mean.X.Gam(0.004, 1)

## Variance de Xi
var.X.Gam <- function(lam, gam) lam * gam * (1 + gam) * 1e6

var1 <- var.X.Gam(0.003, 2)
var2 <- var.X.Gam(0.004, 1)

## =========== VaR de Xi =====================================================
cdf.X.gam <- function(x, lam, gam, kmax = 100) {
  dpois(0, lam) + sum(dpois(1:kmax, lam) * pgamma(x, (1:kmax) * gam, 1e-3))
}

# cdf.X.gam(x = 100, lam = 0.003, gam = 2, kmax = 100)

VaR.X.Gam <- function(kappa, lam, gam) {
  
  optimise(function(x) abs(cdf.X.gam(x, lam, gam) - kappa), c(0, 5000))$minimum
}

VaR11 <- VaR.X.Gam(kappa = 0.995 , lam = 0.003, gam = 2)
VaR12 <- VaR.X.Gam(kappa = 0.995 , lam = 0.004, gam = 1)

VaR21 <- VaR.X.Gam(kappa = 0.9995 , lam = 0.003, gam = 2)
VaR22 <- VaR.X.Gam(kappa = 0.9995 , lam = 0.004, gam = 1)

##=========== TVaR de Xi ===================================================
TVaR.X.Gam <- function(kappa, lam, gam, kmax = 100) {
  varf <- VaR.X.Gam(kappa, lam, gam)
  
  ss1 <- sapply(1:kmax, function(k) {
    dpois(k, lam) * (k * gam * 1e3) * (1 - pgamma(varf, gam * k + 1, 1e-3))
  })
  
  sum(ss1) / (1 - kappa)
}

TVaR11 <- TVaR.X.Gam(kappa = 0.995 , lam = 0.003, gam = 2)
TVaR12 <- TVaR.X.Gam(kappa = 0.995 , lam = 0.004, gam = 1)

TVaR21 <- TVaR.X.Gam(kappa = 0.9995 , lam = 0.003, gam = 2)
TVaR22 <- TVaR.X.Gam(kappa = 0.9995 , lam = 0.004, gam = 1)


## Obtenir un tableau pour latex
data.Xi <- data.frame(i = c("1...n1", "n1 + 1...n1 + n2"),
                      EX_i = c(mean1, mean2),
                      VX_i = c(var1, var2),
                      VaR_Xi1 = c(0, 0),
                      VaR_Xi2 = round(c(VaR21, VaR22), 3),
                      TVaR_Xi1 = round(c(TVaR11, TVaR12), 3),
                      TVaR_Xi2 = round(c(TVaR21, TVaR22), 3))

xtable(data.Xi)

#### ================ S = X1 + ... + X_N ===================================
#### =======================================================================

##============================= VaR de S ================================================================

aa <- 8
phi.Ki <- function(k.i, m = 2^aa) {
 
  ki.long <- c(k.i, rep(0, m - length(k.i)))
  fft(ki.long)
}

phi.Mstar <- function(lam, k.i, m = 2^aa) exp(lam * (phi.Ki(k.i, m) - 1)) 

f.Mstar <- function(lam, k.i, m = 2^aa) {
  Re(fft(phi.Mstar(lam, k.i, m), inverse = T))/(m)
}

lambda.s.ex12 <- function(alpha_0, n1, n2, remove_i = NULL) {
  n <- n1 + n2
  lam_s <- n1 * vect.lambda.ex12[1] + n2 * vect.lambda.ex12[2] - (n - 1) * alpha_0
   
  if(!is.null(remove_i)){
    to_remove <- vect.lambda.ex12[1] - alpha_0 # On mets pour i = 1 ... n1
    if(remove_i > n1) to_remove <- vect.lambda.ex12[2] - alpha_0 # on corrige pour i = n1+1 ... n
    return(lam_s - to_remove)
  }
  return(lam_s)
}


##==================== Trouver v_k tel que P(K1 + ... + Kn = k) = v_k ================================
v_k_fun_ex12 <- function(n1, n2, m = 2^aa, remove_i = NULL,
                    convol = 1) {
  
  if(convol < 1){
    return(
      c(1, rep(0, m - 1))
    )
  }
  
  ind_15 <- 0   # par d?faut, on ne retire pas de i 
  ind_610 <- 0  
  
  ## On v?rifie si on doit retirer des i
  if(!is.null(remove_i)){
    if(remove_i <= n1) ind_15 <- 1 # i = 1...5
    if(remove_i > n1) ind_610 <- 1 # i = 6...10
  }
  
  ## On fait ce qu'on a ? faire
  phi_sum_k <- phi.Ki(k1.ex12, m)^(n1 - ind_15) * phi.Ki(k2.ex12, m)^(n2 - ind_610)
  Re(fft(phi_sum_k^convol, inverse = T))/(m)
}

########################
 
k.names <- sapply(1:2, function(p) paste0("k_", p)) 
#l_names <- sapply(1:10, function(p) paste0("l_", p))


## Ok
kl_table_ex12 <- function(n1, n2) {
 n <- n1 + n2
 tab <- sapply(1:2, function(k) {
  		u1 <- c(0, 1)  
  		u2 <- c(1, 0)
 
  		sapply(1:n, function(l) u1[k] * (l <= n1) + u2[k] * (l > n1))
       })	

  colnames(tab) <- k.names
  tab <- as.data.frame(tab)
  tab$lambda_l <- c(rep(0.003, n1), rep(0.004, n2))
  tab			
}

### ========= zeta ===============================================================
zeta.ex12 <- function(k, i, n1, n2, convol = 1,
                 maxval = 2^aa,
                 full_vec = FALSE){
  n <- n1 + n2
  
  if(convol < 1){
    return(
      c(1, rep(0, maxval - 1))
    )
  }
  
  ## Pour i = 1...n1
  zet <- c(k1.ex12, rep(0, maxval - length(k1.ex12)))
  ## Pour i = n1+1..n on corrige
  if(i %in% (n1 + 1):n) zet <- c(k2.ex12, rep(0, maxval - length(k2.ex12)))
  
  ## On convolue au besoin
  if(convol > 1){
    zet <- Re(fft(fft(zet)^convol, inverse = TRUE))/maxval
  }
  
  ## Si on veut le vecteur entier, on retourne le vecteur entier
  if(full_vec == TRUE) return(zet)
  
  ## sinon, seulement la valeur demand?e
  zet[k + 1]
} 

### =================== trouver les thau_k =============================================
thau_vec_ex12 <- function(alpha_0, n1, n2,  maxval = 2^aa,
                     convol = 1, remove_i = NULL) {

  table <- kl_table_ex12(n1, n2)
 
  if(convol < 1){
    return(
      c(1, rep(0, maxval - 1))
    )
  }
  
  v_k <- v_k_fun_ex12(n1 = n1, n2 = n2, remove_i = remove_i, m = maxval)          # nu_k
  ls <- lambda.s.ex12(alpha_0 = alpha_0, n1 = n1, n2 = n2, remove_i = remove_i)   # lambda_s
  

  
  ## Si on veut retirer un des "i" (\nu^{(-i)})
  if(!is.null(remove_i)){
   table <- table[-1 * remove_i, ] # on retire la colonne que l'on ne veut plus
  }
  
  ## zeta_vec (on ajoute beaucoup de colonnes de 0)
  zeta <- unname(cbind(rep(0, nrow(table)),
                table[, -3],
                matrix(rep(0, (maxval - 3) * nrow(table)),
                       nrow = nrow(table))))
  
  
  
  ## On retourne un vecteur complet de tau (peut ?tre tr?s grand)
  full_vec <- sapply(0:(maxval - 1), function(k){
    ifelse(is.null(remove_i), alpha_0/ls * v_k[k + 1], 0) + 
      sum((table$lambda_l - alpha_0)/ls * zeta[, k + 1])
  })
  
  ## Si on doit convoluer, on convolue
  if(convol > 1){
    full_vec <- Re(fft(fft(full_vec)^convol, inverse = TRUE))/maxval
  }
  
  return(full_vec)
}

## =============== Poids associés à la distribution S de mélange d'erlangs ===================
## ===========================================================================================
coef.vs.ex12 <- function(alpha_0, n1, n2) {
  # table <- kl_table_ex12(n1, n2)
  ft <- fft(thau_vec_ex12(alpha_0, n1 = n1, n2 = n2))
  ls <- lambda.s.ex12(alpha_0, n1 = n1, n2 = n2)
  exp.s <- exp(ls * (ft - 1))

  Re(fft(exp.s, inverse = T)) / (length(exp.s))
  
}

# sum(coef.vs.ex12(alpha_0 = 0.001, n1 = 10, n2 = 10))

##=============== Fonction de répartition de S ================================================
F.S.ex12 <- function(x, alpha_0, n1, n2) {
  pk <- coef.vs.ex12(alpha_0, n1, n2)
   
  sum(pk[-1] * pgamma(x, 1:(length(pk) - 1), beta.ex12)) + pk[1]
}

##============ VaR de S ====================================================================

VaR.S.ex12 <- function(alpha_0, kappa, n1, n2) {
  fm <- coef.vs.ex12(alpha_0, n1, n2)
  if (fm[1] > kappa) return(0)
  
  optimise(function(x) abs(F.S.ex12(x, alpha_0, n1, n2) - kappa), c(0, 100000))$minimum
}

VaR.S.ex12(alpha_0 = 0, kappa = 0.995, n1 = 10, n2 = 10)
VaR.S.ex12(alpha_0 = 0, kappa = 0.995, n1 = 100, n2 = 100)

############ bon 

seq.alpha_0 <- seq(0, 0.003, by = 0.0001)
seq.kappa <- c(seq(0.9, 0.99, by = 0.01), 0.995, 0.999)

col4 <- hcl.colors(8, "Batlow")
## "Viridis", "Plasma, "Purple-Orange"
## "Zissou1", "SunsetDark",  "Spectral"

##  ============== Graphiques =====================================================================
plot.var <- function(n1, n2, seq.alpha_0 = seq.alpha_0, seq.kappa = seq.kappa) {
  n1 = 500; n2 = 500
  data.var.S <- sapply(seq.alpha_0, function(al) {
      sapply(seq.kappa, function(ka) VaR.S.ex12(alpha_0 = al, kappa = ka, n1 = n1, n2 = n2))
   })
  
  names(seq.alpha_0) <- sapply(seq.alpha_0, function(p) paste0("alpha_", p))
  names(seq.kappa) <- sapply(seq.kappa, function(p) paste0("kappa_", p))

  colnames(data.var.S) <- names(seq.alpha_0)
  rownames(data.var.S) <- names(seq.kappa) 

  data.var.S2 <- as.data.frame(data.var.S) 
  data.var.S2$kappa <- seq.kappa 

  data.var.S2_long <- reshape2::melt(data.var.S2, id.vars = "kappa")

  data.var.S2_long$alpha <- as.numeric(str_replace(data.var.S2_long$variable, "alpha_", ""))
  jjj <- sample(1:100, 1)	
  
  ggsave(paste0("graph2_Var.S.ex12_num_", jjj, ".png"),
         
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
              #subtitle = paste0(TeX("$n_1$ = "), n1)) + 
              subtitle = TeX(paste("$n_1 = ", n1, "$", 
                                   "$n_2 = ", n2,"$", 
                                   "$\\gamma_1 = 2$",
                                   "$\\gamma_{n_1 + 1} = 1$", 
                                   "$\\lambda_1 = 0.003$",
                                   "$\\lambda_{n_1 + 1} = 0.004$",
                                   sep = "  "))) +
              #subtitle = paste(TeX(r'($n_1 = $)'), n1) +
         scale_y_continuous(labels = scales::dollar) +
         scale_x_continuous(labels = scales::dollar)
       
       )

}

#plot.var(n1 = 10, n2 = 10)

## ===== TVaR de S =================================================================================
TVaR.S.ex12 <- function(alpha_0, kappa, n1, n2) {
  fm <- coef.vs.ex12(alpha_0, n1, n2)
  vv <- VaR.S.ex12(alpha_0, kappa, n1, n2)
  
  pgamma.bar <- function(k) 1 - pgamma(vv, k + 1, beta.ex12)
  
  k.sum <- sapply(1:(length(fm) - 1), function(k) {
    fm[k + 1] * k * pgamma.bar(k) / beta.ex12
  })
  
  sum(k.sum) / (1 - kappa)
}

tv1 <- TVaR.S.ex12(0, 0.995, 10, 10)
tv2 <- TVaR.S.ex12(0, 0.995, 100, 100)
tv3 <- TVaR.S.ex12(0, 0.995, 500, 500)

tv4 <- TVaR.S.ex12(0.001, 0.995, 10, 10)
tv5 <- TVaR.S.ex12(0.001, 0.995, 100, 100)
tv6 <- TVaR.S.ex12(0.001, 0.995, 500, 500)

vr1 <- VaR.S.ex12(0, 0.995, 10, 10)
vr2 <- VaR.S.ex12(0, 0.995, 100, 100)
vr3 <- VaR.S.ex12(0, 0.995, 500, 500)

vr4 <- VaR.S.ex12(0.001, 0.995, 10, 10)
vr5 <- VaR.S.ex12(0.001, 0.995, 100, 100)
vr6 <- VaR.S.ex12(0.001, 0.995, 500, 500)

#TVaR.S.ex12(0.002, 0.995, 10, 10)
#TVaR.S.ex12(0.002, 0.995, 100, 100)
#TVaR.S.ex12(0.002, 0.995, 500, 500)

## Obtenir un tableau pour latex
data.ex12 <- data.frame(alpha_0 = c(rep(0, 3), rep(0.001, 3)),
                      n1 = rep(c(10, 100, 500), 2),
                      n2 = rep(c(10, 100, 500), 2),
                      VaR = round(c(vr1, vr2, vr3, vr4, vr5, vr6), 3),
                      TVaR  = round(c(tv1, tv2, tv3, tv4, tv5, tv6), 3),
                      TVaR_X1.S = round(c(contri.11, contri.21, contri.31, contri.41, contri.51, contri.61), 3),
                      TVaR_X2.S = round(c(contri.12, contri.22, contri.32, contri.42, contri.52, contri.62), 3))

xtable(data.ex12)

## ================= Contributions ================================================================
## ================================================================================================

TVaR_kap_Xi_S_ex12 <- function(kap, i, alpha0, n1, n2, maxval_k = 50,
                          maxval_gen = 2^8){
  
  kl_table2 <- kl_table_ex12(n1, n2)
  VaR <- VaR.S.ex12(alpha0, kappa = kap, n1 = n1, n2 = n2)
  
  
  terms <- sapply(0:maxval_k, function(ki){ ## sum ki = 0,  ... infty
    
    zet_i_ki <- zeta.ex12(k = 3, ## dummy value when full_vec = TRUE
                     i = i, n1 = n1, n2 = n2, maxval = maxval_gen,
                     convol = ki, full_vec = TRUE)
    
    sum(sapply(0:maxval_k, function(ni){ ## sum n-i = 0 , ... infty
      value <- sum(sapply(0:min(ki, ni), function(j){ ## sum j = 0... min(ki, ni)
        
        ## on calcule nos poids pour les deux autres boucles (ne dépendent pas de k ni de l)
        v <- v_k_fun_ex12(n1 = n1, n2 = n2, m = maxval_gen, convol = j, remove_i = i)
        tau <- thau_vec_ex12(alpha_0 = alpha0, 
                             n1 = n1, n2 = n2,
                             maxval = maxval_gen,
                             convol = ni - j,
                             remove_i = i)
        convol_v_tau <- Re(fft(fft(v)*fft(tau),inverse = TRUE))/maxval_gen
        
        
        ## densité des poissons
        dpois(ki - j, kl_table2$lambda_l[i] - alpha0) *
          dpois(ni - j, lambda.s.ex12(alpha0, n1 = n1, n2 = n2, remove_i = i)) * 
          dpois(j, alpha0) * 
          
          sum(sapply(1:(maxval_gen - 1), function(k){ # sum k = 1 ... infty
            
            sum(sapply(1:k, function(l){
              return(
                zet_i_ki[l + 1] * convol_v_tau[k - l + 1] * l/beta.ex12 *
                  pgamma(VaR,
                         k + 1, beta.ex12,
                         lower.tail = FALSE))
            }))
            
          }))
      }))
      
      to_print <- (maxval_k*ki + ni)/((maxval_k + 1)^2)
      print(paste0("progression...", scales::percent(to_print, 0.1)))
      
      return(value)
    }))
  })
  
  return(
    1/(1 - kap) * sum(terms)
  )
}

## alpha_0  = 0
contri.11 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0, n1 = 10, n2 = 10, maxval_k = 25)
contri.12 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 11, alpha0 = 0, n1 = 10, n2 = 10, maxval_k = 40)

contri.21 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0, n1 = 100, n2 = 100, maxval_k = 25)
contri.22 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 101, alpha0 = 0, n1 = 100, n2 = 100, maxval_k = 40)

contri.31 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0, n1 = 500, n2 = 500, maxval_k = 25)
contri.32 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 501, alpha0 = 0, n1 = 500, n2 = 500, maxval_k = 40)

## alpha_0 = 0.001
contri.41 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.001, n1 = 10, n2 = 10, maxval_k = 25)
contri.42 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 11, alpha0 = 0.001, n1 = 10, n2 = 10, maxval_k = 40)

contri.51 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.001, n1 = 100, n2 = 100, maxval_k = 25)
contri.52 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 101, alpha0 = 0.001, n1 = 100, n2 = 100, maxval_k = 40)

contri.61 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.001, n1 = 500, n2 = 500, maxval_k = 25)
contri.62 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 501, alpha0 = 0.001, n1 = 500, n2 = 500, maxval_k = 40)

## alpha_0 = 0.002
contri.71 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.002, n1 = 10, n2 = 10, maxval_k = 25)
contri.72 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 11, alpha0 = 0.002, n1 = 10, n2 = 10, maxval_k = 40)

contri.81 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.002, n1 = 100, n2 = 100, maxval_k = 25)
contri.82 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 101, alpha0 = 0.002, n1 = 100, n2 = 100, maxval_k = 40)

contri.91 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 1, alpha0 = 0.002, n1 = 500, n2 = 500, maxval_k = 25)
contri.92 <- TVaR_kap_Xi_S_ex12(kap = 0.995, i = 501, alpha0 = 0.002, n1 = 500, n2 = 500, maxval_k = 40)
