#Master thesis (Final)
library(tidyverse)
library(deSolve)
source("M_Theme setting.R", encoding = 'CP950', echo = T)
# Function setting----
{
  Jacobian_full = function(parms, E) {
    with(c(parms, E), {
      matrix(data =
               c(
                 -c1*b1*P1-c2*b2*P2-d, h1*o1, h2*o2, -c1*b1*H, -c2*b2*H, 0,
                 b1*P1, -(m1+o1), 0, b1*H, 0, 0,
                 b2*P2, 0, -(m2+o2), 0, b2*H, 0,
                 -b1*P1, e1H*psi1*a1*S, 0, e1*a1*S - b1*H - m1, 0, e1*a1*P1 + e1H*psi1*a1*P1H,
                 -b2*P2, 0, e2H*psi2*a2*S, 0, e2*a2*S - b2*H - m2, e2*a2*P2 + e2H*psi2*a2*P2H,
                 0, -psi1*a1*S, -psi2*a2*S, -a1*S, -a2*S, r*(1-S/K) - a1*P1 - psi1*a1*P1H - a2*P2 - psi2*a2*P2H - S*r/K
               ), nrow = 6, byrow = TRUE)
    })
  } ##for full model (S, P1, P2, P1H, P2H, and H)
  
  Jacobian_mts = function(parms, E) {
    with(c(parms, E), {
      matrix(data =
               c(
                 e1*a1*S - m1, 0, e1*a1*P1,
                 0, e2*a2*S - m2, e2*a2*P2,
                 -a1*S, -a2*S, r*(1-(S/K)) - a1*P1 - a2*P2 - (S*r)/K
               ), nrow = 3, byrow = TRUE)
    })
  } ##for multi strain w/o H (S, P1, P2)
  
  Jacobain_sgs = function(parms, E) {
    with(c(parms, E), {
      J11 = -c1*b1*P1-d
      J12 = h1*o1 
      J13 = -c1*b1*H
      J14 = 0
      J21 = b1*P1
      J22 = -(m1+o1)
      J23 = b1*H
      J24 = 0
      J31 = -b1*P1 
      J32 = e1H*psi1*a1*S
      J33 = e1*a1*S-b1*H-m1
      J34 = e1*a1*P1+e1H*psi1*a1*P1H
      J41 = 0
      J42 = -psi1*a1*S
      J43 = -a1*S
      J44 = r*(1-S/K) - a1*P1 - psi1*a1*P1H - S*r/K
      matrix(data = 
               c(J11, J12, J13, J14,
                 J21, J22, J23, J24,
                 J31, J32, J33, J34,
                 J41, J42, J43, J44), 
             byrow = T, nrow = 4, ncol = 4)
    })} ##for single strain (S, P1, P1H, and H)
  
  Jacobain_sgsN = function(parms, E){
    with(c(parms, E), {
      J11 = e1*a1*S-m1
      J12 = e1*a1*P1
      J21 = -a1*S
      J22 = r*(1-S/K) - a1*P1 - S*r/K
      matrix(data = 
               c(J11, J12, 
                 J21, J22), 
             byrow = T, nrow = 2, ncol = 2)
    })
  } ##for single strain (S, P1)
  
  Eigen = function(J) {
    if (any(is.na(J)) || any(is.infinite(J))){
      return(NA)
    }else{
      return(max(Re(eigen(J)$values)))
    }
  }
  
  f_E_C = function(parms){
    with(parms, {
      D1 = (m1+o1)
      D2 = (m2+o2)
      A = D1*b1*e2H*psi2*a2*b2 - D2*b2*e1H*psi1*a1*b1
      B = D1*b1*e2*a2*D2 + m1*D1*e2H*psi2*a2*b2 - D2*b2*e1*a1*D1 - m2*D2*e1H*psi1*a1*b1
      C = m1*e2*a2*D1*D2 - m2*e1*a1*D2*D1
      E = B^2 - 4*A*C
      if(E < 0 || is.na(E)){
        H = NA
      }else{
        H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
        H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
        H_12 = c(H1,H2)
        H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
        if(length(H_12) == 0){
          H = NA
        } else {
          H = min(H_12) #Choose the smallest one which is more biologically meaningful
        }
      }
      A1 = (b1 * H) / (m1 + o1)
      A2 = (b2 * H) / (m2 + o2)
      B1 = h1 * o1 * A1 - c1 * b1 * H
      B2 = h2 * o2 * A2 - c2 * b2 * H
      D1 = (1 + psi1 * A1) * a1
      D2 = (1 + psi2 * A2) * a2
      S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      #S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
      P1 = (d * H) / B1 - (B2 / B1) * (((B1 * r * (1 - S / K)) - D1 * d * H) / (D2 * B1 - D1 * B2))
      P2 = (B1 * r * (1 - S / K) - D1 * d * H) / (D2 * B1 - D1 * B2)
      P1H = A1 * P1
      P2H = A2 * P2
      
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_S = function(parms){
    with(parms,{
      H = 0
      P1H = 0
      P2H = 0
      P1 = 0
      P2 = 0
      S = K
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P1 = function(parms){
    with(parms, {
      H = 0
      P1H = 0
      P2H = 0
      S = m1/(e1*a1)
      P1 = (r/a1)*(1-(S/K))
      P2 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P2 = function(parms){
    with(parms, {
      H = 0
      P1H = 0
      P2H = 0
      P1 = 0
      S = m2/(e2*a2)
      P2 = (r/a2)*(1-(S/K))
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P1H = function(parms){
    with(parms,{
      P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
      A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
      B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
      C = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
      E = B^2 - 4*A*C
      if(E < 0 || is.na(E)){
        H = NA
      }else{
        H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
        H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
        H_12 = c(H1,H2)
        H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
        if(length(H_12) == 0){
          H = NA
        } else {
          H = min(H_12) #Choose the smallest one which is more biologically meaningful
        }
      }
      S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      P1H = P1 * (b1 * H) / (m1 + o1)
      P2H = 0
      P2 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P2H = function(parms){
    with(parms, {
      P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
      A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
      B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
      C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
      E = B^2 - 4*A*C
      if(E < 0 || is.na(E)){
        H = NA
      }else{
        H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
        H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
        H_12 = c(H1,H2)
        H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
        if(length(H_12) == 0){
          H = NA
        } else {
          H = min(H_12) #Choose the smallest one which is more biologically meaningful
        }
      }
      S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
      P2H = P2 * (b2 * H) / (m2 + o2)
      P1H = 0
      P1 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
    
  }
}

#Result 1: bifurcation plot of W/ and w/o hyperparasite ----

parms = list(
  r = 1, K = 10,
  a1 = 0.25, psi1 = 1, e1 = 0.5, 
  b1 = 0.2, e1H = 0.5, 
  o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.01)

comp_out = expand.grid(m1 = seq(0.001, 0.2, by = 0.001))

## Single strain model with no H (i.e., S and P) ----
comp_out_SP = as.data.frame(cbind(comp_out,
                               matrix(0, 
                                      nrow = dim(comp_out)[1],
                                      ###dim(data)[1] is the number of row of data; [2] is col
                                      ncol = 2+1)))
names(comp_out_SP) = c("m1", "P1", "S", "Stability")


for (i in 1:dim(comp_out_SP)[1]) {
  temp_parms = parms
  temp_parms["m1"]  = comp_out_SP[i, "m1"]
  E_P1 = f_E_P1(temp_parms)
  E_S = f_E_S(temp_parms)
  Lambda_E_P1 = Eigen(Jacobain_sgsN(temp_parms, E_P1))
  Lambda_E_S = Eigen(Jacobain_sgsN(temp_parms, E_S))
  
  Stable_E = c()
  
  if(!is.na(Lambda_E_P1) && is.finite(Lambda_E_P1) && Lambda_E_P1 < 0 && all(E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P1")
    comp_out_SP[i, "Stability"] = "Stable"
    comp_out_SP[i, c("P1", "S")] = f_E_P1(temp_parms)[c("P1", "S")]
  }
  if(!is.na(Lambda_E_S) && is.finite(Lambda_E_S) && Lambda_E_S < 0 && all(E_S[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "S")
    comp_out_SP[i, "Stability"] = "Stable"
    comp_out_SP[i, c("P1", "S")] = f_E_P1(temp_parms)[c("P1", "S")]
  }
  
  comp_out_SP[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  
  if(length(Stable_E) == 0){
    comp_out_SP[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out_SP[i, "Stability"] = "ASS"
  }
}

comp_out_SP = comp_out_SP[,c(1:3)]
names(comp_out_SP) = c("m1", "SP_P", "SP_S")

## Single strain model with H (i.e., S, P and H) ----
comp_out_SPH = as.data.frame(cbind(comp_out,
                               matrix(0, 
                                      nrow = dim(comp_out)[1],
                                      ###dim(data)[1] is the number of row of data; [2] is col
                                      ncol = 4)))
names(comp_out_SPH) = c("m1", "H", "P1H", "P1", "S")

comp_out_SPH

for (i in 1:dim(comp_out_SPH)[1]) {
  temp_parms = parms
  temp_parms["m1"]  = comp_out_SPH[i, "m1"]
  
  E_P1H = f_E_P1H(temp_parms)
  Lambda_P1H = Eigen(Jacobain_sgs(temp_parms, E_P1H))
  E_P1 = f_E_P1(temp_parms)
  Lambda_P1 = Eigen(Jacobain_sgs(temp_parms, E_P1))
  E_S = f_E_S(temp_parms)
  Lambda_S = Eigen(Jacobain_sgs(temp_parms, E_S))
  
  Stable_E = c()
  if (!is.na(Lambda_P1H) == T && is.finite(Lambda_P1H) == T && Lambda_P1H < 0 && all(!is.na(E_P1H) == T)){
    Stable_E = c(Stable_E, "P1H")
    comp_out_SPH[i, "Stability"] = "Stable"
    comp_out_SPH[i, c("H", "P1H", "P1", "S")] = E_P1H[c("H", "P1H", "P1", "S")]
  }
  if(!is.na(Lambda_P1) == T && is.finite(Lambda_P1) == T && Lambda_P1 < 0 && all(!is.na(E_P1) == T)){
    Stable_E = c(Stable_E, "P1")
    comp_out_SPH[i, "Stability"] = "Stable"
    comp_out_SPH[i, c("H", "P1H", "P1", "S")] = E_P1[c("H", "P1H", "P1", "S")]
  }
  if(!is.na(Lambda_S) == T && is.finite(Lambda_S) == T && Lambda_S < 0 && all(!is.na(E_S) == T)){
    Stable_E = c(Stable_E, "S")
    comp_out_SPH[i, "Stability"] = "Stable"
    comp_out_SPH[i, c("H", "P1H", "P1", "S")] = E_S[c("H", "P1H", "P1", "S")]
  }
  
  comp_out_SPH[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  if(length(Stable_E) == 0){
    comp_out_SPH[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out_SPH[i, "Stability"] = "ASS"
  }
}

###It is possible that data contains NA because m1 can be large enough that P1 would no longer persist H.
comp_out_SPH = mutate(comp_out_SPH, P.total = P1+P1H)
comp_out_SPH[,c(1:4, 8, 5)] -> comp_out_SPH

names(comp_out_SPH) = c("m1", "SPH_H", "SPH_PH", "SPH_P", "SPH_P.total", "SPH_S")


#comp_out = left_join(comp_out_SP, comp_out_SPH, by = "m1")
comp_out = cbind(comp_out_SP, comp_out_SPH)[,-4]

## Bifurcation plot ----
D = 
  comp_out %>%
  #select(c(m1, SP_S, SPH_S)) %>% #P1, P2, P1H, P2H, H, S
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1))
#gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
D$Group = 0
D[D[, "Species"] == "SP_S", ]$Group = 1
D[D[, "Species"] == "SPH_H", ]$Group = 1
D[D[, "Species"] == "SPH_S", ]$Group = 1

D[D[, "Species"] == "SP_P", ]$Group = 2
D[D[, "Species"] == "SPH_P", ]$Group = 2
D[D[, "Species"] == "SPH_PH", ]$Group = 2
D[D[, "Species"] == "SPH_P.total", ]$Group = 2


#D$Species = factor(D$Species, levels = c("SPH_S", "SP_S"))
ggplot(D, aes(x = m1, y = Abundance, color = Species)) +
  geom_line(mapping = aes(x = m1, y = Abundance, color = Species), lwd = 1) +
  labs(x = "Intrinsic mortality rate of pathogen (m)", y = "Abundance", color = "Species")+
  geom_vline(xintercept = 0.0755, color = "black", linetype = 2, linewidth = 0.7) +
  scale_colour_manual(labels = c("SP_S" = expression(S^"*" ~ "|" [E[P]]),
                                 "SP_P" = expression(P^"*" ~ "|" [E[P]]),
                                 "SPH_S" = expression(S^"*" ~ "|" [E[PH]]),
                                 "SPH_P" = expression(P^"*" ~ "|" [E[PH]]),
                                 "SPH_PH" = expression(P[H]^"*" ~ "|" [E[PH]]),
                                 "SPH_P.total" = expression(P.total^"*" ~ "|" [E[PH]]),
                                 "SPH_H" = expression(H^"*" ~ "|" [E[PH]])),
                      
                       values = c("SP_S" = alpha("#00AF66", 0.1),
                                  "SP_P" = alpha("#a50f15", 0.1),
                                  "SPH_S" = "#00AF66",
                                  "SPH_P" = alpha("#a50f15", 0.4),
                                  "SPH_PH" = "#a50f15",
                                  "SPH_P.total" = "grey50",
                                  "SPH_H" = "#D6B701"))+
  theme(axis.title.y.right = element_text(angle = 90), legend.position = "bottom", strip.background = element_blank(), strip.text = element_blank())


ggsave("Result_1 SS Bifurcation plot.png", width = 15, height = 12, units = "cm", dpi = 800) #15,10


#Result 2+3: parameter space of m1 and m2 W/ and w/o hyperparasite ----
parms <- list(r = 1, K = 10,
              a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
              b1 = 0.2, b2 = 0.42, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
              o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

comp_out = expand.grid(m1 = seq(0.001, 0.1, by = 0.001), 
                       m2 = seq(0.001, 0.1, by = 0.001))


## Multi-strain model with H (i.e., S, P1 and P2) ----
comp_out_SPP <- as.data.frame(cbind(comp_out,
                                    matrix(0, 
                                           nrow = dim(comp_out)[1],
                                           ###dim(data)[1] is the number of row of data; [2] is col
                                           ncol = 3 + 1)))

names(comp_out_SPP) <- c("m1", "m2", "P1", "P2", "S", "Stability")

start_time <- Sys.time()
for(i in 1:dim(comp_out_SPP)[1]){
  temp_parms = parms
  temp_parms$m1 = comp_out_SPP[i, "m1"]
  temp_parms$m2 = comp_out_SPP[i, "m2"]
  
  #Calculate each equilibrium point
  E_C = f_E_C(temp_parms)
  E_S = f_E_S(temp_parms)
  E_P1 = f_E_P1(temp_parms)
  E_P2 = f_E_P2(temp_parms)
  
  Lambda_E_S = Eigen(Jacobian_mts(temp_parms, E_S))
  Lambda_E_P1 = Eigen(Jacobian_mts(temp_parms, E_P1))
  Lambda_E_P2 = Eigen(Jacobian_mts(temp_parms, E_P2))
  
  Stable_E = c()
  
 
  if(!is.na(Lambda_E_S) && is.finite(Lambda_E_S) && Lambda_E_S < 0 && all(E_S[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "S")
    comp_out_SPP[i, "Stability"] = "Stable"
    comp_out_SPP[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_S
  }
  if(!is.na(Lambda_E_P1) && is.finite(Lambda_E_P1) && Lambda_E_P1 < 0 && all(E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P1")
    comp_out_SPP[i, "Stability"] = "Stable"
    comp_out_SPP[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P1
  }
  if(!is.na(Lambda_E_P2) && is.finite(Lambda_E_P2) && Lambda_E_P2 < 0 && all(E_P2[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P2")
    comp_out_SPP[i, "Stability"] = "Stable"
    comp_out_SPP[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P2
  }
  
  
  
  comp_out_SPP[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  
  if(length(Stable_E) == 0){
    comp_out_SPP[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out_SPP[i, "Stability"] = "ASS"
  }
}

end_time <- Sys.time()
end_time - start_time

### Plot the result ----
ggplot(filter(comp_out_SPP), aes(x = m1, y = m2, fill = Stable_E)) +
  geom_raster() +
  labs(title = expression(), x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = expression("Intrinsic mortality rate of pathogen strain 2"~(m[2])))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_fill_manual("Equilibrium", values = final_colors, labels = outcome_labels) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  coord_fixed(ratio = 1)

ggsave("Result_2 Multi strain parameterspace with no H.png", width = 16, height = 11, units = "cm", dpi = 1600)

ggplot(comp_out_SPP, aes(x = m1, y = m2, fill = S)) +
  geom_tile() +
  labs(title = expression(), x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = expression("Intrinsic mortality rate of pathogen strain 2"~(m[2])))+
  scale_fill_gradient(high = "#013320",
                       low = "white") + #F9DDC8 #013320 ->S #006BA6 ->P2 #a50f15 ->P1
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  coord_fixed(ratio = 1)

ggsave("Result_2_2 Multi strain parameterspace with no H.png", width = 16, height = 11, units = "cm", dpi = 1600)

## Multi-strain model with H (i.e., S, P1, P2 and H) ----
comp_out_SPPH <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 6 + 1)))

names(comp_out_SPPH) <- c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S", "Stability")

start_time <- Sys.time()
for(i in 1:dim(comp_out_SPPH)[1]){
  temp_parms = parms
  temp_parms$m1 = comp_out_SPPH[i, "m1"]
  temp_parms$m2 = comp_out_SPPH[i, "m2"]
  # temp_parms$psi1 = comp_out_SPPH[i, "psi1"]
  # temp_parms$psi2 = comp_out_SPPH[i, "psi2"]
  
  #Calculate each equilibrium point
  E_C = f_E_C(temp_parms)
  E_S = f_E_S(temp_parms)
  E_P1 = f_E_P1(temp_parms)
  E_P2 = f_E_P2(temp_parms)
  E_P1H = f_E_P1H(temp_parms)
  E_P2H = f_E_P2H(temp_parms)
  
  Lambda_E_C = Eigen(Jacobian_full(temp_parms, E_C))
  Lambda_E_S = Eigen(Jacobian_full(temp_parms, E_S))
  Lambda_E_P1 = Eigen(Jacobian_full(temp_parms, E_P1))
  Lambda_E_P2 = Eigen(Jacobian_full(temp_parms, E_P2))
  Lambda_E_P1H = Eigen(Jacobian_full(temp_parms, E_P1H))
  Lambda_E_P2H = Eigen(Jacobian_full(temp_parms, E_P2H))
  
  Stable_E = c()
  
  if(!is.na(Lambda_E_C) && is.finite(Lambda_E_C) && Lambda_E_C < 0 && all(E_C[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "C")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_C
  }
  if(!is.na(Lambda_E_S) && is.finite(Lambda_E_S) && Lambda_E_S < 0 && all(E_S[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "S")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_S
  }
  if(!is.na(Lambda_E_P1) && is.finite(Lambda_E_P1) && Lambda_E_P1 < 0 && all(E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P1")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P1
  }
  if(!is.na(Lambda_E_P2) && is.finite(Lambda_E_P2) && Lambda_E_P2 < 0 && all(E_P2[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P2")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P2
  }
  if(!is.na(Lambda_E_P1H) && is.finite(Lambda_E_P1H) && Lambda_E_P1H < 0 && all(E_P1H[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P1H")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P1H
  }
  if(!is.na(Lambda_E_P2H) && is.finite(Lambda_E_P2H) && Lambda_E_P2H < 0 && all(E_P2H[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P2H")
    comp_out_SPPH[i, "Stability"] = "Stable"
    comp_out_SPPH[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = E_P2H
  }
  
  
  comp_out_SPPH[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  
  if(length(Stable_E) == 0){
    comp_out_SPPH[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out_SPPH[i, "Stability"] = "ASS"
  }
}

end_time <- Sys.time()
end_time - start_time

comp_out_SPPH[which(comp_out_SPPH$Stable_E == ""), "Stable_E"] = "U"

###如果前面都不想跑....
comp_out_SPPH = readRDS("Pre4A1_d002_psi08_rev")
comp_out_SPPH = readRDS("Pre4A1_d002_psi08_rev")
### Draw the boundary line ----
parms <- list(r = 1, K = 10,
              a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
              b1 = 0.2, b2 = 0.42, e1H = 0.5, e2H = 0.5,
              o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)
#When calculating IGR1 and 2, H must persist.
f_IGR1 = function(m1, m2){
  Leslie <- matrix(
    data = 
      with(parms, {
        P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
        A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
        B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
        C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
        E = B^2 - 4*A*C
        if(E < 0 || is.na(E)){
          H = 0
          S = m2/(e2*a2)
          P2 = (r/a2)*(1-(S/K))
        }else{
          H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
          H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
          H_12 = c(H1,H2)
          H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
          if(length(H_12) == 0){
            H = 0
            S = m2/(e2*a2)
            P2 = (r/a2)*(1-(S/K))
          } else {
            H = min(H_12) #Choose the smallest one which is more biologically meaningful
            S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
          }
        }
        #I set H = 0, because the S* from E1 = S*
        return(c(-(m1+o1), b1*H,
                 e1H*psi1*a1*S, e1*a1*S-b1*H-m1))
      }),
    nrow = 2,
    ncol = 2,
    byrow = T)
  EIGEN <- eigen(Leslie)
  max(EIGEN$values)
}
f_IGR2 = function(m1, m2){
  Leslie <- matrix(
    data = 
      with(parms, {
        P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
        A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
        B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
        C = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
        E = B^2 - 4*A*C
        if(E < 0 || is.na(E)){
          H = 0
          S = m1/(e1*a1)
          P1 = (r/a1)*(1-(S/K))
        }else{
          H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
          H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
          H_12 = c(H1,H2)
          H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
          if(length(H_12) == 0){
            H = 0
            S = m1/(e1*a1)
            P1 = (r/a1)*(1-(S/K))
          } else {
            H = min(H_12) #Choose the smallest one which is more biologically meaningful
            S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
          }
        }
        return(c(-(m2+o2), b2*H,
                 e2H*psi2*a2*S, e2*a2*S-b2*H-m2))
      }),
    nrow = 2,
    ncol = 2,
    byrow = T)
  EIGEN <- eigen(Leslie)
  max(EIGEN$values)
}
f_IGRH1 = function(m1){
  Leslie <- matrix(
    data = with(parms, {
      P1 = (r/a1)*(1-(m1/(e1*a1*K)))
      M = c(-(c1 * b1 * P1 + d), h1 * o1, b1 * P1, -(m1 + o1))
      return(M)
    }),
    nrow = 2,
    ncol = 2,
    byrow = T)
  EIGEN <- eigen(Leslie)
  max(EIGEN$values)
}
f_IGRH2 = function(m2){
  Leslie <- matrix(
    data = with(parms, {
      P2 = (r/a2)*(1-(m2/(e2*a2*K)))
      M = c(-(c2 * b2 * P2 + d), h2 * o2, b2 * P2, -(m2 + o2))
      return(M)
    }),
    nrow = 2,
    ncol = 2,
    byrow = T)
  EIGEN <- eigen(Leslie)
  max(EIGEN$values)
}

comp_out_SPPH$IGR1 = mapply(f_IGR1, comp_out_SPPH$m1, comp_out_SPPH$m2)
comp_out_SPPH$IGR2 = mapply(f_IGR2, comp_out_SPPH$m1, comp_out_SPPH$m2)
comp_out_SPPH$IGRH1 = mapply(f_IGRH1, comp_out_SPPH$m1)
comp_out_SPPH$IGRH2 = mapply(f_IGRH2, comp_out_SPPH$m2)


#To see the state of coexistence equilibrium (is stable/unstable, feasible/infeasible)
#Funciton setting
Jacobian = function(m1, m2, parms, E) {
  with(c(parms, E), {
    matrix(data =
             c(
               -c1*b1*P1-c2*b2*P2-d, h1*o1, h2*o2, -c1*b1*H, -c2*b2*H, 0,
               b1*P1, -(m1+o1), 0, b1*H, 0, 0,
               b2*P2, 0, -(m2+o2), 0, b2*H, 0,
               -b1*P1, e1H*psi1*a1*S, 0, e1*a1*S - b1*H - m1, 0, e1*a1*P1 + e1H*psi1*a1*P1H,
               -b2*P2, 0, e2H*psi2*a2*S, 0, e2*a2*S - b2*H - m2, e2*a2*P2 + e2H*psi2*a2*P2H,
               0, -psi1*a1*S, -psi2*a2*S, -a1*S, -a2*S, r*(1-S/K) - a1*P1 - psi1*a1*P1H - a2*P2 - psi2*a2*P2H - S*r/K
             ), nrow = 6, byrow = TRUE)
  })
}
Eigen = function(J) {
  if (any(is.na(J)) || any(is.infinite(J))){
    return(NA)
  }else{
    return(max(Re(eigen(J)$values)))
  }
}
f_C_State = function(m1, m2){
  E_C = with(parms, {
    D1 = (m1+o1)
    D2 = (m2+o2)
    A = D1*b1*e2H*psi2*a2*b2 - D2*b2*e1H*psi1*a1*b1
    B = D1*b1*e2*a2*D2 + m1*D1*e2H*psi2*a2*b2 - D2*b2*e1*a1*D1 - m2*D2*e1H*psi1*a1*b1
    C = m1*e2*a2*D1*D2 - m2*e1*a1*D2*D1
    E = B^2 - 4*A*C
    if(E < 0 || is.na(E)){
      H = NA
    }else{
      H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
      H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
      H_12 = c(H1,H2)
      H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
      if(length(H_12) == 0){
        H = NA
      } else {
        H = min(H_12) #Choose the smallest one which is more biologically meaningful
      }
    }
    A1 = (b1 * H) / (m1 + o1)
    A2 = (b2 * H) / (m2 + o2)
    B1 = h1 * o1 * A1 - c1 * b1 * H
    B2 = h2 * o2 * A2 - c2 * b2 * H
    D1 = (1 + psi1 * A1) * a1
    D2 = (1 + psi2 * A2) * a2
    S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
    #S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
    P1 = (d * H) / B1 - (B2 / B1) * (((B1 * r * (1 - S / K)) - D1 * d * H) / (D2 * B1 - D1 * B2))
    P2 = (B1 * r * (1 - S / K) - D1 * d * H) / (D2 * B1 - D1 * B2)
    P1H = A1 * P1
    P2H = A2 * P2
    
    return(setNames(
      c(H, P1H, P2H, P1, P2, S),
      c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
    )) 
  })
  if(any(is.na(E_C))){
    return("Infeasible")
  }else if(Eigen(Jacobian(m1, m2,parms, E_C)) > 0){
    return("Unstable")
  }else{
    return("Stable")
  }
}
comp_out_SPPH$C_State = mapply(f_C_State, comp_out_SPPH$m1, comp_out_SPPH$m2)
comp_out_SPPH$C_State_V = ifelse(comp_out_SPPH$C_State == "Stable", 1, 
                            ifelse(comp_out_SPPH$C_State == "Unstable", 0, -1))

#Stability of E_2H
f_E2H_State = function(m1, m2){
  E_2H = with(parms, {
    P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
    A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
    B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
    C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
    E = B^2 - 4*A*C
    if(E < 0 || is.na(E)){
      H = NA
    }else{
      H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
      H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
      H_12 = c(H1,H2)
      H_12 = H_12[H_12 > 0 & !is.na(H_12) & is.finite(H_12)]
      if(length(H_12) == 0){
        H = NA
      } else {
        H = min(H_12) #Choose the smallest one which is more biologically meaningful
      }
    }
    S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
    P2H = P2 * (b2 * H) / (m2 + o2)
    P1H = 0
    P1 = 0
    return(setNames(
      c(H, P1H, P2H, P1, P2, S),
      c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
    )) 
  })
  
  if(any(is.na(E_2H))){
    return("Infeasible")
  }else if(Eigen(Jacobian(m1, m2,parms, E_2H)) > 0){
    return("Unstable")
  }else{
    return("Stable")
  }
}
comp_out_SPPH$E2H_State = mapply(f_E2H_State, comp_out_SPPH$m1, comp_out_SPPH$m2)
comp_out_SPPH$E2H_State_V = ifelse(comp_out_SPPH$E2H_State == "Stable", 1, 
                              ifelse(comp_out_SPPH$E2H_State == "Unstable", 0, -1))

### Plot the result ----
ggplot(filter(comp_out_SPPH)) +
  geom_raster(mapping = aes(x = m1, y = m2, fill = Stable_E)) +
  geom_contour(filter(comp_out_SPPH, m2 < 0.065), mapping = aes(x = m1, y = m2, z = IGR2), breaks = 0, color = "grey", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 < 0.066), mapping = aes(x = m1, y = m2, z = IGR1), breaks = 0, color = "darkred", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 <0.043), mapping = aes(x = m1, y = m2, z = C_State_V), breaks = 0, color = "darkblue", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.065), mapping = aes(x = m1, y = m2, z = IGR2), breaks = 0, color = "purple", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.065), mapping = aes(x = m1, y = m2, z = IGR1), breaks = 0, color = "purple", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.064), mapping = aes(x = m1, y = m2, z = IGRH1), breaks = 0, color = "darkgreen", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m1 > 0.048), mapping = aes(x = m1, y = m2, z = IGRH2), breaks = 0, color = "darkgreen", linewidth = 1.2)+
  labs(title = expression(), x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = expression("Intrinsic mortality rate of pathogen strain 2"~(m[2])), fill = "Equilibrium")+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_fill_manual(values = final_colors, labels = outcome_labels) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  coord_fixed(ratio = 1)

ggsave("Result_3 Multi strain parameterspace with H.png", width = 16, height = 11, units = "cm", dpi = 1600)

ggplot(comp_out_SPPH, aes(x = m1, y = m2, fill = S)) +
  geom_tile() +
  geom_contour(filter(comp_out_SPPH, m2 < 0.065), mapping = aes(x = m1, y = m2, z = IGR2), breaks = 0, color = "grey", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 < 0.066), mapping = aes(x = m1, y = m2, z = IGR1), breaks = 0, color = "darkred", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 <0.043), mapping = aes(x = m1, y = m2, z = C_State_V), breaks = 0, color = "darkblue", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.065), mapping = aes(x = m1, y = m2, z = IGR2), breaks = 0, color = "purple", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.065), mapping = aes(x = m1, y = m2, z = IGR1), breaks = 0, color = "purple", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m2 > 0.064), mapping = aes(x = m1, y = m2, z = IGRH1), breaks = 0, color = "darkgreen", linewidth = 1.2)+
  geom_contour(filter(comp_out_SPPH, m1 > 0.048), mapping = aes(x = m1, y = m2, z = IGRH2), breaks = 0, color = "darkgreen", linewidth = 1.2)+
  labs(title = expression(), x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = expression("Intrinsic mortality rate of pathogen strain 2"~(m[2])))+
  scale_fill_gradient(high = "#013320",
                      low = "white") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.1, by = 0.02)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  coord_fixed(ratio = 1)

ggsave("Result_3_3 Multi strain parameterspace with H.png", width = 16, height = 11, units = "cm", dpi = 1600)

#Result 4: the relationship between delta_H and m ----
M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)# * epsilon
    return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
  })
}

## H_press to flip equilibrium along the parameter space ----
### Read the data to further subset different initial states and parameters
P.space = 
  readRDS("Pre4A1_d002_psi08_rev") %>%
  #readRDS("Pre4A1_forPSpace") %>%
  filter(Stable_E == "P2H,C")
P.space = P.space[,-c(9:10)]
P.space$H_press = 0
head(P.space)

### Model parameters
times <- c(0, 16000)
state <- c(H = 0, P1H = 0, P2H = 0, P1 = 0, P2 = 0, S = 0)

parms <- c(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
# parameter = expand_grid(P.space[, c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S")])

#function that will output the competition outcome
f_outcome = function(state, parms){
  pop_size = ode(func = M2,
                 times = times,
                 y = state,
                 parms = parms,
                 maxsteps = 80000)
  Outcome =
    paste(ifelse(pop_size[nrow(pop_size), c("H", "P1H", "P2H", "P1", "P2")] > 1e-7, "T", "F"), collapse = "")
  return(Outcome)
}

#function of bisection method
Bisection_H = function(min_H, max_H, state, parms){
  ###The goal is to find a maximum value of H that make competition outcome == E2H
  ###Min
  Min_state = state
  Min_state["H"] = min_H
  min_E = f_outcome(state = Min_state, parms = parms)
  
  ###Max
  Max_state = state
  Max_state["H"] = max_H
  max_E = f_outcome(state = Max_state, parms = parms)
  
  if(min_E == max_E){
    print(paste0("The equilibrium of minimum H is equal to maximum H (", min_H, ")"))
    return(NA)
  }
  
  loop_count = 0
  while(max_H - min_H > 1e-9 && loop_count < 100){
    ###Mid
    mid_H = (min_H + max_H)/2
    Mid_state = state
    Mid_state["H"] = mid_H
    mid_E = f_outcome(state = Mid_state, parms = parms)
    
    ###determine the value of mid_H
    if(mid_E == min_E){
      min_H = mid_H
    }else if(mid_E == max_E){
      max_H = mid_H
    }else{
      return(NA)
    }
    loop_count = loop_count + 1
  }
  return((min_H + max_H)/2)
}

Start_time = Sys.time()
for (i in 1:dim(P.space)[1]) {
  #setting the initial condition
  temp_state = unlist(P.space[i,  c("H", "P1H", "P2H", "P1", "P2", "S")])
  
  #setting the parameters
  temp_parms = parms
  temp_parms[c("m1", "m2")] = P.space[i, c("m1", "m2")]
  
  #run with the bisection method
  P.space[i,]$H_press = Bisection_H(min_H = min(P.space[i, "H"]),
                                    max_H = 8,
                                    state = temp_state,
                                    parms = temp_parms)
  if(i %% 10 == 0) print(paste0("Total: ", nrow(P.space), ", Now: ", i))
}
Ending_time = Sys.time()
Ending_time - Start_time
saveRDS(P.space, file = "")


readRDS("Result 4 (H_press to flip equilibrium)") %>%
  filter(round(m2, 3) == 0.02) %>%
  mutate(delta_H = H_press - H) %>%
  ggplot()+
  geom_line(mapping = aes(x = m1, y = delta_H), linewidth = 0.8, color = "darkgreen")+
  #scale_y_continuous(limits = c(0, 7)) +
  #scale_x_continuous(breaks = c(seq(0, 0.015, by = 0.003))) +
  labs(title = expression(m[2] == 0.02), x = expression(m[1]), y = expression(Delta~H))

ggsave("Result_4 H_press to flip equilibrium_2.png", width = 15, height = 10, units = "cm", dpi = 800)


## H_press to flip equilibrium along the time series ----
times = seq(0, 10000, by = 0.1)

###state----
ini_High_H <- c(H = 0.5, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_1
ini_Low_H <- c(H = 0.01, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_2

###parms----
parms_E1H_E2H =
  c(r = 1, K = 10,
    a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.42, m1 = 0.04, m2 = 0.01, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)
parms_EC_E2H =
  c(r = 1, K = 10,
    a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.01, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)


###Function of time Series----
TimeSeries = function(times, state, parms){
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
      dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
  }
  pop_size = ode(func = M2, times = times, y = state, parms = parms)
  pop_size %>%
    as.data.frame() #%>%
  #filter(time %% 1 == 0)
}


### Run the data set----

High_H_EC_E2H = TimeSeries(times, ini_High_H, parms_EC_E2H)

#function that will output the competition outcome
f_outcome = function(state, parms){
  pop_size = ode(func = M2,
                 times = times,
                 y = state,
                 parms = parms,
                 maxsteps = 20000)
  Outcome =
    paste(ifelse(pop_size[nrow(pop_size), c("H", "P1H", "P2H", "P1", "P2")] > 1e-7, "T", "F"), collapse = "")
  return(Outcome)
}

#function of bisection method
Bisection_H = function(min_H, max_H, state, parms){
  ###The goal is to find a maximum value of H that make competition outcome == E2H
  ###Min
  Min_state = state
  Min_state["H"] = min_H
  min_E = f_outcome(state = Min_state, parms = parms)
  
  ###Max
  Max_state = state
  Max_state["H"] = max_H
  max_E = f_outcome(state = Max_state, parms = parms)
  
  if(min_E == max_E){
    print(paste0("The equilibrium of minimum H is equal to maximum H (", min_H, ")"))
    return(NA)
  }
  
  loop_count = 0
  while(max_H - min_H > 1e-6 && loop_count < 80){
    ###Mid
    mid_H = (min_H + max_H)/2
    Mid_state = state
    Mid_state["H"] = mid_H
    mid_E = f_outcome(state = Mid_state, parms = parms)
    
    ###determine the value of mid_H
    if(mid_E == min_E){
      min_H = mid_H
    }else if(mid_E == max_E){
      max_H = mid_H
    }else{
      return(NA)
    }
    loop_count = loop_count + 1
  }
  return((min_H + max_H)/2)
}

#The timing of invasion
x = seq(1, 4, by = 0.01)
Invasion_timing = 10^x

#Set up
High_H_EC_E2H$H_press = 0

#Run a for-loop
Start_time = Sys.time()
for (i in Invasion_timing) {
  #setting the initial condition
  #temp_state = round(unlist(High_H_EC_E2H[i,  c("H", "P1H", "P2H", "P1", "P2", "S")]), 8) #use unlist() to make values become vectors
  temp_state = round(unlist(High_H_EC_E2H[round(High_H_EC_E2H$time, 2) == round(i, 2),  c("H", "P1H", "P2H", "P1", "P2", "S")]), 8)
  
  #run with the bisection method
  High_H_EC_E2H[i,]$H_press = 
    Bisection_H(min_H = High_H_EC_E2H[i, "H"],
                max_H = 25,
                state = temp_state,
                parms = parms_EC_E2H)
  print(paste0("Total: ", length(Invasion_timing), ", Now: ", which(Invasion_timing == i)))
}
End_time = Sys.time()
End_time - Start_time

#saveRDS(High_H_EC_E2H, "")
##Plotting----
#High_H_EC_E2H = readRDS("Result 4 High_H_EC_E2H_invade")
High_H_EC_E2H %>%
  filter(H_press != 0) %>%
  ggplot()+
  geom_line(mapping = aes(x = time, y = H_press-H))+
  geom_point(mapping = aes(x = time, y = H_press-H))+
  labs(x = "Adding time", y = "Addtional H")+
  scale_x_log10()
ggsave("High_H_EC_E2H_invasion.png", width = 14, height = 10, units = "cm", dpi = 800)
