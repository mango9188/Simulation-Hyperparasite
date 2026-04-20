#### This script is for master thesis.
# Read the package and plot setting
library(tidyverse)
source("M_Theme setting.R", encoding = 'CP950', echo = T)
library(patchwork)
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
0
# Single strain model with no H (i.e., S and P1)----
parms = list(
  r = 1, K = 10, a1 = 0.35, e1 = 0.5, m1 = 0.05)

#Create data frame to expand m1
comp_out = expand.grid(m1 = seq(0.01, 0.5, by = 0.01))
comp_out = as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 2+1)))
names(comp_out) = c("m1", "P1", "S", "Stability")

comp_out


for (i in 1:dim(comp_out)[1]) {
  temp_parms = parms
  temp_parms["m1"]  = comp_out[i, "m1"]
  E_P1 = f_E_P1(temp_parms)
  E_S = f_E_S(temp_parms)
  Lambda_E_P1 = Eigen(Jacobain_sgsN(temp_parms, E_P1))
  Lambda_E_S = Eigen(Jacobain_sgsN(temp_parms, E_S))
  
  Stable_E = c()
  
  if(!is.na(Lambda_E_P1) && is.finite(Lambda_E_P1) && Lambda_E_P1 < 0 && all(E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "P1")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("P1", "S")] = f_E_P1(temp_parms)[c("P1", "S")]
  }
  if(!is.na(Lambda_E_S) && is.finite(Lambda_E_S) && Lambda_E_S < 0 && all(E_S[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    Stable_E = c(Stable_E, "S")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("P1", "S")] = f_E_P1(temp_parms)[c("P1", "S")]
  }

  comp_out[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  
  if(length(Stable_E) == 0){
    comp_out[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out[i, "Stability"] = "ASS"
  }
}

comp_out$System = "SP"

#Bifurcation plot for S and P1
D = 
  comp_out %>%
  select(c(m1, P1, S)) %>% #P1, P2, P1H, P2H, H, S
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1))
#gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()

ggplot(D, aes(x = m1, y = Abundance, color = Species)) +
  geom_line(mapping = aes(x = m1, y = Abundance, color = Species), lwd = 1) +
  #geom_line(filter(D, Species == "total"), mapping = aes(x = m1, y = Abundance, color = Species), lwd = 0.8, linetype = 2) +
  #scale_linetype_manual(values = c("Stable" = "solid", "Unstable" = "dashed")) +
  labs(title = "No natural enemies", x = "Pathogen mortality", y = "Abundance", color = "Species")+
  scale_y_continuous() +
  scale_x_continuous() + #breaks = c(seq(0.2, 1, by = 0.2))
  scale_colour_manual(labels = 
                        c("P1" = "Pathogen", "S" = "Host"),
                      values = c("P1" = "#BCAAA4", "S" = "#00AF66"))+
  theme(axis.title.y.right = element_text(angle = 90))

ggsave("5min No natural enemies.png", width = 15, height = 11, units = "cm", dpi = 800, bg = "transparent")

# Single strain model with H (i.e., S, P1, P1H, and H)----
# parms = list(
#   r = 1, K = 10,
#   a1 = 0.38, psi1 = 0.8, e1 = 0.5, 
#   b1 = 0.2, m1 = 0.05, e1H = 0.5, 
#   o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.028, DL = 0) #c1 = 0.9

parms = list(
  r = 1, K = 10,
  a1 = 0.25, psi1 = 1, e1 = 0.5, 
  b1 = 0.2, e1H = 0.5, 
  o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.01)

## Create data frame to expand m1----
comp_out = expand.grid(m1 = seq(0.001, 0.1, by = 0.001))
comp_out = as.data.frame(cbind(comp_out,
                               matrix(0, 
                                      nrow = dim(comp_out)[1],
                                      ###dim(data)[1] is the number of row of data; [2] is col
                                      ncol = 4)))
names(comp_out) = c("m1", "H", "P1H", "P1", "S")

comp_out

for (i in 1:dim(comp_out)[1]) {
  temp_parms = parms
  temp_parms["m1"]  = comp_out[i, "m1"]
  
  E_P1H = f_E_P1H(temp_parms)
  Lambda_P1H = Eigen(Jacobain_sgs(temp_parms, E_P1H))
  E_P1 = f_E_P1(temp_parms)
  Lambda_P1 = Eigen(Jacobain_sgs(temp_parms, E_P1))
  E_S = f_E_S(temp_parms)
  Lambda_S = Eigen(Jacobain_sgs(temp_parms, E_S))
  
  Stable_E = c()
  if (!is.na(Lambda_P1H) == T && is.finite(Lambda_P1H) == T && Lambda_P1H < 0 && all(!is.na(E_P1H) == T)){
    Stable_E = c(Stable_E, "P1H")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_P1H[c("H", "P1H", "P1", "S")]
  }
  if(!is.na(Lambda_P1) == T && is.finite(Lambda_P1) == T && Lambda_P1 < 0 && all(!is.na(E_P1) == T)){
    Stable_E = c(Stable_E, "P1")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_P1[c("H", "P1H", "P1", "S")]
  }
  if(!is.na(Lambda_S) == T && is.finite(Lambda_S) == T && Lambda_S < 0 && all(!is.na(E_S) == T)){
    Stable_E = c(Stable_E, "S")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_S[c("H", "P1H", "P1", "S")]
  }
  
  comp_out[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  if(length(Stable_E) == 0){
    comp_out[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out[i, "Stability"] = "ASS"
  }
}

###It is possible that data contains NA because m1 can be large enough that P1 would no longer persist H.
comp_out = mutate(comp_out, P.total = P1+P1H)

## Bifurcation plot----
D = 
  comp_out %>%
  select(c(m1, P1, P1H, H, S, P.total)) %>% #P1, P2, P1H, P2H, H, S
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1))
#gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()
D$Group = 0
D[D[, "Species"] == "H", ]$Group = 1
D[D[, "Species"] == "S", ]$Group = 1
D[D[, "Species"] == "P1", ]$Group = 2
D[D[, "Species"] == "P1H", ]$Group = 2
D[D[, "Species"] == "P.total", ]$Group = 2



ggplot(D, aes(x = m1, y = Abundance, color = Species)) +
  geom_line(mapping = aes(x = m1, y = Abundance, color = Species), lwd = 1) +
  labs(x = "Intrinsic mortality rate of pathogen (m)", y = "Abundance", color = "Species")+
  geom_vline(xintercept = 0.0755, color = "black", linetype = 2, linewidth = 1) +
  #geom_vline(xintercept = 0.6045, color = "darkgreen", linetype = 2, linewidth = 1)+
  #breaks = c(seq(0.2, 1, by = 0.2))
  #ylim(0,3.1)+
  #labs(title = expression(r + "5%"))+
  #xlim(0, 0.08)+
  #ylim(0, 5)+
  facet_grid(Group ~ ., scales = "free_y")+
  scale_colour_manual(values = State_values, labels = State_labels_SS) +
  theme(axis.title.y.right = element_text(angle = 90), legend.position = "bottom", strip.background = element_blank(), strip.text = element_blank())


ggsave("Single strain bifurcation (Dash line).png", width = 15, height = 11, units = "cm", dpi = 800)
# ggsave("", width = 20, height = 11, units = "cm", dpi = 800)
# ggsave("", width = 14, height = 10, units = "cm", dpi = 800)


## To find the boundary of H invasion----
parms = list(
  r = 1, K = 10,
  a1 = 0.51, psi1 = 0.8, e1 = 0.5, 
  b1 = 0.42, e1H = 0.5, 
  o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.02, DL = 0)
#Remember to use the parameter without m1!
f_Max.m = function(m1){
  with(parms, {
    P_a = (r / a1) * (1 - (m1 / (e1 * a1 * K)))
    P_b = (m1 + o1) * d /(h1 * o1 * b1 - c1 * b1 * (m1 + o1))
    return(P_a - P_b)
    })
}


(Max.m = uniroot(f_Max.m, interval = c(0.05, 0.07), tol = 1e-16)$root)



## Find dP*/dm to see the effectness of hydra effect----
### For single strain model----
parms = list(
  r = 1, K = 10,
  a1 = 0.5, psi1 = 1, e1 = 0.5, 
  b1 = 0.45, e1H = 0.5, m1 = 0.05643642,
  o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.02, DL = 0)

f_E_P1H = function(parms){
  with(parms,{
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

###Method 1
E1H_P1 = expression((d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1)))
E1H_dP1_dm = D(E1H_P1, "m1") #Output -> d/(b1 * h1 * o1 - c1 * b1 * (o1 + m1)) + (d * (o1 + m1)) * (c1 * b1)/(b1 * h1 * o1 - c1 * b1 * (o1 + m1))^2
with(parms,{
  eval(E1H_dP1_dm)
})

E1_P1 = expression((r/a1)*(1-(m1/(e1*a1)/K)))
E1_dP1_dm = D(E1_P1, "m1") #Output -> -((r/a1) * (1/(e1 * a1)/K)) is always <0
with(parms,{
  eval(E1_dP1_dm)
})#This will always be -0.8, which means that there is no hydra effect at E1.

E1H_P1H = expression(P1 * (b1 * H) / (m1 + o1))
E1H_dP1H_dm = D(E1H_P1H, "m1")
with(parms,{
  eval(E1H_dP1H_dm)
})

###Method 2
library(Deriv)
E1H_P1 = function(m1) {(d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))}
E1H_dP_dm = Deriv(E1H_P1, "m1")


###Add into the for-loop with method 1
#function setting
E1H_P1 = expression((d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1)))
E1H_dP1_dm = D(E1H_P1, "m1")
E1_P1 = expression((r/a1)*(1-(m1/(e1*a1)/K)))
E1_dP1_dm = D(E1_P1, "m1")


comp_out = expand.grid(m1 = seq(0.001, 0.1, by = 0.001))
comp_out = as.data.frame(cbind(comp_out,
                               matrix(0, 
                                      nrow = dim(comp_out)[1],
                                      ###dim(data)[1] is the number of row of data; [2] is col
                                      ncol = 5)))
names(comp_out) = c("m1", "H", "P1H", "P1", "S", "dP1_dm")
comp_out
for (i in 1:dim(comp_out)[1]) {
  temp_parms = parms
  temp_parms["m1"]  = comp_out[i, "m1"]
  
  E_P1H = f_E_P1H(temp_parms)
  Lambda_P1H = Eigen(Jacobain_sgs(temp_parms, E_P1H))
  E_P1 = f_E_P1(temp_parms)
  Lambda_P1 = Eigen(Jacobain_sgs(temp_parms, E_P1))
  E_S = f_E_S(temp_parms)
  Lambda_S = Eigen(Jacobain_sgs(temp_parms, E_S))
  
  Stable_E = c()
  if (!is.na(Lambda_P1H) == T && is.finite(Lambda_P1H) == T && Lambda_P1H < 0 && all(!is.na(E_P1H) == T)){
    Stable_E = c(Stable_E, "P1H")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_P1H[c("H", "P1H", "P1", "S")]
    comp_out[i, "dP1_dm"] = with(temp_parms,{eval(E1H_dP1_dm)})
  }
  if(!is.na(Lambda_P1) == T && is.finite(Lambda_P1) == T && Lambda_P1 < 0 && all(!is.na(E_P1) == T)){
    Stable_E = c(Stable_E, "P1")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_P1[c("H", "P1H", "P1", "S")]
    comp_out[i, "dP1_dm"] = with(temp_parms,{eval(E1_dP1_dm)})
  }
  if(!is.na(Lambda_S) == T && is.finite(Lambda_S) == T && Lambda_S < 0 && all(!is.na(E_S) == T)){
    Stable_E = c(Stable_E, "S")
    comp_out[i, "Stability"] = "Stable"
    comp_out[i, c("H", "P1H", "P1", "S")] = E_S[c("H", "P1H", "P1", "S")]
    comp_out[i, "dP1_dm"] = 0
  }
  
  comp_out[i, "Stable_E"] = paste(Stable_E, collapse = ",")
  if(length(Stable_E) == 0){
    comp_out[i, "Stability"] = "Unstable"
  }else if(length(Stable_E) != 1){
    comp_out[i, "Stability"] = "ASS"
  }
}

#plot the result
D = 
  comp_out %>%
  select(c(m1, dP1_dm)) %>% #P1, P1H, H, S, P.total
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1))
#gather(key = Species, value = Abundance, -c(a1, r)) %>% #using gather()

ggplot(D, aes(x = m1, y = Abundance, color = Species)) +
  geom_line(mapping = aes(x = m1, y = Abundance, color = Species), lwd = 1) +
  labs(title = expression(d == 0.02), x = "Mortality (m)", y = "Abundance", color = "Species")+
  #geom_vline(xintercept = 0.1462919, color = "black", linetype = 2, linewidth = 1) +
  #breaks = c(seq(0.2, 1, by = 0.2))
  #ylim(0,3)+
  scale_colour_manual(labels = 
                        c("P1" = "Pathogen", "P1H" = "Pathogen/H",
                          "S" = "Host", "H" = "Hyperparasite",
                          "dP1_dm" = expression(frac(dP[1], dm))),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "S" = "#00AF66", "H" = "#C03728",
                                 "P.total" = "black",
                                 "dP1_dm" = "black")) +
  #scale_shape_manual(values = c("Stable" = 16, "Unstable" = 3)) + 
  theme(axis.title.y.right = element_text(angle = 90))

### For multiple strain model----
library(numDeriv)
parms <- list(r = 1, K = 10,
              a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
              b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
              o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

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

func_for_m1 = function(m1_value) {
  local_parms = parms
  local_parms$m1 = m1_value
  return(f_E_C(local_parms))
}
deriv_results = jacobian(func = func_for_m1, x = 0.01)
result_df = data.frame(
  Variable = c('H', 'P1H', 'P2H', 'P1', 'P2', 'S'),
  Value = as.numeric(f_E_C(parms)), 
  Derivative_wrt_m1 = as.vector(deriv_results)
)
print(result_df)

func_for_m2 = function(m2_value) {
  local_parms = parms
  local_parms$m2 = m2_value
  return(f_E_C(local_parms))
}
deriv_results = jacobian(func = func_for_m2, x = 0.01)

result_df = data.frame(
  Variable = c('H', 'P1H', 'P2H', 'P1', 'P2', 'S'),
  Value = as.numeric(f_E_C(parms)), 
  Derivative_wrt_m2 = as.vector(deriv_results)
)
print(result_df)


Data = expand.grid(m2 = seq(0, 1, by = 0.01))
Data = data.frame(Data,
                  Value = 0,
                  Derivative_wrt_m2 = 0)
for (i in dim(Data)[1]) {
  deriv_results = jacobian(func = func_for_m2, x = 0.01)
  
    Variable = c('H', 'P1H', 'P2H', 'P1', 'P2', 'S'),
    Value = as.numeric(f_E_C(parms)), 
    Derivative_wrt_m2 = as.vector(deriv_results)
  )
}

# how parameters affect pathogen resurgence? (sensitivity analysis)-----
parms = list(
  r = 1, K = 10,
  a1 = 0.25, psi1 = 0.8, e1 = 0.5, 
  b1 = 0.2, e1H = 0.5, 
  o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.03)
0
## Method 1, simple for-loop -----
#To understand how hyperparasite induced mortality affect pathogen resurgence
#Remember to use the parameter without m1!
m1.min = 0.001

comp_out = data.frame(o1 = c(parms$o1 - parms$o1*0.1, parms$o1, parms$o1 + parms$o1*0.1),
                      m1.min = m1.min,
                      m1.max = 0,
                      P.total.min = 0,
                      P.total.max = 0)

for (i in 1:dim(comp_out)[1]) {
  temp_parms = parms
  temp_parms$o1 = comp_out[i, "o1"]
  f_Max.m = function(m1){
    with(temp_parms, {
      P_a = (r / a1) * (1 - (m1 / (e1 * a1 * K)))
      P_b = (m1 + o1) * d /(h1 * o1 * b1 - c1 * b1 * (m1 + o1))
      return(P_a - P_b)
    })
  }
  comp_out[i, "m1.max"] = uniroot(f_Max.m, interval = c(0.01, 0.08), tol = 1e-8)$root
  #Get the critical value of m that makes hyperparasite persist.
  temp_parms$m1 = comp_out[i, "m1.min"]
  comp_out[i, "P.total.min"] = f_E_P1H(temp_parms)["P1H"] + f_E_P1H(temp_parms)["P1"]
  temp_parms$m1 = comp_out[i, "m1.max"]
  comp_out[i, "P.total.max"] = f_E_P1(temp_parms)["P1"]
}

comp_out = mutate(comp_out, slope = (P.total.max-P.total.min) / (m1.max-m1.min))

comp_out
ggplot(comp_out, mapping = aes(x = factor(o1), y = slope))+
  geom_col(position = "identity")+
  coord_flip()

## Method 2, using function----
Get_slope = function(parmsA){
  # Define minimum m1
  m1.min = 0.001
  
  # Define maximum m1
  f_Max.m = function(parmsA){
    with(parmsA, {
      A = (c1*b1*r)/(e1*a1*K)
      B = -c1*b1*r - a1*d - (h1*o1*b1*r)/(e1*a1*K) + (c1*b1*r*o1)/(e1*a1*K)
      C = h1*o1*b1*r - a1*d*o1 - c1*b1*r*o1
      E = B^2 - 4*A*C
      if(E < 0 || is.na(E)){
        m1.max = NA
      }else{
        R1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
        R2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
        R = c(R1,R2)
        R = R[R > 0 & !is.na(R) & is.finite(R)]
        if(length(R) == 0){
          m1.max = NA
        } else {
          m1.max = min(R) #Choose the smallest one which is biologically meaningful
        }
      }
      return(m1.max)
    })
  }
  m1.max = f_Max.m(parmsA)
  
  # Define parameters
  temp_parms = parmsA
  
  # Calculate state variables
  temp_parms$m1 = m1.min
  H.min = f_E_P1H(temp_parms)["H"]
  P.min = f_E_P1H(temp_parms)["P1"]
  PH.min = f_E_P1H(temp_parms)["P1H"]
  S.min = f_E_P1H(temp_parms)["S"]
  P.total.min = f_E_P1H(temp_parms)["P1H"] + f_E_P1H(temp_parms)["P1"]
  
  temp_parms$m1 = m1.max
  P.total.max = f_E_P1(temp_parms)["P1"]
  S.max = f_E_P1(temp_parms)["S"]

  # Calculate slope
  P.slope = (P.total.max - P.total.min) / (m1.max - m1.min)
  S.slope = (S.max - S.min) / (m1.max - m1.min)

  
  #Return value
  return(setNames(
    c(m1.min, m1.max, H.min, P.min, PH.min, S.min, S.max, P.total.min, P.total.max, P.slope, S.slope),
    c("m1.min", "m1.max", "H.min", "P.min", "PH.min", "S.min", "S.max", "P.total.min", "P.total.max", "P.slope", "S.slope")))
}

result_list = list()
start_time = Sys.time()
for (i in names(unlist(parms))) {
  # Set and reset the parameters
  parms = parms
  parms.in = parms
  parms.de = parms
  
  parms.de[i] = unlist(parms[i])*0.95 #specific parm - 5%
  parms.in[i] = unlist(parms[i])*1.05 #specific parm + 5%
  
  Get_slope(parms.de)
  Get_slope(parms.in)
  
  result_list[[length(result_list) + 1]] = 
    data.frame(parameter = i,
               value = "-5%",
               m1.min = Get_slope(parms.de)[1],
               m1.max = Get_slope(parms.de)[2],
               H.min = Get_slope(parms.de)[3],
               P.min = Get_slope(parms.de)[4],
               PH.min = Get_slope(parms.de)[5],
               S.min = Get_slope(parms.de)[6],
               S.max = Get_slope(parms.de)[7],
               P.total.min = Get_slope(parms.de)[8],
               P.total.max = Get_slope(parms.de)[9],
               S.slope = Get_slope(parms.de)[11],
               P.slope = Get_slope(parms.de)[10],
               delta_slope = (Get_slope(parms.de)[10] - Get_slope(parms)[10])/ Get_slope(parms)[10])
  
  # result_list[[length(result_list) + 1]] = 
  #   data.frame(parameter = i,
  #              value = "Origin",
  #              m1.min = Get_slope(parms)[1],
  #              m1.max = Get_slope(parms)[2],
  #              P.total.min = Get_slope(parms)[3],
  #              P.total.max = Get_slope(parms)[4],
  #              slope = Get_slope(parms)[5])
  
  result_list[[length(result_list) + 1]] = 
    data.frame(parameter = i,
               value = "+5%",
               m1.min = Get_slope(parms.in)[1],
               m1.max = Get_slope(parms.in)[2],
               H.min = Get_slope(parms.in)[3],
               P.min = Get_slope(parms.in)[4],
               PH.min = Get_slope(parms.in)[5],
               S.min = Get_slope(parms.in)[6],
               S.max = Get_slope(parms.in)[7],
               P.total.min = Get_slope(parms.in)[8],
               P.total.max = Get_slope(parms.in)[9],
               S.slope = Get_slope(parms.in)[11],
               P.slope = Get_slope(parms.in)[10],
               delta_slope = (Get_slope(parms.in)[10] - Get_slope(parms)[10])/ Get_slope(parms)[10])

}
end_time = Sys.time()
end_time - start_time



result_df = do.call(rbind, result_list)

####Plotting relative slope (compared to the original parameter set)----
result_df %>%
  mutate(parameter = paste(parameter, value)) %>%
  ggplot()+
  geom_col(mapping = aes(x = factor(parameter), y = delta_slope))+
  coord_flip()+
  labs(x = "Parameters", y = expression(Delta*"Slope"))

####Plotting slope value----

param_labels = c("r" = "r",
                  "psi1" = expression(psi),
                  "o1" = expression(omega),
                  "K" = "K",
                  "h1" = "h",
                  "e1H" = expression(e[H]),
                  "e1" = expression(e),
                  "d" = "d",
                  "c1" = "c",
                  "b1" = expression(beta),
                  "a1" = expression(alpha))

Sensitivity_P = 
result_df %>%
  # filter(value == "+5%") %>%
  #mutate(parameter = paste(parameter, value)) %>%
  ggplot()+
  geom_col(mapping = aes(x = factor(parameter), y = P.slope, fill = value), position = "dodge")+
  geom_hline(yintercept = Get_slope(parms)["P.slope"], color = "black", linetype = 2, linewidth = 1) +
  labs(x = "Parameters", y = "Slope of pathogen biomass", fill = "Variation")+
  scale_x_discrete(labels = param_labels)+
  coord_flip()+
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank()) 

Sensitivity_S = 
result_df %>%
  # filter(value == "+5%") %>%
  #mutate(parameter = paste(parameter, value)) %>%
  ggplot()+
  geom_col(mapping = aes(x = factor(parameter), y = S.slope, fill = value), position = "dodge")+
  geom_hline(yintercept = Get_slope(parms)["S.slope"], color = "black", linetype = 2, linewidth = 1) +
  labs(x = "Parameters", y = "Slope of plant biomass", fill = "Variation")+
  scale_x_discrete(labels = param_labels)+
  coord_flip()
  
####Plotting initial abundance of H-----
result_df %>%
  filter(value == "+5%") %>%
  #mutate(parameter = paste(parameter, value)) %>%
  ggplot()+
  geom_col(mapping = aes(x = factor(parameter), y = H.min))+
  geom_hline(yintercept = Get_slope(parms)["H.min"], color = "black", linetype = 2, linewidth = 1) +
  labs(x = "Parameters + 5%", y = "Initial abundance ratio of H")+
  scale_x_discrete(labels = c("r" = "r",
                              "psi1" = expression(psi),
                              "o1" = expression(omega),
                              "K" = "K",
                              "h1" = "h",
                              "e1H" = expression(e[H]),
                              "e1" = expression(e),
                              "d" = "d",
                              "c1" = "c",
                              "b1" = expression(beta),
                              "a1" = expression(alpha)))


####Patch work-----

Sensitivity_S + Sensitivity_P + 
  plot_layout(widths = c(0.2)) & 
  theme(plot.tag = element_text(size = 10, face = "bold"))

Sensitivity_S + Sensitivity_P + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("Result_1_2 Sensitivity analysis.png", width = 20, height = 12, units = "cm", dpi = 800) #15, 16
