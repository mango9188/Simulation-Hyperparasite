library(tidyverse)
library(deSolve)
library(patchwork)
source("M_Theme setting.R", encoding = 'CP950', echo = T)

###m_i == (m_i + a * ep_i)
###where m_i is natural mortality, a is concentration of pesticide, and ep_i is the sensitivity of P_i to the pesticide.

# Function setting----
{
  Jacobian_full = function(parms, E) {
    with(c(parms, E), {
      matrix(data =
               c(
                 -c1*b1*P1-c2*b2*P2-d, h1*o1, h2*o2, -c1*b1*H, -c2*b2*H, 0,
                 b1*P1, -((m1 + a * ep1)+o1), 0, b1*H, 0, 0,
                 b2*P2, 0, -((m2 + a * ep2)+o2), 0, b2*H, 0,
                 -b1*P1, e1H*psi1*a1*S, 0, e1*a1*S - b1*H - (m1 + a * ep1), 0, e1*a1*P1 + e1H*psi1*a1*P1H,
                 -b2*P2, 0, e2H*psi2*a2*S, 0, e2*a2*S - b2*H - (m2 + a * ep2), e2*a2*P2 + e2H*psi2*a2*P2H,
                 0, -psi1*a1*S, -psi2*a2*S, -a1*S, -a2*S, r*(1-S/K) - a1*P1 - psi1*a1*P1H - a2*P2 - psi2*a2*P2H - S*r/K
               ), nrow = 6, byrow = TRUE)
    })
  } ##for full model (S, P1, P2, P1H, P2H, and H)
  
  Jacobian_mts = function(parms, E) {
    with(c(parms, E), {
      matrix(data =
               c(
                 e1*a1*S - (m1 + a * ep1), 0, e1*a1*P1,
                 0, e2*a2*S - (m2 + a * ep2), e2*a2*P2,
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
      J22 = -((m1 + a * ep1)+o1)
      J23 = b1*H
      J24 = 0
      J31 = -b1*P1 
      J32 = e1H*psi1*a1*S
      J33 = e1*a1*S-b1*H-(m1 + a * ep1)
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
      J11 = e1*a1*S-(m1 + a * ep1)
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
      D1 = ((m1 + a * ep1)+o1)
      D2 = ((m2 + a * ep2)+o2)
      A = D1*b1*e2H*psi2*a2*b2 - D2*b2*e1H*psi1*a1*b1
      B = D1*b1*e2*a2*D2 + (m1 + a * ep1)*D1*e2H*psi2*a2*b2 - D2*b2*e1*a1*D1 - (m2 + a * ep2)*D2*e1H*psi1*a1*b1
      C = (m1 + a * ep1)*e2*a2*D1*D2 - (m2 + a * ep2)*e1*a1*D2*D1
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
      A1 = (b1 * H) / ((m1 + a * ep1) + o1)
      A2 = (b2 * H) / ((m2 + a * ep2) + o2)
      B1 = h1 * o1 * A1 - c1 * b1 * H
      B2 = h2 * o2 * A2 - c2 * b2 * H
      D1 = (1 + psi1 * A1) * a1
      D2 = (1 + psi2 * A2) * a2
      S = ((b1 * H + (m1 + a * ep1)) * ((m1 + a * ep1) + o1))/ (e1 * a1 * ((m1 + a * ep1) + o1) + e1H * psi1 * a1 * b1 * H)
      #S2 = ((b2 * H + (m2 + a * ep2)) * ((m2 + a * ep2) + o2))/ (e2 * a2 * ((m2 + a * ep2) + o2) + e2H * psi2 * a2 * b2 * H)
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
      S = (m1 + a * ep1)/(e1*a1)
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
      S = (m2 + a * ep2)/(e2*a2)
      P2 = (r/a2)*(1-(S/K))
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P1H = function(parms){
    with(parms,{
      P1 = (d*(o1+(m1 + a * ep1)))/(b1*h1*o1 - c1*b1*(o1+(m1 + a * ep1)))
      A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+(m1 + a * ep1))^2) )
      B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+(m1 + a * ep1)) - (e1*psi1*a1^2*b1*P1) / (r*(o1+(m1 + a * ep1))) )*K - b1
      C = (e1*a1 - (e1*a1^2*P1)/r)*K - (m1 + a * ep1)
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
      S = ((b1 * H + (m1 + a * ep1)) * ((m1 + a * ep1) + o1))/ (e1 * a1 * ((m1 + a * ep1) + o1) + e1H * psi1 * a1 * b1 * H)
      P1H = P1 * (b1 * H) / ((m1 + a * ep1) + o1)
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
      P2 = (d*(o2+(m2 + a * ep2)))/(b2*h2*o2 - c2*b2*(o2+(m2 + a * ep2)))
      A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+(m2 + a * ep2))^2) )
      B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+(m2 + a * ep2)) - (e2*psi2*a2^2*b2*P2) / (r*(o2+(m2 + a * ep2))) )*K - b2
      C = (e2*a2 - (e2*a2^2*P2)/r)*K - (m2 + a * ep2)
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
      S = ((b2 * H + (m2 + a * ep2)) * ((m2 + a * ep2) + o2))/ (e2 * a2 * ((m2 + a * ep2) + o2) + e2H * psi2 * a2 * b2 * H)
      P2H = P2 * (b2 * H) / ((m2 + a * ep2) + o2)
      P1H = 0
      P1 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
    
  }
}