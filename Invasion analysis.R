library(tidyverse)

#For H and Pi/H's to invade E_SPi----
parms <- list(
  r = 1, K = 10,
  a1 = 0.25, psi1 = 1, e1 = 0.5,
  b1 = 0.2, m1 = 0.05, e1H = 0.5,
  o1 = 0.8, h1 = 1,  c1 = 0.9, d = 0.03, DL = 0)
## Leslie matrix for H and PiH's invasion growth rate----
###For P1H+H invade P1----
Leslie <- matrix(
  data = with(parms, {
    P1 = (r/a1)*(1-(m1/(e1*a1*K)))
    M = 
      c(-(c1 * b1 * P1 + d), h1 * o1, b1 * P1, -(m1 + o1))
    return(M)
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

###For P2H+H invade P2----
Leslie <- matrix(
  data = with(parms, {
    P2 = (r/a2)*(1-(m2/(e2*a2*K)))
    M = 
      c(-(c2 * b2 * P2 + d), h2 * o2, b2 * P2, -(m2 + o2))
    return(M)
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

EIGEN <- eigen(Leslie)
### Dominant eigenvalue -> invasion growth rate
EIGEN$values

##Analytic sol of eigenvalue----
with(parms, {
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  A = 1
  # B = (c1*b1*P1 + d + m1*o1)
  B = (c1*b1*P1 + d + m1+o1)
  # C = (c1*b1*P1 + d) * (m1*o1) - h1*o1*b1*P1
  C = (c1*b1*P1 + d) * (m1+o1) - h1*o1*b1*P1
  eigen1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
  eigen2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
  return(c(eigen1, eigen2))
})

##Approach 1----
with(parms,{
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  # A = b1*P1 -c1*b1*P1 - d + h1*o1 - m1*o1
  A = -b1*P1 +c1*b1*P1 + d + h1*o1 - m1-o1
  # B = -(h1*o1 - m1*o1) - (c1*b1*P1 + d) - h1 * o1
  B = -(h1*o1 - m1-o1) - (c1*b1*P1 + d) - h1 * o1
  C = h1 * o1
  phiH1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
  phiH2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
  #note that phiH must be greater than zero and smaller than 1.
  IGR1 = 
    ifelse(phiH1 > 0, 
           ifelse(phiH1 < 1, 
                  # (h1*o1 - m1*o1) * (1-phiH1) + (b1*P1- c1*b1*P1 - d) * phiH1,
                  (h1*o1 - m1-o1) * (1-phiH1) + (b1*P1- c1*b1*P1 - d) * phiH1,
                  "Not meaningful"),
           "Not meaningful")
  IGR2 =
    ifelse(phiH2 > 0,
          ifelse(phiH2 < 1,
                # (h1*o1 - m1*o1) * (1-phiH2) + (b1*P1- c1*b1*P1 - d) * phiH2,
                  (h1*o1 - m1-o1) * (1-phiH2) + (b1*P1- c1*b1*P1 - d) * phiH2,
                  "Not meaningful"),
          "Not meaningful")
  print(c(IGR1, IGR2))
  #print(c(phiH1, phiH2))
})

##Approach2----
with(parms,{
  P1 = (r/a1)*(1-(m1/(e1*a1*K)))
  A = h1*o1
  B = m1+o1 - c1*b1*P1 - d
  C = -b1*P1
  Z1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
  Z2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
  #note that phiH must be greater than zero
  #it is possible that both roots are positive.
  IGR1 = 
    ifelse(Z1 > 0,
           (h1*o1*Z1 - c1*b1*P1 - d),
           "Not meaningful")
  IGR2 =
    ifelse(Z2 > 0,
           (h1*o1*Z2 - c1*b1*P1 - d),
           "Not meaningful")
  print(c(IGR1, IGR2))
  #print(c(Z1, Z2))
})



#For P_i and P_{i/H} to invade E_{S P_j P_j/H}----
## Leslie matrix for Pi and Pi/H's invasion growth rate----
parms = list(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.04, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#For P1+P1H invade P2+P2H equilibrium
Leslie <- matrix(
  data = 
    with(parms, {
      P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
      A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
      B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
      C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
      H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
      H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
      H = 
        ifelse(H1 > 0, H1, 
               ifelse(H2 > 0, H2, break))
      S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
      
    return(c(-(m1+o1), b1*H,
           e1H*psi1*a1*S, e1*a1*S-b1*H-m1))
  }),
  nrow = 2,
  ncol = 2,
  byrow = T)

#For P2+P2H invade P1+P1H equilibrium
Leslie <- matrix(
  data = 
    with(parms, {
      P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
      A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
      B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
      C = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
      H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
      H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
      H = 
        ifelse(H1 > 0, H1, 
                 ifelse(H2 > 0, H2, break))
      S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      
      return(c(-(m2+o2), b2*H,
               e2H*psi2*a2*S, e2*a2*S-b2*H-m2))
    }),
  nrow = 2,
  ncol = 2,
  byrow = T)

EIGEN <- eigen(Leslie)
### Dominant eigenvalue (The largest and positive one) -> invasion growth rate
max(EIGEN$values)

## Analytic sol of eigenvalue----
#For P1+P1H invade P2+P2H equilibrium
with(parms,{
  P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
  A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
  B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
  C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
  H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
  H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
  H = ifelse(H1 > 0, H1, H2)
  S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  A = 1
  B = m1+o1-(e1*a1*S-b1*H-m1)
  C = -b1*H*e1H*psi1*a1*S-(m1+o1)*(e1*a1*S-b1*H-m1)
  eigen1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
  eigen2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
  return(c(eigen1, eigen2))
})

##Quasi equilibrium Approach 1.1 (For P1+P1H invade P2+P2H equilibrium)----
with(parms,{
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
  #A = -e1*a1*S + e1H*psi1*a1*S - (m1+o1) + m1
  A = -e1*a1*S + e1H*psi1*a1*S + o1
  #B = e1*a1*S - b1*H - m1 - e1H*psi1*a1*S + (m1+o1)
  B = e1*a1*S - 2*e1H*psi1*a1*S - b1*H + o1
  #C = e1H*psi1*a1*S
  C = e1H*psi1*a1*S
  
  phiPa = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
  phiPb = (-B - sqrt(B^2 - 4*A*C)) / (2*A)
  #note that phiP must be greater than zero and smaller than 1.
  phiP = 
    ifelse(phiPa > 0 && phiPa < 1, phiPa,
           ifelse(phiPb > 0 && phiPb < 1, phiPb, "Not meaningful"))
  print(c(phiPa, phiPb))
  IGR = e1*a1*S*phiP + e1H*psi1*a1*S*(1-phiP) - (m1+o1)*(1-phiP) - m1*phiP
  print(IGR)
})

##Quasi equilibrium Approach 1.2 (For P1+P1H invade P2+P2H equilibrium)----
with(parms,{
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

  PHI2 = b1*H

  PHI1 = b1*H + m1 - (m1+o1) - e1*a1*S

  PHI0 = -e1H*psi1*a1*S
  
  phiPa = (-PHI1 + sqrt(PHI1^2 - 4*PHI2*PHI0)) / (2*PHI2)
  phiPb = (-PHI1 - sqrt(PHI1^2 - 4*PHI2*PHI0)) / (2*PHI2)
  #note that phiP must be greater than zero and smaller than 1.
  phiP = 
    ifelse(phiPa > 0, phiPa,
           ifelse(phiPb > 0, phiPb, "Not meaningful"))
  print(c(phiPa, phiPb))
  IGR = b1*H*phiP - (m1+o1)
  print(IGR)
})


##Quasi equilibrium Approach 2 (shorcut) (For P1+P1H invade P2+P2H equilibrium)----
with(parms, {
  P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
  A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
  B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
  C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
  H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
  H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
  H = ifelse(H1 > 0, H1, H2)
  S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  IGR = e1*a1*S - b1*H - m1 + (b1*e1H*a1*S*H) / (m1+o1)
  return(IGR)
})
#This will not be the actual IGR

#Data analysis----
parms = list(
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

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

comp_out$IGR1 = mapply(f_IGR1, comp_out$m1, comp_out$m2)
comp_out$IGR2 = mapply(f_IGR2, comp_out$m1, comp_out$m2)
comp_out$IGRH1 = mapply(f_IGRH1, comp_out$m1)
comp_out$IGRH2 = mapply(f_IGRH2, comp_out$m2)
