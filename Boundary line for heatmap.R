####This is the script for the boundary line.
#####Remember to define the "parms" before calling the data!
#comp_out = readRDS("Pre5A1")
#comp_out = filter(comp_out, d == 0.02)
#comp_out = comp_out[,-3]
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
comp_out$C_State = mapply(f_C_State, comp_out$m1, comp_out$m2)
comp_out$C_State_V = ifelse(comp_out$C_State == "Stable", 1, 
                            ifelse(comp_out$C_State == "Unstable", 0, -1))

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
comp_out$E2H_State = mapply(f_E2H_State, comp_out$m1, comp_out$m2)
comp_out$E2H_State_V = ifelse(comp_out$E2H_State == "Stable", 1, 
                            ifelse(comp_out$E2H_State == "Unstable", 0, -1))
View(comp_out)
