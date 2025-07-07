
parms = as.list(parms)
parms <- list(#H = 0.133791930, #0.133791930
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.01, m2 = 0.05, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#Calculate the monoculture equilibrium point
parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = 
  with(parms, {
  D1 = (m1+o1)
  D2 = (m2+o2)
  A = D1*b1*e2H*psi2*a2*b2 - D2*b2*e1H*psi1*a1*b1
  B = D1*b1*e2*a2*D2 + m1*D1*e2H*psi2*a2*b2 - D2*b2*e1*a1*D1 - m2*D2*e1H*psi1*a1*b1
  C = m1*e2*a2*D1*D2 - m2*e1*a1*D2*D1
  H1 = (-B + sqrt(B^2-4*A*C)) / (2*A)
  H2 = (-B - sqrt(B^2-4*A*C)) / (2*A)
  H =
    ifelse(H1 > 0, 
           ifelse(H2 > 0,
                  ifelse(H1 > H2, H2, H1),
                  H2),
           H2)
  
  
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
  
  return(c(H, P1H, P2H, P1, P2, S))
}) #Coexist equilibrium

#Stability analysis, see if Coexistence equilibrium is stable.
if (with(parms, {
  J11 = -c1*b1*P1-c2*b2*P2-d
  J12 = h1 * o1
  J13 = h2 * o2
  J14 = -c1*b1*H
  J15 = -c2*b2*H
  J16 = 0
  J21 = b1*P1
  J22 = -(m1+o1)
  J23 = 0
  J24 = b1 * H
  J25 = 0
  J26 = 0
  J31 = b2*P2
  J32 = 0
  J33 = -(m2+o2)
  J34 = 0
  J35 = b2*H
  J36 = 0
  J41 = -b1*P1
  J42 = e1H * psi1 * a1 * S
  J43 = 0
  J44 = e1*a1*S - b1*H - m1
  J45 = 0
  J46 = e1*a1*P1 + e1H*psi1*a1*P1H
  J51 = -b2*P2
  J52 = 0
  J53 = e2H * psi2 * a2 * S
  J54 = 0
  J55 = e2*a2*S - b2*H - m2
  J56 = e2*a2*P2 + e2H*psi2*a2*P2H
  J61 = 0
  J62 = -psi1*a1*S
  J63 = -psi2*a2*S
  J64 = -a1*S
  J65 = -a2*S
  J66 = r*(1-S/K) - a1*P1 - psi1*a1*P1H - a2*P2 - psi2*a2*P2H - S*r/K
  
  Jacobian = matrix(data = 
                      c(J11, J12, J13, J14, J15, J16,
                        J21, J22, J23, J24, J25, J26,
                        J31, J32, J33, J34, J35, J36,
                        J41, J42, J43, J44, J45, J46,
                        J51, J52, J53, J54, J55, J56,
                        J61, J62, J63, J64, J65, J66), 
                    byrow = T, nrow = 6, ncol = 6)
  
  return(max(Re(eigen(Jacobian)$value)))
}) > 0 ){
  #if unstable (> 0), try E_Pi
  parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = 
  with(parms, {
    H = 0
    P1H = 0
    P2H = 0
    P1 = (r/a1)*(1-(S/K))
    P2 = 0
    S = m1/(e1*a1)
    return(c(H, P1H, P2H, P1, P2, S))
  }) #E_P1
  if (with(parms, {
    J11 = -c1*b1*P1-c2*b2*P2-d
    J12 = h1 * o1
    J13 = h2 * o2
    J14 = -c1*b1*H
    J15 = -c2*b2*H
    J16 = 0
    J21 = b1*P1
    J22 = -(m1+o1)
    J23 = 0
    J24 = b1 * H
    J25 = 0
    J26 = 0
    J31 = b2*P2
    J32 = 0
    J33 = -(m2+o2)
    J34 = 0
    J35 = b2*H
    J36 = 0
    J41 = -b1*P1
    J42 = e1H * psi1 * a1 * S
    J43 = 0
    J44 = e1*a1*S - b1*H - m1
    J45 = 0
    J46 = e1*a1*P1 + e1H*psi1*a1*P1H
    J51 = -b2*P2
    J52 = 0
    J53 = e2H * psi2 * a2 * S
    J54 = 0
    J55 = e2*a2*S - b2*H - m2
    J56 = e2*a2*P2 + e2H*psi2*a2*P2H
    J61 = 0
    J62 = -psi1*a1*S
    J63 = -psi2*a2*S
    J64 = -a1*S
    J65 = -a2*S
    J66 = r*(1-S/K) - a1*P1 - psi1*a1*P1H - a2*P2 - psi2*a2*P2H - S*r/K
    
    Jacobian = matrix(data = 
                        c(J11, J12, J13, J14, J15, J16,
                          J21, J22, J23, J24, J25, J26,
                          J31, J32, J33, J34, J35, J36,
                          J41, J42, J43, J44, J45, J46,
                          J51, J52, J53, J54, J55, J56,
                          J61, J62, J63, J64, J65, J66), 
                      byrow = T, nrow = 6, ncol = 6)
    
    return(max(Re(eigen(Jacobian)$value)))
  }) > 0 ){
    #if unstable (> 0), try E_P2
    parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = 
      with(parms, {
        H = 0
        P1H = 0
        P2H = 0
        P1 = 0
        P2 = (r/a2)*(1-(S/K))
        S = m2/(e2*a2)
        return(c(H, P1H, P2H, P1, P2, S))
      })#E_P2
  }else{
    print("Both are unstable, H should persist.")
    
    parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = 
      with(parms,{
        P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
        P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
        #When mono-cul, system follow S* and H* rule,
        #i.e., the one with small S* and large H*  will win.
        
        #calculate H1*
        A1 = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
        B1 = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
        C1 = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
        
        H1_1 = (-B1-sqrt(B1^2 - 4*A1*C1)) / (2*A1)
        H1_2 = (-B1+sqrt(B1^2 - 4*A1*C1)) / (2*A1)
        H1 = ifelse(H1_1 > 0, H1_1, H1_2)
        
        #calculate H2*
        A2 = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
        B2 = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
        C2 = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
        H2_1 = (-B2-sqrt(B2^2 - 4*A2*C2)) / (2*A2)
        H2_2 = (-B2+sqrt(B2^2 - 4*A2*C2)) / (2*A2)
        H2 = ifelse(H2_1 > 0, H2_1, H2_2)
        
        #compare H1* and H2*
        if (H1 > H2){
          H = H1
          S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
          P1H = P1 * (b1 * H) / (m1 + o1)
          return(c(H, P1H, P1, S))
        }else if(H1 < H2){
          H = H2
          S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
          P2H = P2 * (b2 * H) / (m2 + o2)
          
          return(c(H, P2H, P2, S))
        }else{
          print("You fucked up")
        }
      })
  }
}

parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')]



if (H1>H2){
  print("H1>H2")
}else if (H1<H2){
  print("H1 < H2")
}else{
  print("H1 = H2")
}


