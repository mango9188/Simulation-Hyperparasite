###Verify low H model with analytic solution----
comp_out = readRDS("Pre3Sim1")
comp_out2 = readRDS("Pre1Sim3NC3")
View(comp_out2)

comp_out = 
mutate(comp_out,
       Ratio = P1H/P1,
       True.Ratio = (b1 * H) / (m1 + o1))  #P1H*P1 = (b1 * H) / (m1 + o1)


parms <- list(#H = 0.133791930, #0.133791930
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms = parms %>% as.list()
#Mono-culture equilibrium (need H to solve)----
Equilibrium = 
  with(parms, {
    S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
    S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
    P1 = r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
    P2 = r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
    P1H = P1 * (b1 * H) / (m1 + o1)
    P2H = P2 * (b2 * H) / (m2 + o2)
    return(c(P1H, P2H, P1, P2, S1, S2))
  })

parms[c('P1H', 'P2H', 'P1', 'P2', 'S1', 'S2')] = Equilibrium

f_H = function(H, parms){
  with(parms,{
    S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
    P1 = r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
    P1H = P1 * (b1 * H) / (m1 + o1)
    approx_H = (h1 * o1 * P1H) / (c1 * b1 * P1 + d)
    return(approx_H)
  })
}

H.value = 
  uniroot(f_H, interval = c(0, 1), tol = 1e-16, parms = parms)$root

f <- function(H, parms) {
  with(parms, {
    S1_num <- (b1 * H + m1) * (m1 + o1)
    S1_den <- e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H
    S1 <- S1_num / S1_den
    
    P1 <- r * (1 - S1 / K) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
    
    lhs <- P1 * (h1 * o1 * b1 / (m1 + o1) - c1 * b1)
    result <- lhs - d
    return(result)
  })
}
root <- uniroot(f, interval = c(0, 1), parms = parms)$root



#Mono-culture equilibrium (Full analytic)----
E1 = 
  with(parms, {
    P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
    A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
    B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
    C = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
    H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
    H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
    H = ifelse(H1 > 0, H1, H2)
    S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
    P1a = r * (1 - (S/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
    P1H = P1 * (b1 * H) / (m1 + o1)
    return(setNames(
      c(H1, H2, P1H, P1, P1a, S),
      c("H1", "H2", "P1H", "P1", "P1a", "S")
      #c(H, P1H, P1, S),
      #c("H", "P1H", "P1", "S")
    )
      )
  })

E1

E2 = 
  with(parms, {
    P2 = (d*(o2+m2))/(b2*h2*o2 - c2*b2*(o2+m2))
    A = -( (e2H*psi2^2*a2^2*b2^2*P2*K) / (r*(o2+m2)^2) )
    B = ( (e2H*psi2*a2*b2 - e2H*psi2*a2^2*b2*P2) / (o2+m2) - (e2*psi2*a2^2*b2*P2) / (r*(o2+m2)) )*K - b2
    C = (e2*a2 - (e2*a2^2*P2)/r)*K - m2
    H1 = (-B-sqrt(B^2 - 4*A*C)) / (2*A)
    H2 = (-B+sqrt(B^2 - 4*A*C)) / (2*A)
    H = ifelse(H1 > 0, H1, H2)
    S = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
    P2a = r * (1 - (S/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
    P2H = P2 * (b2 * H) / (m2 + o2)
    
    return(setNames(
      c(H1, H2, P2H, P2, P2a, S),
      c("H1", "H2", "P2H", "P2", "P2a", "S")))
  })
E2

E = 
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
      return("Coexist")
    }
  })

#Coexist equilibrium (need H to solve)----
Equilibrium =
  with(parms, {
    A1 = (b1 * H) / (m1 + o1)
    A2 = (b2 * H) / (m2 + o2)
    B1 = h1 * o1 * A1 - c1 * b1 * H
    B2 = h2 * o2 * A2 - c2 * b2 * H
    D1 = (1 + psi1 * A1) * a1
    D2 = (1 + psi2 * A2) * a2
    S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
    S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
    P1 = (d * H) / B1 - (B2 / B1) * (((B1 * r * (1 - S1 / K)) - D1 * d * H) / (D2 * B1 - D1 * B2))
    P2 = (B1 * r * (1 - S1 / K) - D1 * d * H) / (D2 * B1 - D1 * B2)
    P1H = A1 * P1
    P2H = A2 * P2
    return(c(P1H, P2H, P1, P2, S1, S2))
  })

names(Equilibrium) = c('P1H', 'P2H', 'P1', 'P2', 'S1', 'S2')
Equilibrium
parms[c('P1H', 'P2H', 'P1', 'P2', 'S1', 'S2')] = Equilibrium

with(parms, {
  round ((h1 * o1 * A1 - c1 * b1 *H) * P1, 7) == round (d * H - (h1 * o1 * A1 - c1 * b1 *H) * P2, 7)
})

Equilibrium["P1"]


#Coexist equilibrium (Full analytic)----
E =
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
  })
names(E) = c("H", "P1H", "P2H", "P1", "P2", "S")
E


#Mono-culture equilibrium (with no H)----
parms_E_P1 = parms
parms_E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')]  =
  with(parms, {
    H = 0
    P1H = 0
    P2H = 0
    P2 = 0
    S = m1/(e1*a1)
    P1 = (r/a1)*(1-(S/K))
    return(c(H, P1H, P2H, P1, P2, S))
  })
E_P1

parms_E_P2 = parms
parms_E_P2[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] =
  with(parms, {
    H = 0
    P1H = 0
    P2H = 0
    P1 = 0
    S = m2/(e2*a2)
    P2 = (r/a2)*(1-(S/K))
    return(c(H, P1H, P2H, P1, P2, S))
  })
E_P2
