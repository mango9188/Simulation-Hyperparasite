###Verify low H model with analytic solution----
comp_out = readRDS("Pre3Sim1")
comp_out2 = readRDS("Pre1Sim3NC3")
View(comp_out2)

comp_out = 
mutate(comp_out,
       Ratio = P1H/P1,
       True.Ratio = (b1 * H) / (m1 + o1))  #P1H*P1 = (b1 * H) / (m1 + o1)


parms <- list(H = 0.133791930, #0.133791930
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.0536, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms = list(
  H = 0.122431816,
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms = list(
  H = 1.141915e-01, #6.657334e-02
  r = 1, K = 10,
  a1 = 0.5, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.5, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

###Mono-culture equilibrium----
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

S1 =
  with(parms, {
    ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
  })
  

S2 = 
  with(parms, {
    ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  })

P1 = 
  with(parms, {
    r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
  })

P1 = 
  with(parms, {
    (d*(o1+m1))/(b1*h1*o1-c1*b1*(o1+m1))
  })
  
P2 = 
  with(parms, {
    r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
  })

P2 = 
  with(parms, {
    (d*(o2+m2))/(b2*h2*o2-c2*b2*(o2+m2))
  })

P1H = 
  with(parms, {
    P1 * (b1 * H) / (m1 + o1)
  })

P2H =
  with(parms, {
    P2 * (b2 * H) / (m2 + o2)
  })

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

####Coexist equilibrium----

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
