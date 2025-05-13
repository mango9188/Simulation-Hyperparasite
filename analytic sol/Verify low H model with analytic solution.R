###Verify low H model with analytic solution----
comp_out = readRDS("Pre3Sim1")
comp_out2 = readRDS("Pre1Sim3NC3")

comp_out = 
mutate(comp_out,
       Ratio = P1H/P1,
       True.Ratio = (b1 * H) / (m1 + o1))  #P1H*P1 = (b1 * H) / (m1 + o1)


parms <- list(H = 0.0122447840,
           r = 1, K = 10,
           a1 = 0.4, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.25, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms = list(
  H = 0.06657334,
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms = list(
  H = 0.181164775, #6.657334e-02
  r = 1, K = 10,
  a1 = 0.3, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.15, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)


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
  
P2 = 
  with(parms, {
    r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
  })

P1H = 
  with(parms, {
    P1 * (b1 * H) / (m1 + o1)
  })

P2H =
  with(parms, {
    P2 * (b2 * H) / (m2 + o2)
  })
