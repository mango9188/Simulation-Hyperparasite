#Stability analisis

####Monoculture equilibrium----
#1 for full model simulation
parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = pop_size[length(times), -1]
parms = as.list(parms)
##for P1 wins
with(parms, {
  J11 = -c1*b1*P1*H-d
  J12 = h1*o1 
  J13 = -c1*b1*H
  J14 = 0
  J21 = b1*P1
  J22 = -(m1+o1)
  J23 = b1*P1
  J24 = 0
  J31 = -b1*P1 
  J32 = e1H*psi1*a1*S
  J33 = e1*a1*S-b1*H-m1
  J34 = e1*a1*P1+e1H*psi1*a1*P1H
  J41 = 0
  J42 = -psi1*a1*S
  J43 = -a1*S
  J44 = r*(1-S/K) - a1*P1 - psi1*a1*P1H - S*r/K
  
  Jacobian = matrix(data = 
             c(J11, J12, J13, J14,
               J21, J22, J23, J24,
               J31, J32, J33, J34,
               J41, J42, J43, J44), 
    byrow = T, nrow = 4, ncol = 4)
  return(Re(eigen(Jacobian)$value))
  })

##for P2 wins
with(parms, {
  J11 = -c2*b2*P2*H-d
  J12 = h2*o2 
  J13 = -c2*b2*H
  J14 = 0
  J21 = b2*P2
  J22 = -(m2+o2)
  J23 = b2*P2
  J24 = 0
  J31 = -b2*P2
  J32 = e2H*psi2*a2*S
  J33 = e2*a2*S-b2*H-m2
  J34 = e2*a2*P2+e2H*psi2*a2*P2H
  J41 = 0
  J42 = -psi2*a2*S
  J43 = -a2*S
  J44 = r*(1-S/K) - a2*P2 - psi2*a2*P2H - S*r/K
  
  Jacobian = matrix(data = 
                      c(J11, J12, J13, J14,
                        J21, J22, J23, J24,
                        J31, J32, J33, J34,
                        J41, J42, J43, J44), 
                    byrow = T, nrow = 4, ncol = 4)
  return(Re(eigen(Jacobian)$value))
})



####Coexistence equilibrium----
with(parms, {
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
  
  return(Re(eigen(Jacobian)$value))
})


