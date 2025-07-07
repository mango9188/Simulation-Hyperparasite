#Stability analisis

####Monoculture equilibrium----
parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = pop_size[length(times), -1]
parms[c('H', 'P1H', 'P1', 'S')] = pop_size[length(times), -1]
parms = as.list(parms)

##for single strain (P1+P1H)
with(parms, {
  J11 = -c1*b1*P1*H-d
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
  
  Jacobian = matrix(data = 
             c(J11, J12, J13, J14,
               J21, J22, J23, J24,
               J31, J32, J33, J34,
               J41, J42, J43, J44), 
    byrow = T, nrow = 4, ncol = 4)
  return(max(Re(eigen(Jacobian)$value)))
  })

##for single strain (P2+P2H)
with(parms, {
  J11 = -c2*b2*P2*H-d
  J12 = h2*o2 
  J13 = -c2*b2*H
  J14 = 0
  J21 = b2*P2
  J22 = -(m2+o2)
  J23 = b2*H
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

##for full model
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
  
  return(max(Re(eigen(Jacobian)$value)))
})

##for P1
with(parms, {
  J11 = e1*a1*S-m1
  J12 = e1*a1*P1
  J21 = -a1*S
  J22 = r*(1-S/K) - a1*P1 - S*r/K
  
  Jacobian = matrix(data = 
                      c(J11, J12, 
                        J21, J22
                        ), 
                    byrow = T, nrow = 2, ncol = 2)
  return(max(Re(eigen(Jacobian)$value)))
})


##for P2
with(parms, {
  J11 = e2*a2*S-m2
  J12 = e2*a2*P2
  J21 = -a2*S
  J22 = r*(1-S/K) - a2*P2 - S*r/K
  
  Jacobian = matrix(data = 
                      c(J11, J12, 
                        J21, J22
                      ), 
                    byrow = T, nrow = 2, ncol = 2)
  return(max(Re(eigen(Jacobian)$value)))
})

##for P1+P2
with(parms, {
  J11 = e1*a1*S-m1
  J12 = 0
  J13 = e1*a1*P1
  J21 = 0
  J22 = e2*a2*S-m2
  J23 = e2*a2*P2
  J31 = -a1*S
  J32 = -a2*S
  J33 = r*(1-S/K) - a1*P1 - a2*P2 - S*r/K
  
  Jacobian = matrix(data = 
                      c(J11, J12, J13,
                        J21, J22, J23,
                        J31, J32, J33
                      ), 
                    byrow = T, nrow = 3, ncol = 3)
  return(max(Re(eigen(Jacobian)$value)))
})


####Stability analysis for for-loop-----
parms <- list(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

for(i in 1:dim(comp_out)[1]){
  parms["a1"] = comp_out[i, "a1"]
  parms["b1"] = comp_out[i, "b1"]
  parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = comp_out[i, c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')]
  
  ifelse(with(parms, {
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
      }) < 0, 
      comp_out[i ,"Stability"] <- "Stable", 
      comp_out[i ,"Stability"] <- "Unstable")
}
View(comp_out)
filter(comp_out, Stability == "Stable")

#long-term average not always equals to equilibrium point!
#It is not appropriate to plug in these value to test stability