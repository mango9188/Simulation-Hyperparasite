#Find the minimun H under different parameter sets.

parms = list(
  r = 1, K = 10,
  a1 = 0.4, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#parms = list(r = 1, K = 10, a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5, b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5, o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

library(rootSolve)

####Get the cross point----
CP =
  function(H){
    with(parms, {
      S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
      S1-S2
    })
  }

#Min.H = uniroot(CP, lower = 0, upper = 5)
#Min.H$root
#curve(CP(x), from = 0, to = 10)
#abline(h = 0, lty = 2)

comp_out = expand.grid(a1 = seq(0.05, 1, by = 0.05), b1 = seq(0.05, 1, by = 0.05))

####Input each parameter sets----
comp_out = as.data.frame(cbind(comp_out, 
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 3)))
names(comp_out) = c("a1", "b1", "S1_H0", "S2_H0", "Min.H")
#The value of Si here is the value when H = 0.

comp_out[, "S2_H0"] = 
  with(parms, {
    S2 = m2/(e2*a2)
    return(S2)
  })

###Create a parameter space with analytic way----
for(i in 1:dim(comp_out)[1]){
  temp_parms = parms
  temp_parms["a1"] = comp_out[i, "a1"]
  temp_parms["b1"] = comp_out[i, "b1"]
  
  comp_out[i, "S1_H0"] = 
    with(temp_parms, {
      S1 = m1/(e1*a1)
      return(S1)
    })
  
  CP =
    function(H){
      with(temp_parms, {
        S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
        S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
        S1-S2
      })
    }
  comp_out[i, "Min.H"] = 
    tryCatch(
      {
        uniroot(CP, lower = 0, upper = 1, tol = 1e-16)$root
      }, error = function(e){
        NA #S1 and S2 do not cross over, system will predictably follow exploitative competition that only consist S, P1 and P2 no matter how large the H is.
      })
  
  if(is.na(comp_out[i, "Min.H"]) == FALSE){
      temp_parms["H"] = comp_out[i,"Min.H"]
      comp_out[i, c('P1H', 'P2H', 'P1', 'P2', 'S')] = 
        with(temp_parms, {
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
      return(c(P1H, P2H, P1, P2, S1))
    })
      temp_parms[c('P1H', 'P2H', 'P1', 'P2', 'S')] = 
        comp_out[i, c('P1H', 'P2H', 'P1', 'P2', 'S')]
      ifelse(with(temp_parms, {
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
  
  
  #If there's Min.H, calculate the coexistence equilibrium and its stability

}


comp_out$Outcome.ZeroH = 
  ifelse(comp_out[, 3] < comp_out[, 4], "P1(+H)", 
         ifelse(comp_out[, 3] > comp_out[, 4], "P2(+H)", "Coexist"))
View(comp_out)


#Find the corresponding S under different parameter sets----


parms = list(
  r = 1, K = 10,
  a1 = 0.4, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#parms = list(r = 1, K = 10, a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5, b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5, o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

library(rootSolve)

##Get the cross point
CP =
  function(S){
    with(parms, {
      H1 = (m1 - e1 * a1 * S)/((e1H * psi1 * a1 * b1 * S)/ (m1 + o1) - b1)
      H2 = (m2 - e2 * a2 * S)/((e2H * psi2 * a2 * b2 * S)/ (m2 + o2) - b2)
      H1-H2
    })
  }

#Min.H = uniroot(CP, lower = 0, upper = 5)
#Min.H$root
#curve(CP(x), from = 0, to = 10)
#abline(h = 0, lty = 2)

comp_outS = expand.grid(a1 = seq(0.05, 1, by = 0.05), b1 = seq(0.05, 1, by = 0.05))

####Input each parameter sets----
comp_outS = as.data.frame(cbind(comp_outS, 
                               matrix(0, 
                                      nrow = dim(comp_outS)[1],
                                      ###dim(data)[1] is the number of row of data; [2] is col
                                      ncol = 3)))
names(comp_outS) = c("a1", "b1", "H1", "H2", "S")

###Create a parameter space with analytic way----
for(i in 1:dim(comp_outS)[1]){
  temp_parms = parms
  temp_parms["a1"] = comp_outS[i, "a1"]
  temp_parms["b1"] = comp_outS[i, "b1"]
  
  CP =
    function(S){
      with(temp_parms, {
        H1 = (m1 - e1 * a1 * S)/((e1H * psi1 * a1 * b1 * S)/ (m1 + o1) - b1)
        H2 = (m2 - e2 * a2 * S)/((e2H * psi2 * a2 * b2 * S)/ (m2 + o2) - b2)
        H1-H2
      })
    }
  
  comp_outS[i, "S"] = 
    tryCatch(
      {
        uniroot(CP, lower = 0, upper = 1, tol = 1e-16)$root
      }, error = function(e){
        NA #S1 and S2 do not cross over, system will predictably follow exploitative competition that only consist S, P1 and P2 no matter how large the H is.
      })
  if (is.na(comp_outS[i, "S"]) == FALSE){
    temp_parms$S = comp_outS[i, "S"]
    comp_outS[i, c("H1", "H2")] =
      with(temp_parms, {
        H1 = (m1 - e1 * a1 * S)/((e1H * psi1 * a1 * b1 * S)/ (m1 + o1) - b1)
        H2 = (m2 - e2 * a2 * S)/((e2H * psi2 * a2 * b2 * S)/ (m2 + o2) - b2)
      })
  }else{
    next
  }
  
}

View(comp_outS)
 