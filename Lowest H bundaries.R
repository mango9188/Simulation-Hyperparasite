#Find the minimun H under different parameter sets.

parms = list(
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
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

Jacobian = function(parms) {
  with(parms, {
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

#Min.H = uniroot(CP, lower = 0, upper = 5)
#Min.H$root
#curve(CP(x), from = 0, to = 10)
#abline(h = 0, lty = 2)

comp_out = expand.grid(a1 = seq(0.05, 0.5, by = 0.05), b1 = seq(0.05, 0.5, by = 0.05))
comp_out = expand.grid(a1 = seq(0.35, 0.35), b1 = seq(0.2, 0.2))
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
      J = Eigen(Jacobian(temp_parms))
      ifelse(!is.na(J) < 0 && is.finite(J) && J < 0,
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
 