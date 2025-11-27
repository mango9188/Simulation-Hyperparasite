#Using analytic method to expand parameter space----
#####!!!!!This method is used to draw the bifurcation plot!
##Expanding parameter space----
#### Parameter sets----
parms <- list(r = 1, K = 10,
              a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
              b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
              o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

comp_out = expand.grid(m1 = seq(0.001, 0.1, by = 0.001), 
                       m2 = seq(0.001, 0.1, by = 0.001))
#### Function setting----
{
  Jacobian = function(parms, E) {
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
  
  f_E_C = function(parms){
    with(parms, {
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
  }
  
  f_E_S = function(parms){
    with(parms,{
      H = 0
      P1H = 0
      P2H = 0
      P1 = 0
      P2 = 0
      S = K
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P1 = function(parms){
    with(parms, {
      H = 0
      P1H = 0
      P2H = 0
      S = m1/(e1*a1)
      P1 = (r/a1)*(1-(S/K))
      P2 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P2 = function(parms){
    with(parms, {
      H = 0
      P1H = 0
      P2H = 0
      P1 = 0
      S = m2/(e2*a2)
      P2 = (r/a2)*(1-(S/K))
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P1H = function(parms){
    with(parms,{
      P1 = (d*(o1+m1))/(b1*h1*o1 - c1*b1*(o1+m1))
      A = -( (e1H*psi1^2*a1^2*b1^2*P1*K) / (r*(o1+m1)^2) )
      B = ( (e1H*psi1*a1*b1 - e1H*psi1*a1^2*b1*P1) / (o1+m1) - (e1*psi1*a1^2*b1*P1) / (r*(o1+m1)) )*K - b1
      C = (e1*a1 - (e1*a1^2*P1)/r)*K - m1
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
      S = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      P1H = P1 * (b1 * H) / (m1 + o1)
      P2H = 0
      P2 = 0
      return(setNames(
        c(H, P1H, P2H, P1, P2, S),
        c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')
      )) 
    })
  }
  
  f_E_P2H = function(parms){
    with(parms, {
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
    
  }
}

#### Solve equilibrium across the parameter space ----
comp_list = list()
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  temp_parms = parms
  temp_parms$m1 = comp_out[i, "m1"]
  temp_parms$m2 = comp_out[i, "m2"]
  # temp_parms$psi1 = comp_out[i, "psi1"]
  # temp_parms$psi2 = comp_out[i, "psi2"]
  
  #Calculate each equilibrium point
  E_C = f_E_C(temp_parms)
  E_S = f_E_S(temp_parms)
  E_P1 = f_E_P1(temp_parms)
  E_P2 = f_E_P2(temp_parms)
  E_P1H = f_E_P1H(temp_parms)
  E_P2H = f_E_P2H(temp_parms)
  
  Lambda_E_C = Eigen(Jacobian(temp_parms, E_C))
  Lambda_E_S = Eigen(Jacobian(temp_parms, E_S))
  Lambda_E_P1 = Eigen(Jacobian(temp_parms, E_P1))
  Lambda_E_P2 = Eigen(Jacobian(temp_parms, E_P2))
  Lambda_E_P1H = Eigen(Jacobian(temp_parms, E_P1H))
  Lambda_E_P2H = Eigen(Jacobian(temp_parms, E_P2H))
  
  Num.E = 0
  
  if(!is.na(Lambda_E_C) && is.finite(Lambda_E_C) && Lambda_E_C < 0 && all(E_C[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
        data.frame(m1 = temp_parms$m1,
                   m2 = temp_parms$m2,
                   as.list(E_C),
                   Stable_E = "C",
                   eigen_value = Lambda_E_C)
    Num.E = Num.E + 1
  }
  
  if(!is.na(Lambda_E_S) && is.finite(Lambda_E_S) && Lambda_E_S < 0 && all(E_S[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 as.list(E_S),
                 Stable_E = "S",
                 eigen_value = Lambda_E_S)
    Num.E = Num.E + 1
  }
  
  if(!is.na(Lambda_E_P1) && is.finite(Lambda_E_P1) && Lambda_E_P1 < 0 && all(E_P1[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 as.list(E_P1),
                 Stable_E = "P1",
                 eigen_value = Lambda_E_P1)
    Num.E = Num.E + 1
  }
  
  if(!is.na(Lambda_E_P2) && is.finite(Lambda_E_P2) && Lambda_E_P2 < 0 && all(E_P2[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 as.list(E_P2),
                 Stable_E = "P2",
                 eigen_value = Lambda_E_P2)
    Num.E = Num.E + 1
  }
  
  if(!is.na(Lambda_E_P1H) && is.finite(Lambda_E_P1H) && Lambda_E_P1H < 0 && all(E_P1H[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 as.list(E_P1H),
                 Stable_E = "P1H",
                 eigen_value = Lambda_E_P1H)
    Num.E = Num.E + 1
  }
  
  if(!is.na(Lambda_E_P2H) && is.finite(Lambda_E_P2H) && Lambda_E_P2H < 0 && all(E_P2H[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] >= 0)){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 as.list(E_P2H),
                 Stable_E = "P2H",
                 eigen_value = Lambda_E_P2H)
    Num.E = Num.E + 1
  }
  
  if(Num.E == 0){
    comp_list[[length(comp_list) + 1]] = 
      data.frame(m1 = temp_parms$m1,
                 m2 = temp_parms$m2,
                 H = NA, P1H = NA, P2H = NA, P1 = NA, P2 = NA, S = NA,
                 Stable_E = "U",
                 eigen_value = NA)
  }
}

end_time <- Sys.time()
end_time - start_time

comp_list
comp_out = do.call(rbind, comp_list)

#### Data analysis ----


View(comp_out)
saveRDS(comp_out, "Pre4A1_FullASS")
