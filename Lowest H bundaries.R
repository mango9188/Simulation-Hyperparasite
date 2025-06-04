#Find the minimun H under different parameter sets.

parms = list(
  r = 1, K = 10,
  a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
  b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.06, e1H = 0.5, e2H = 0.5,
  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

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

comp_out = expand.grid(a1 = seq(0.05, 0.5, by = 0.05), b1 = seq(0.05, 0.5, by = 0.05))

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
}

comp_out$Outcome.ZeroH = 
  ifelse(comp_out[, 3] < comp_out[, 4], "P1(+H)", 
         ifelse(comp_out[, 3] > comp_out[, 4], "P2(+H)", "Coexist"))

View(comp_out)

